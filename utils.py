import pandas as pd
import numpy as np

from schema import TRGTDeNovoSchema


def al_in_parents(row: pd.Series) -> bool:
    """check if the child's de novo allele length is equal to
    either of their two parents' diploid allele lengths.

    Args:
        row (pd.Series): pandas Series object.

    Returns:
        bool: simple True/False boolean
    """
    # access the de novo allele length in the child
    child_als = row["child_AL"].split(",")
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]

    # access the two parents' allele lengths
    mom_als, dad_als = row["mother_AL"].split(","), row["father_AL"].split(",")

    # if we're on a male sex chromosome, we only have to check to see
    # if the child's AL is present in a single parent.
    if row["suffix"].startswith("S"):
        if row["trid"].startswith("chrX") and denovo_al in mom_als: return True
        elif row["trid"].startswith("chrY") and denovo_al in dad_als: return True
        else: return False

    # otherwise, we'll check to see if the de novo AL is in one parent, and
    # the non de novo AL is in the other.
    else:
        non_denovo_idx = 1 - denovo_idx
        non_denovo_al = child_als[non_denovo_idx]
        if (denovo_al in mom_als) and (non_denovo_al in dad_als):
            return True
        elif (denovo_al in dad_als) and (non_denovo_al in mom_als):
            return True
        else: return False


def add_likely_de_novo_size(row: pd.Series, use_parsimony: bool = False,) -> int:
    """determine the most likely size of the de novo expansion or contraction
    in the child (measured in base pairs).

    Args:
        row (pd.Series): pandas Series.
        use_parsimony (bool, optional): whether to use a simple parsimony
            based approach to determine the likely DNM size. if False, we'll
            try to leverage both parent-of-origin and haplotype-of-origin information
            to determine size. Defaults to False.

    Returns:
        int: the likely de novo expansion or contraction size, with respect to the 
            parental allele that likely mutated.
    """

    # get length of child's de novo AL
    child_als = list(map(int, row["child_AL"].split(",")))
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]

    # figure out the inferred phase of the site, if not unknown
    phase = row["phase_summary"].split(":")[0]

    # if phase is unknown, we *could* just use the minimum difference between the
    # de novo AL and any of the parental ALs, assuming that the smallest
    # difference is the most parsimonious expansion/contraction. 
    parent_als = []
    for parent in ("father", "mother"):
        parent_als.extend(list(map(int, row[f"{parent}_AL"].split(","))))
    # calculate the absolute difference between the child's AL and each of
    # the parents' four ALs.
    abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
    denovo_diffs = np.array([denovo_al - al for al in parent_als])
    min_diff_idx = np.argmin(abs_denovo_diffs)

    if use_parsimony:
        return denovo_diffs[min_diff_idx]
    
    # however, the parsimony approach is a little naive. if we have information 
    # about the parent-of-origin (and the haplotype-of-origin), we can do much better.
    else:
        # get de novo size using precise haplotype information
        parent_al = row["allele_length_summary"]
        if parent_al != -1: # condition when we don't know it
            return denovo_al - parent_al
        else:
            if phase == "unknown":
                return denovo_diffs[min_diff_idx]
            # if we know the phase of the de novo mutation, we can
            # do a parsimony-based check, but limit ourselves to the
            # allele lengths *just* in that parent-of-origin
            else:
                if phase == "dad":
                    parent_als = list(map(int, row["father_AL"].split(",")))
                elif phase == "mom":
                    parent_als = list(map(int, row["mother_AL"].split(",")))

                abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
                denovo_diffs = np.array([denovo_al - al for al in parent_als])
                min_diff_idx = np.argmin(abs_denovo_diffs)
                return denovo_diffs[min_diff_idx]


def determine_simplified_motif_size(row: pd.Series) -> str:
    """determine the "simplified" motif size at each locus. here,
    we distinguish between loci that comprise two different motif "classes" --
    e.g., both STRs and VNTRs, vs pure STRs/VNTRs

    Args:
        row (pd.Series): pandas Series

    Returns:
        str: simplified motif label
    """

    is_compound = row["n_motifs"] > 1
    # if the locus isn't compound (e.g., if it doesn't comprise
    # more than one motif), just report whether the locus contains
    # a single STR or VNTR. 
    if not is_compound:
        assert row["min_motiflen"] == row["max_motiflen"]
        return "STR" if row["min_motiflen"] <= 6 else "VNTR"
    # if the locus *is* compound, check to see if the smallest
    # and larger motifs are all <= 6 bp (i.e., STRs).
    else:
        min_motif_len = row["min_motiflen"]
        max_motif_len = row["max_motiflen"]
        # all STR
        if max_motif_len <= 6: return "STR"
        # mix of STR and VNTR
        elif min_motif_len <= 6 and max_motif_len > 6: return "complex"
        # all VNTR
        elif min_motif_len > 6: return "VNTR"
        else: return "?"


def filter_mutation_dataframe(
    mutations: pd.DataFrame,
    remove_inherited: bool = False,
    parental_overlap_frac_max: float = 1,
    denovo_coverage_min: int = 1,
    child_ratio_min: float = 0.,
    depth_min: int = 10,
) -> pd.DataFrame:
    """filter a "raw" TRGT-denovo output dataframe.

    Args:
        mutations (pd.DataFrame): pd.DataFrame with TRGT-denovo output for all
            assayed TR loci.
        remove_inherited (bool, optional): whether to remove TR loci at which the 
            de novo TR AL is present in either parents' ALs. Defaults to False.
        parental_overlap_frac_max (float, optional): a threshold for the fraction
            of reads in either parent that support the de novo AL in the child.. Defaults to 1.
        denovo_coverage_min (int, optional): number of reads in the child that must
            support the de novo AL. Defaults to 1.
        child_ratio_min (float, optional): fraction of reads in the child that must
            support the de novo AL. Defaults to 0..
        depth_min (int, optional): number of reads that must support the TRGT genotype
            in each member of the trio. Defaults to 10.

    Returns:
        pd.DataFrame: filtered pd.DataFrame object.
    """

    # validate schema of input dataframe
    TRGTDeNovoSchema.validate(mutations)

    # if we're looking for de novos, do a quick initial filter
    # of the dataframe to reduce the subsequent search space
    if denovo_coverage_min > 0:
        mutations = mutations[
            (mutations["denovo_coverage"] >= denovo_coverage_min)
            & (mutations["child_ratio"] > child_ratio_min)
        ]
        if remove_inherited:
            # figure out if the de novo AL is observed in either parent
            mutations["denovo_al_in_parents"] = mutations.apply(lambda row: al_in_parents(row), axis=1)
            mutations = mutations[~mutations["denovo_al_in_parents"]]

    # only look at autosomes + X + Y
    mutations = mutations[~mutations["trid"].str.contains("Un|random", regex=True)]

    # remove sites where there's evidence of parental dropout
    mutations = mutations[
        (mutations["father_dropout"] == "N") & (mutations["mother_dropout"] == "N")
    ]

    # annotate with depth and filter on depth
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(
        lambda a: sum(list(map(int, a.split(","))))
    )
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(
        lambda a: sum(list(map(int, a.split(","))))
    )
    mutations = mutations.query(
        f"child_coverage >= {depth_min} and mom_dp >= {depth_min} and dad_dp >= {depth_min}"
    )

    # filter on the fraction of parental reads that support the de novo allele.
    # calculate the total number of reads on either parental haplotype that support the
    # de novo AL in the child
    mutations["dad_ev"] = mutations["father_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))
    mutations["mom_ev"] = mutations["mother_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))

    # calculate the total amount of parental evidence for the de novo AL, aggregated
    # across both parents. calculate the total amount of parental depth at the site.
    mutations["parental_ev"] = mutations["dad_ev"] + mutations["mom_ev"]
    mutations["parental_dp"] = mutations["dad_dp"] + mutations["mom_dp"]

    mutations["parental_ev_frac"] = mutations["parental_ev"].values / mutations["parental_dp"].values
    mutations = mutations[mutations["parental_ev_frac"] <= parental_overlap_frac_max]

    return mutations

def infer_combined_phase(row: pd.Series):
    phase_3gen = row["phase"]
    phase_2gen = row["phase_summary"]
    n_upstream, n_downstream = row["n_upstream"], row["n_downstream"]
    consistency = row["most_common_freq"]

    combined_phase = "unknown"
    # if we have phases for both...
    if phase_2gen != "unknown" and phase_3gen != "unknown":
        if phase_2gen.split(":")[0] == phase_3gen: combined_phase = phase_2gen.split(":")[0]
        else: combined_phase = "conflicting"
    elif phase_2gen == "unknown" and phase_3gen != "unknown":
        if consistency >= 0.75: combined_phase = phase_3gen
        else: combined_phase = "unknown"
    elif phase_2gen != "unknown" and phase_3gen == "unknown":
        if n_upstream + n_downstream >= 1: return phase_2gen.split(":")[0]
        else: combined_phase = "unknown"
    
    return combined_phase
