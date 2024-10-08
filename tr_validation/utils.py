from collections import Counter
from typing import Tuple, Union, List
import pandas as pd

import pysam
import numpy as np


MATCH, INS, DEL = range(3)
OP2DIFF = {MATCH: 0, INS: 1, DEL: -1}
# from SAM spec (https://samtools.github.io/hts-specs/SAMv1.pdf)
# page 8 -- 1 if the operation consumes reference sequence
CONSUMES_REF = dict(zip(range(9), [1, 0, 1, 1, 0, 0, 0, 1, 1]))

def al_in_parents(row: pd.Series):
    child_als = row["child_AL"].split(",")
    # denovo idx
    denovo_idx = int(row["index"])

    mom_als, dad_als = row["mother_AL"].split(","), row["father_AL"].split(",")

    # if we're on a male X, there's only one haplotype so the "non-denovo idx"
    # is not applicable. we only have to ask if the denovo AL is among the mom's ALs
    if row["suffix"].startswith("S") and row["trid"].split("_")[0] in ("chrX", "chrY"):
        denovo_al = child_als[denovo_idx]
        if denovo_al in mom_als: return True
        else: return False

    else:
        non_denovo_idx = 1 - denovo_idx
        try:
            denovo_al, non_denovo_al = child_als[denovo_idx], child_als[non_denovo_idx]
        except IndexError: print (child_als, row["trid"], row["suffix"])
        if (denovo_al in mom_als) and (non_denovo_al in dad_als):
            return True
        elif (denovo_al in dad_als) and (non_denovo_al in mom_als):
            return True
        else: return False


def add_likely_de_novo_size(row: pd.Series, use_phase: bool = True):

    child_als = list(map(int, row["child_AL"].split(",")))
    # denovo idx
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]
    non_denovo_al = None

    if not row["trid"].startswith("chrX"):
        non_denovo_idx = 1 - denovo_idx
        non_denovo_al = child_als[non_denovo_idx]

    # figure out the inferred phase of the site, if not unknown
    phase = "unknown"
    if "phase_summary" in row:
        phase = row["phase_summary"].split(":")[0]

    # if phase is unknown, we *could* just use the minimum difference between the
    # de novo AL and any of the four parental ALs, assuming that the smallest
    # difference is the most parsimonious expansion/contraction. however, this is too
    # naive. for example, if the kid's non-denovo AL perfectly matches one of the two parents,
    # then we can assume inheritance of that allele from taht parent. we can call this "fuzzy"
    # phasing.

    # what if the kid perfectly inherited an AL form one parent?
    if use_phase:
        if phase == "unknown":
            parent_als = []
            for parent in ("father", "mother"):
                parent_als.extend(list(map(int, row[f"{parent}_AL"].split(","))))

            abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
            denovo_diffs = np.array([denovo_al - al for al in parent_als])
            min_diff_idx = np.argmin(abs_denovo_diffs)
            # likely original AL
            return denovo_diffs[min_diff_idx]
        else:
            if phase == "dad":
                parent_als = list(map(int, row["father_AL"].split(",")))
            elif phase == "mom":
                parent_als = list(map(int, row["mother_AL"].split(",")))

            abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
            denovo_diffs = np.array([denovo_al - al for al in parent_als])
            min_diff_idx = np.argmin(abs_denovo_diffs)

            # likely original AL
            return denovo_diffs[min_diff_idx]
    else:
        parent_als = []
        for parent in ("father", "mother"):
            parent_als.extend(list(map(int, row[f"{parent}_AL"].split(","))))

        abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
        denovo_diffs = np.array([denovo_al - al for al in parent_als])
        min_diff_idx = np.argmin(abs_denovo_diffs)
        # likely original AL
        return denovo_diffs[min_diff_idx]
    
def determine_motif_size(row: pd.Series):
    """returns index of the motif in a complex locus that is mutated"""
    
    min_motif_len = row["min_motiflen"]

    return min_motif_len

def determine_simplified_motif_size(row: pd.Series):
    """returns index of the motif in a complex locus that is mutated"""

    is_complex = row["n_motifs"] > 1
    if not is_complex:
        return "STR" if row["min_motiflen"] <= 6 else "VNTR"
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
    remove_complex: bool = False,
    remove_duplicates: bool = False,
    remove_gp_ev: bool = False,
    remove_inherited: bool = False,
    parental_overlap_frac_max: float = 1,
    denovo_coverage_min: int = 1,
    child_ratio_min: float = 0.,
    depth_min: int = 10,
):

    if denovo_coverage_min > 0:
        mutations = mutations[
            (mutations["denovo_coverage"] >= denovo_coverage_min)
           & (mutations["allele_ratio"] > 0)
           & (mutations["child_ratio"] > 0)
        ]
    else:
        mutations = mutations[
            (mutations["denovo_coverage"] >= denovo_coverage_min)
        ]
    

    mutations = mutations[mutations["child_ratio"] >= child_ratio_min]
    # only look at autosomes + X + Y
    mutations = mutations[~mutations["trid"].str.contains("Un|random", regex=True)]

    mutations["denovo_al"] = mutations.apply(lambda row: row["child_AL"].split(",")[row["index"]], axis=1)

    # remove sites where de novo AL is observed in a parent if we're interested
    # in filtering DNMs
    if denovo_coverage_min > 0:
        mutations["denovo_al_in_parents"] = mutations.apply(lambda row: al_in_parents(row), axis=1)
        if remove_inherited:
            mutations = mutations[~mutations["denovo_al_in_parents"]]

        # no parental dropout, only when we filter DNMs
        mutations = mutations[(mutations["father_dropout"] == "N") & (mutations["mother_dropout"] == "N")]

    # annotate with depth and filter on depth
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(lambda a: sum(list(map(int, a.split(",")))))
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(lambda a: sum(list(map(int, a.split(",")))))
    mutations = mutations.query(f"child_coverage >= {depth_min} and mom_dp >= {depth_min} and dad_dp >= {depth_min}")


    if remove_complex:
        mutations = mutations[~mutations["child_MC"].str.contains("_")]
    if remove_duplicates:
        dup_cols = ["trid"]
        if "sample_id" in mutations.columns:
            dup_cols.append("sample_id")
        mutations = mutations.drop_duplicates(dup_cols, keep="first")

    if parental_overlap_frac_max < 1:


        mutations["dad_ev"] = mutations["father_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))
        mutations["mom_ev"] = mutations["mother_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))

        mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(lambda o: sum(list(map(int, o.split(",")))))
        mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(lambda o: sum(list(map(int, o.split(",")))))

        mutations["parental_ev"] = mutations["dad_ev"] + mutations["mom_ev"]
        mutations["parental_dp"] = mutations["dad_dp"] + mutations["mom_dp"]


        mutations["parental_ev_frac"] = mutations["parental_ev"].values / mutations["parental_dp"].values

        mutations = mutations[mutations["parental_ev_frac"] <= parental_overlap_frac_max]

    # if we want to remove grandparental evidence
    if remove_gp_ev and "gp_ev" in mutations.columns:
        mutations = mutations[~mutations["gp_ev"].str.contains("Y")]


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