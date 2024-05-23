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


def add_likely_de_novo_size(row: pd.Series):

    child_als = list(map(int, row["child_AL"].split(",")))
    # denovo idx
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]

    if not row["trid"].startswith("chrX"):
        non_denovo_idx = 1 - denovo_idx
        non_denovo_al = child_als[non_denovo_idx]

    # figure out the inferred phase of the site, if not unknown
    phase = row["phase_summary"].split(":")[0]

    # if phase is unknown, we just use the minimum difference between the
    # de novo AL and any of the four parental ALs, assuming that the smallest
    # difference is the most parsimonious expansion/contraction
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
    # only look at autosomes + X
    mutations = mutations[~mutations["trid"].str.contains("Un|random", regex=True)]
    # remove pathogenics
    mutations = mutations[~mutations["trid"].str.contains("pathogenic")]
    # remove non-TRF/TRSOLVE
    mutations = mutations[mutations["trid"].str.contains("trsolve|TRF", regex=True)]


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
        mutations = mutations.drop_duplicates(dup_cols)

    mutations["dad_ev"] = mutations["father_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))
    mutations["mom_ev"] = mutations["mother_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))

    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(lambda o: sum(list(map(int, o.split(",")))))
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(lambda o: sum(list(map(int, o.split(",")))))

    # calculate fraction of reads supporting the denovo in the parents
    mutations["parental_overlap_coverage_total"] = mutations["dad_ev"] + mutations["mom_ev"]
    mutations["parental_coverage_total"] = mutations["mom_dp"] + mutations["dad_dp"]
    mutations["parental_overlap_coverage_frac"] = mutations["parental_overlap_coverage_total"] / mutations["parental_coverage_total"]
    mutations = mutations.query(f"parental_overlap_coverage_frac <= {parental_overlap_frac_max}")

    # if we want to remove grandparental evidence
    if remove_gp_ev and "gp_ev" in mutations.columns:
        mutations = mutations[~mutations["gp_ev"].str.contains("Y")]

    # mutations["denovo_size"] = mutations.apply(lambda row: add_likely_de_novo_size(row), axis=1)

    return mutations
