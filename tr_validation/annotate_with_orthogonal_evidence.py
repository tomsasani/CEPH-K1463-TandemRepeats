import pandas as pd
import argparse
from typing import List
from sklearn.mixture import GaussianMixture
import numpy as np


def assign_allele(row: pd.Series):

    # figure out whether the observed diff is closer
    # to the expected ALT or REF allele size
    denovo_distance = abs(row["diff"] - row["exp_allele_diff_denovo"])
    non_denovo_distance = abs(row["diff"] - row["exp_allele_diff_non_denovo"])

    # if the diff is closer to the REF allele length, we'll assume
    # that it is a REF observation. and vice versa
    #if row["sample"] in ("mom", "dad"):
    if denovo_distance < non_denovo_distance:
        is_denovo = True
    elif non_denovo_distance < denovo_distance:
        is_denovo = False
    else:
        return "unk"

    return "denovo" if is_denovo else "non_denovo"

def check_for_parental_evidence(df_sub: pd.DataFrame):

    # check that every member of trio has sufficient read
    # support
    if not all([p in df_sub["sample"].unique() for p in ("mom", "dad", "kid")]):
        return "incomplete_trio"

    # keep track of what fraction of parental reads overlap the denovo
    # range in the child's reads. shape is (n_alleles, n_parents)
    allele_overlap = np.ones((2, 2))

    # fit gaussian mixture to child's reads
    diffs = df_sub[df_sub["sample"] == "kid"]["diff"].values
    reps = df_sub[df_sub["sample"] == "kid"]["Read support"].values
    X = np.repeat(diffs, reps)
    if X.shape[0] <= 1:
        return np.max(allele_overlap, axis=1)#"insufficient_child_depth"
    if np.unique(X).shape[0] == 1:
        return np.max(allele_overlap, axis=1)#"child_not_heterozygous"

    gm = GaussianMixture(n_components=2).fit(X.reshape(-1, 1))
    mixture_means = gm.means_ 
    mixture_cov = gm.covariances_
    mixture_std = np.array([np.sqrt(np.trace(mixture_cov[i]) / 2) for i in range(2)])
    # figure out if either parent has a read that overlaps a normal distribution
    # centered around the child's means + stdevs
    lo, hi = mixture_means - (2 * mixture_std), mixture_means + (2 * mixture_std)

    for allele_i in range(2):
        mean, std = mixture_means[allele_i], mixture_std[allele_i]
        lo, hi = mean - std, mean + std
        for parent_i, parent in enumerate(("mom", "dad")):
            p_diffs = df_sub[df_sub["sample"] == parent]["diff"].values
            p_reps = df_sub[df_sub["sample"] == parent]["Read support"].values
            X_p = np.repeat(p_diffs, p_reps)
            in_parent = np.sum((X_p >= lo) & (X_p <= hi))
            allele_overlap[allele_i, parent_i] += in_parent / np.sum(p_reps)
    total_allele_overlap = np.sum(np.sum(allele_overlap, axis=1) > 0)
    return np.max(allele_overlap, axis=1)#"pass" if total_allele_overlap <= 1 else "both_child_alleles_in_parents"

def is_validated(df: pd.DataFrame):


    status = []
    # make sure we have Element read evidence 
    # in every member of the trio
    if not all([m in df["sample"].unique() for m in ("mom", "dad", "kid")]): 
        status.append("trio_coverage")
    
    for sample, sample_df in df.groupby("sample"):
        # check if there's read support for the denovo in mom or dad
        if sample in ("mom", "dad"):
            if "denovo" in sample_df["allele_assignment"].to_list():
                status.append(f"dnm_in_{sample}")
        # make sure there's read support for the denovo in the kid
        if sample == "kid":
            if "denovo" not in sample_df["allele_assignment"].to_list():
                status.append("no_dnm_in_kid")

    if len(status) == 0: is_validated = "Y"
    else: is_validated = "|".join(list(set(status)))
    return is_validated


def main(args):

    mutations = pd.read_csv(
        args.mutations,
        sep="\t",
        dtype={"sample_id": str},
    )

    ortho_evidence = pd.read_csv(args.orthogonal_evidence)

    # subset element read depth information to the child, and calculate
    # total read support
    GROUP_COLS = [
        "trid",
        "sample", # whether it's mom, dad or kid
        "genotype",
    ]

    # total number of reads in mom, dad and the child
    ortho_totals = (
        ortho_evidence.groupby(GROUP_COLS)
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    ortho_evidence = ortho_evidence.merge(ortho_totals, how="left")

    ortho_evidence["Read support frac"] = ortho_evidence["Read support"] / ortho_evidence["Total read support"]

    # for every observed number of net CIGAR operations in a read,
    # figure out if that net CIGAR is closer to the expected net CIGAR
    # operation for a denovo or non-denovo allele.
    ortho_evidence["allele_assignment"] = ortho_evidence.apply(
        lambda row: assign_allele(row),
        axis=1,
    )

    # for every TRID, ask if the site is "validated." this requires at least
    # one read supporting the denovo allele in the kid and zero reads supporting the
    # denovo allele in the parents
    res: List[pd.DataFrame] = []

    for (trid, genotype), sub_df in ortho_evidence.groupby(["trid", "genotype"]):
        #is_valid = is_validated(sub_df)
        parental_overlap = check_for_parental_evidence(sub_df)
        sub_df["parental_overlap"] = parental_overlap
        sub_df = sub_df[sub_df["sample"] == "kid"].drop_duplicates(["trid", "genotype"])
        res.append(sub_df)

    ortho_validation = pd.concat(res)
    
    # merge mutations with orthogonal validation
    mutations = mutations.merge(ortho_validation, how="left")
    mutations.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--orthogonal_evidence")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
