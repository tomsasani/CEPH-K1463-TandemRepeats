import pandas as pd
import argparse
from typing import List
from sklearn.mixture import GaussianMixture
import numpy as np
import math

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

def check_for_parental_evidence(row: pd.Series):

    mom_diffs, dad_diffs, kid_diffs = [], [], []
    for col in ("mom_evidence", "dad_evidence", "kid_evidence"):
        for diff_count in row[col].split("|"):
            diff, count = list(map(int, diff_count.split(":")))
            if col == "mom_evidence":
                mom_diffs.extend([diff] * count)
            elif col == "dad_evidence":
                dad_diffs.extend([diff] * count)
            else:
                kid_diffs.extend([diff] * count)

    n_alleles = len(set(kid_diffs))
    if n_alleles > 2: n_alleles = 2
    if len(set(kid_diffs)) == 1:
        return "child_not_heterozygous"
    allele_overlap = np.ones((n_alleles, 2))

    X = np.array(kid_diffs)

    gm = GaussianMixture(n_components=n_alleles).fit(X.reshape(-1, 1))
    mixture_means = gm.means_ 
    mixture_cov = gm.covariances_
    # calculate standard deviation from covariance matrix
    # standard deviation (sigma) is sqrt of the product of the residuals with themselves 
    # (i.e., sqrt of the mean of the trace of the covariance matrix)
    mixture_std = np.array([np.sqrt(np.trace(mixture_cov[i]) / n_alleles) for i in range(n_alleles)])
    
    denovo_al = int(row["exp_allele_diff_denovo"])
    non_denovo_al = int(row["exp_allele_diff_non_denovo"])

    denovo_overlap, non_denovo_overlap = False, False
    for allele_i in range(n_alleles):
        mean, std = mixture_means[allele_i][0], 2 * mixture_std[allele_i]
        lo, hi = mean - std, mean + std
        # check that at least one of the allele observations in the kid
        # overlaps the expected de novo size
        if lo <= denovo_al <= hi:
            denovo_overlap = True
        if lo <= non_denovo_al <= hi:
            non_denovo_overlap = True
        for parent_i, parent_diffs in enumerate((mom_diffs, dad_diffs)):
            X_p = np.array(parent_diffs)
            in_parent = np.sum((X_p >= lo) & (X_p <= hi))
            allele_overlap[allele_i, parent_i] = in_parent / X_p.shape[0]

    #if denovo_overlap: return "pass"
    #else: return "no_denovo_overlap"
    # for each of the two alleles in the kid, return the maximum overlap
    # with a parent's reads. the minimum value of this array will indicate
    # the minimum amount of read overlap with a parent for one of the two
    # alleles in the kid.
    return f"pass|{str(round(np.min(np.max(allele_overlap, axis=1)), 2))}"


def main(args):

    mutations = pd.read_csv(
        args.mutations,
        sep="\t"
    )
    ortho_evidence = pd.read_csv(args.orthogonal_evidence)

    res: List[pd.DataFrame] = []

    for i, row in ortho_evidence.iterrows():
        row_dict = row.to_dict()
        parental_overlap = check_for_parental_evidence(row)
        row_dict.update({"validation_status": parental_overlap})
        res.append(row_dict)

    ortho_validation = pd.DataFrame(res)

    # merge mutations with orthogonal validation
    mutations = mutations.merge(ortho_validation, how="left").fillna({"validation_status": "no_element_data"})
    
    mutations.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--orthogonal_evidence")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
