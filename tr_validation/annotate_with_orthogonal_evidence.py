import pandas as pd
import argparse
from typing import List
from sklearn.mixture import GaussianMixture
import numpy as np
import math
import utils
from collections import defaultdict


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

def annotate_with_concordance(row: pd.Series) -> str:
    """given orthogonal evidence at a locus, figure out whether
    the orthogonal evidence is concordant with the expected genotype
    in the child.

    Args:
        row (pd.Series): pd.Series object

    Returns:
        str: annotation label
    """

    # gather diffs in all members of the trio
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

    # figure out the difference between the expected
    # sizes of the de novo and non-de novo alleles
    denovo_al = int(row["exp_allele_diff_denovo"])
    non_denovo_al = int(row["exp_allele_diff_non_denovo"])
    exp_diff = denovo_al - non_denovo_al
    # how many alleles do we expect to see in the kid? i.e., is the kid HOM or HET?
    n_alleles_exp = 1 if exp_diff == 0 else 2

    # figure out the observed number of distinct alleles in the kid
    n_alleles = len(set(kid_diffs))
    if n_alleles_exp > 1 and n_alleles == 1:
        return "not_ibs"
    
    # validate by asking if the de novo allele is observed in the
    # child and absent from both parents

    out_dict = defaultdict(int)

    for diffs, label in zip((kid_diffs, mom_diffs, dad_diffs), ["kid", "mom", "dad"]):

        dn_overlap = 0
        # allow for off-by-one errors when asking if kid has de novo.
        # don't allow for off-by-one errors when asking if parents have de novo.
        for d in diffs:
            # if label == "kid":
            if denovo_al == d: dn_overlap += 1
            # else:
            #     if d == denovo_al: dn_overlap += 1
        out_dict[label] = dn_overlap

        # fit a Guassian mixture model to the child's allele length
        # distribution with the expected number of components equal to
        # the expected number of alleles
        # X = np.array(diffs)
        # gm = GaussianMixture(n_components=n_alleles_exp).fit(X.reshape(-1, 1))
        # get the means and covariance sof the two components
        # mixture_means = gm.means_ 
        # mixture_cov = gm.covariances_
        # calculate standard deviation from covariance matrix
        # standard deviation (sigma) is sqrt of the product of the residuals with themselves
        # (i.e., sqrt of the mean of the trace of the covariance matrix)
        # mixture_std = np.array(
        #     [np.sqrt(np.trace(mixture_cov[i]) / n_alleles_exp) for i in range(n_alleles_exp)]
        # )

        # check to see if the kid's mixture components +/- STDEV overlap the expected
        # allele sizes. assume this locus passes by default.
        # overlap = 0
        # for allele_i in range(n_alleles_exp):
        #     # allow for 1bp off
        #     mean, std = round(mixture_means[allele_i][0]), round(max([2 * mixture_std[allele_i], 1]))
        #     lo, hi = mean - std, mean + std
        #     if row["region"] == "chr17:15781239-15781258":
        #         print (lo, hi, denovo_al, non_denovo_al)
        #         print (lo <= denovo_al <= hi, lo <= non_denovo_al <= hi)
        #     # if this component doesn't overlap *either* the de novo or non denovo allele,
        #     # set overlap to False
        #     if (lo <= denovo_al <= hi) or (lo <= non_denovo_al <= hi): overlap += 1
        #     # if not (lo <= denovo_al <= hi) and not (lo <= non_denovo_al <= hi):
        #     #     overlap = False

    if "chr2_1578" in row["trid"] and row["sample_id"] == "2211": print (out_dict)

    return "pass" if (out_dict["kid"] >= 1 and out_dict["mom"] == 0 and out_dict["dad"] == 0) else "fail"


def main(args):

    mutations = pd.read_csv(
        args.mutations,
        sep="\t"
    )
    ortho_evidence = pd.read_csv(args.orthogonal_evidence)

    res: List[pd.DataFrame] = []

    for i, row in ortho_evidence.iterrows():
        row_dict = row.to_dict()
        parental_overlap = annotate_with_concordance(row)
        if "chr2_" in row_dict["trid"] and row_dict["sample_id"] == "2211": print (parental_overlap)
        row_dict.update({"validation_status": parental_overlap})
        res.append(row_dict)

    ortho_validation = pd.DataFrame(res)

    # merge mutations with orthogonal validation
    mutations = mutations.merge(ortho_validation, how="left").fillna({"validation_status": "no_data", "phase_summary": "unknown"})
    mutations["likely_denovo_size"] = mutations.apply(lambda row: utils.add_likely_de_novo_size(row), axis=1)
    mutations.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--orthogonal_evidence")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
