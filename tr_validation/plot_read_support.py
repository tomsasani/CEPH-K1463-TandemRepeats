import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable
import numba

plt.rc("font", size=10)


@numba.njit
def compute_nanmean(a: np.ndarray, row: bool = True) -> np.ndarray:
    """Compute the sum of a 2D numpy array
    on a per-row basis, ignoring nans. Since `numba` does not
    support kwargs in the `np.nansum` function, it's
    necessary to piece this out into its own function
    so that it can be decorated with `numba.njit`.

    Args:
        a (np.ndarray): A 2D numpy array of size (N, M).

        row (bool, optional): Whether to calculate means by row or column. Defaults to True.


    Returns:
        rowsums (np.ndarray): A 1D numpy array of size (N, ) containing \
            sums of every row in the input.
    """
    idx = 0 if row else 1
    empty_a = np.zeros(a.shape[idx])
    for i in np.arange(a.shape[idx]):
        empty_a[i] = np.nanmean(a[i]) if row else np.nanmean(a[:, i])
    return empty_a


@numba.njit(parallel=True)
def bootstrap(X: np.ndarray, n: int = 100):
    means = np.zeros((n, X.shape[1]))

    i = X.shape[0]
    for boot_i in numba.prange(n):
        idxs = np.random.randint(0, i, size=i)
        resampled = X[idxs, :]
        mean = compute_nanmean(resampled, row=False)
        means[boot_i, :] = mean
    return means


def classify_motif(motif: str):
    motif_len = len(motif)
    if motif_len < 5:
        return str(motif_len)
    else:
        return "5+"


def calculate_min_dist(row: pd.Series):
    # figure out whether the observed diff is closer
    # to the expected ALT or REF allele size
    alt_distance = row["diff"] - row["exp_allele_diff_alt"]
    ref_distance = row["diff"] - row["exp_allele_diff_ref"]

    # if the diff is closer to the REF allele length, we'll assume
    # that it is a REF observation. and vice versa
    if abs(ref_distance) < abs(alt_distance):
        is_ref = True
    elif abs(ref_distance) > abs(alt_distance):
        is_ref = False
    # if they're equal, we'll randomly assign the diff to be a REF/ALT obs
    else:
        is_ref = np.random.uniform() > 0.5

    return "REF" if is_ref else "ALT"


def main(args):
    numba.set_num_threads(args.threads)

    res_df = pd.read_csv(args.csv)

    # add a homopolymer column if it's not in there
    if "Is homopolymer?" not in res_df:
        res_df["Is homopolymer?"] = res_df["motif"].apply(lambda m: len(m) == 1)

    # for each read, calculate the minimum distance to either the REF or ALT allele
    # (i.e., the "closest" allele)
    res_df["is_ref"] = res_df.apply(lambda row: calculate_min_dist(row), axis=1)
    res_df["Expected allele length"] = res_df.apply(
        lambda row: row["exp_allele_diff_ref"]
        if row["is_ref"] == "REF"
        else row["exp_allele_diff_alt"],
        axis=1,
    )
    res_df["measured_allele_length_error"] = (
        res_df["diff"] - res_df["Expected allele length"]
    )

    # subset to kid
    kid_df = res_df[res_df["sample"] == "kid"]
    kid_df_totals = (
        kid_df.groupby("region")
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    kid_df = kid_df.merge(kid_df_totals, on="region", how="left")
    # require at least 10 total spanning reads
    kid_df = kid_df[kid_df["Total read support"] >= 10]

    col = "measured_allele_length_error"

    f, axarr = plt.subplots(2, 2, sharex=True, sharey=False, figsize=(12, 6))

    group2ax = {
        ("REF", True): (0, 0),
        ("REF", False): (0, 1),
        ("ALT", True): (1, 0),
        ("ALT", False): (1, 1),
    }

    for sub, kid_df_sub in kid_df.groupby(["is_ref", "Is homopolymer?"]):
        is_ref, is_homo = sub
        ax_i, ax_j = group2ax[(is_ref, is_homo)]

        # figure out min and max of error
        emin, emax = kid_df_sub[col].min(), kid_df_sub[col].max()
        # and number of loci
        n_loci = kid_df_sub["region"].nunique()
        loc2idx = dict(zip(kid_df_sub["region"].unique(), range(n_loci)))
        kid_df_sub["region_idx"] = kid_df_sub["region"].apply(lambda r: loc2idx[r])

        arr = np.zeros((n_loci, emax - emin), dtype=np.int16)

        # increment the array
        region_idxs = kid_df_sub["region_idx"].values
        measured_al_idxs = kid_df_sub[col].values - emin - 1
        vals = kid_df_sub["Read support"].values

        arr[region_idxs, measured_al_idxs] = vals

        arr_sums = np.nansum(arr, axis=1)
        arr_fracs = arr / arr_sums.reshape(-1, 1)

        bs_fracs = bootstrap(arr_fracs)
        bs_fracs_l = np.nanpercentile(bs_fracs, q=2.5, axis=0)
        bs_fracs_u = np.nanpercentile(bs_fracs, q=97.5, axis=0)

        x = np.arange(emin, emax) + 1

        axarr[ax_i, ax_j].plot(x, np.nanmean(arr_fracs, axis=0), c="grey")
        axarr[ax_i, ax_j].fill_between(
            x, bs_fracs_l, bs_fracs_u, alpha=0.75, color="gainsboro"
        )

        axarr[ax_i, ax_j].set_xlabel("Difference between Element and TRGT AL (bp)")
        axarr[ax_i, ax_j].set_ylabel("Fraction of reads at a locus")
        homopol_status = "homopolymers" if is_homo else "non-homopolymers"
        axarr[ax_i, ax_j].set_title(
            f"Reads supporting the\n{is_ref} allele at {homopol_status}"
        )
        sns.despine(ax=axarr[ax_i, ax_j], top=True, right=True)
    f.tight_layout()
    f.savefig(args.out, dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--csv",
        type=str,
        help="""Path to CSV file containing estimated allele lengths.""",
    )
    p.add_argument(
        "--out",
        type=str,
        help="""Path to output image.""",
    )
    p.add_argument(
        "-tech",
        type=str,
        default="illumina",
        help="""Technology used to generate allele length estimates.""",
    )
    p.add_argument(
        "-pctile",
        type=int,
        default=95,
        help="""Percentile to plot surrounding mean estimates using a bootstrap.""",
    )
    p.add_argument(
        "-threads",
        type=int,
        default=1,
        help="""Number of threads to use for bootstrap estimation.""",
    )
    args = p.parse_args()
    main(args)
