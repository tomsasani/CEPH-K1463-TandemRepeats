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
    bootstrapped = np.zeros((n, X.shape[0], X.shape[1]))

    i = X.shape[0]
    for boot_i in numba.prange(n):
        idxs = np.random.randint(0, i, size=i)
        resampled = X[idxs, :]
        bootstrapped[boot_i, :, :] = resampled
    return bootstrapped


def classify_motif(motif: str):
    motifs = motif.split(",")
    if len(motifs) == 1:
        if len(motifs[0]) == 1:
            return "homopolymer"
        else:
            return "non-homopolymer"
    else:
        # if any([len(m) for m in motifs]) == 1:
        #     return "contains homopolymer"
        # else:
        return "non-homopolymer"


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

    # add a transmission column if it's not in there
    if "n_with_denovo_allele_strict" in res_df.columns:
        res_df["Is transmitted?"] = res_df["n_with_denovo_allele_strict"].apply(
            lambda n: "transmitted" if n > 0 else "untransmitted"
        )
    else:
        res_df["Is transmitted?"] = "unknown"

    res_df["exp_diff"] = res_df["exp_allele_diff_ref"] - res_df["exp_allele_diff_alt"]
    res_df = res_df.query("exp_diff != 0")

    res_df["Motif type"] = res_df["motif"].apply(lambda m: classify_motif(m))

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

    val_col = "measured_allele_length_error"

    if args.tech == "ont":
        kid_df = kid_df[(kid_df[val_col] >= -100) & (kid_df[val_col] <= 100)]

    row_col, col_col, group_col = "Is transmitted?", "is_ref", "Motif type"

    row_vals, col_vals = kid_df[row_col].unique(), kid_df[col_col].unique()
    row2idx, col2idx = dict(zip(row_vals, range(len(row_vals)))), dict(
        zip(col_vals, range(len(col_vals))),
    )

    f1, axarr1 = plt.subplots(
        len(row_vals),
        len(col_vals),
        sharex=True,
        sharey=True,
        figsize=(12, 7),
    )

    f2, axarr2 = plt.subplots(
        len(row_vals),
        len(col_vals),
        sharex=True,
        sharey=True,
        figsize=(12, 7),
    )

    for sub, kid_df_sub in kid_df.groupby([row_col, col_col, group_col]):
        row_val, col_val, group_val = sub
        ax_i, ax_j = row2idx[row_val], col2idx[col_val]
        if kid_df_sub.shape[0] <= 1: continue

        
        # figure out min and max of error
        emin, emax = kid_df_sub[val_col].min(), kid_df_sub[val_col].max()

        # and number of loci
        n_loci = kid_df_sub["region"].nunique()
        loc2idx = dict(zip(kid_df_sub["region"].unique(), range(n_loci)))
        kid_df_sub["region_idx"] = kid_df_sub["region"].apply(lambda r: loc2idx[r])

        # if col_val == "ALT" and group_val == "homopolymer":
        #     print (kid_df_sub[(kid_df_sub[val_col] == 0) & (kid_df_sub["Read support"] > 5)])

        arr = np.zeros((n_loci, emax - emin), dtype=np.int16)

        # increment the array
        region_idxs = kid_df_sub["region_idx"].values
        measured_al_idxs = kid_df_sub[val_col].values - emin - 1
        vals = kid_df_sub["Read support"].values

        arr[region_idxs, measured_al_idxs] = vals

        arr_sums = np.nansum(arr, axis=1)
        arr_fracs = arr / arr_sums.reshape(-1, 1)
        bootstrap_fracs = bootstrap(arr_fracs, n=args.n_bootstraps)

        # measure fraction of loci at which X% of reads perfectly support the allele
        n_pass_loci = np.sum(arr_fracs[:, 0 - emin - 1] >= 0.5)
        frac_pass_loci = round(100 * n_pass_loci / arr_fracs.shape[0], 1)

        xvals = np.arange(0, 1.01, 0.01)
        yvals = np.zeros(xvals.shape[0])
        yvals_l, yvals_u = np.zeros(xvals.shape[0]), np.zeros(xvals.shape[0])
        for i, pct in enumerate(xvals):
            n_pass_loci = np.sum(arr_fracs[:, 0 - emin - 1] >= pct)
            yval = n_pass_loci / arr_fracs.shape[0]
            n_pass_loci_bs = np.sum(bootstrap_fracs[:, :, 0 - emin - 1] >= pct, axis=1)
            frac_pass_loci_bs = n_pass_loci_bs / arr_fracs.shape[0]
            yval_l = np.percentile(frac_pass_loci_bs, q=2.5)
            yval_u = np.percentile(frac_pass_loci_bs, q=97.5)
            yvals[i] = yval
            yvals_l[i] = yval_l
            yvals_u[i] = yval_u

        axarr2[ax_i, ax_j].plot(xvals, yvals, label=f"{group_val} (n = {kid_df_sub.shape[0]})")
        axarr2[ax_i, ax_j].fill_between(xvals, yvals_l, yvals_u, alpha=0.25)

        axarr2[ax_i, ax_j].set_xlabel(
            "Fraction of reads that perfectly support the exp. AL"
        )
        axarr2[ax_i, ax_j].set_ylabel("Fraction of loci")
        axarr2[ax_i, ax_j].set_title(
            f"Read support for the {col_val} allele at {row_val} loci"
        )
        axarr2[ax_i, ax_j].legend(title="Motif type", frameon=False)
        sns.despine(ax=axarr2[ax_i, ax_j], top=True, right=True)

        # bootstrap_means = np.nanmean(bootstrap_fracs, axis=(0, 1))
        bs_fracs_l = np.nanpercentile(
            np.nanmean(bootstrap_fracs, axis=1), q=2.5, axis=0
        )
        bs_fracs_u = np.nanpercentile(
            np.nanmean(bootstrap_fracs, axis=1), q=97.5, axis=0
        )

        x = np.arange(emin, emax) + 1

        axarr1[ax_i, ax_j].plot(x, np.nanmean(arr_fracs, axis=0), label=f"{group_val} (n = {kid_df_sub.shape[0]})")
        axarr1[ax_i, ax_j].fill_between(
            x,
            bs_fracs_l,
            bs_fracs_u,
            alpha=0.25,  # color="gainsboro"
        )

        if ax_i == len(row_vals) - 1:
            axarr1[ax_i, ax_j].set_xlabel(
                f"Difference between TRGT AL and AL estimated using {args.tech.upper()} reads (bp)"
            )
        if ax_j == 0:
            axarr1[ax_i, ax_j].set_ylabel("Fraction of reads at a locus\n(+/- 95% bootstrap CI)")
        axarr1[ax_i, ax_j].set_title(
            f"Reads supporting the {col_val} allele at {row_val} loci\n(at least 50% of reads match exp. AL at {frac_pass_loci}% of loci)"
        )
        axarr1[ax_i, ax_j].set_title(
            f"Reads supporting the {col_val} allele at\n{row_val} candidate de novo STRs"
        )
        axarr1[ax_i, ax_j].legend(title="Motif type", frameon=False)
        sns.despine(ax=axarr1[ax_i, ax_j], top=True, right=True)
    f1.tight_layout()
    f1.savefig(args.out, dpi=300)
    new_fh = ".".join(args.out.split(".")[:-1])
    f2.tight_layout()
    f2.savefig(new_fh + ".alt.png", dpi=300)


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
    p.add_argument(
        "-n_bootstraps",
        type=int,
        default=1_000,
        help="""Number of bootstraps to use for bootstrap estimation.""",
    )
    args = p.parse_args()
    main(args)
