import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable
import numba
import tqdm
import matplotlib.patches as patches
from utils import filter_mutation_dataframe
import glob


plt.rc("font", size=10)

def trid_to_region(t):
    try:
        chrom, start, end, method = t.split("_")
        start, end = int(start) - 1, int(end)
        return f"{chrom}:{start}-{end}"
    except ValueError:
        return "UNK"

def classify_motif(motif: str):
    motifs = motif.split(",")
    if len(motifs) == 1:
        if len(motifs[0]) == 1:
            return "homopolymer"
        else:
            return "non-homopolymer"
    else:
        if any([len(m) for m in motifs]) == 1:
            return "contains homopolymer"
        else:
            return "non-homopolymer"


def assign_allele(row: pd.Series):
    # figure out whether the observed diff is closer
    # to the expected ALT or REF allele size
    denovo_distance = abs(row["diff"] - row["exp_allele_diff_denovo"])
    non_denovo_distance = abs(row["diff"] - row["exp_allele_diff_non_denovo"])

    # if the diff is closer to the REF allele length, we'll assume
    # that it is a REF observation. and vice versa
    if denovo_distance < non_denovo_distance:
        is_denovo = True
    elif denovo_distance > non_denovo_distance:
        is_denovo = False
    # if they're equal, we'll randomly assign the diff to be a REF/ALT obs
    else:
        is_denovo = np.random.uniform() > 0.5

    return "denovo" if is_denovo else "non_denovo"


def main():

    mutations = []
    for fh in glob.glob("tr_validation/csv/*.element.dnm.read_support.csv"):
        df = pd.read_csv(fh)
        smp = fh.split("/")[-1].split(".")[0]
        df["sample_id"] = smp
        mutations.append(df)
    mutations = pd.concat(mutations)

    DIFF_COLUMN = "exp_allele_diff_denovo"

    # filter mutation dataframe
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=True,
        remove_duplicates=False,
    )

    mutations = mutations[
        (mutations["motif"].str.len() >= 1)
        & (mutations[DIFF_COLUMN] <= 20)
        & (mutations[DIFF_COLUMN] >= -20)
        & (mutations["denovo_coverage"] >= 2)
        & (mutations["child_coverage"] >= 10)
    ]

    mutations["Is homopolymer?"] = mutations["motif"].apply(lambda m: len(m) == 1)

    # subset element read depth information to the child, and calculate
    # total read support
    GROUP_COLS = [
        "trid",
        "sample_id",
        "region",
        "genotype",
        #"exp_allele_diff_denovo",
        "exp_allele_diff_non_denovo",
        DIFF_COLUMN,
    ]
    mutations = mutations[mutations["sample"] == "kid"]
    mutation_totals = (
        mutations.groupby(GROUP_COLS)
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    mutations = mutations.merge(mutation_totals, on=GROUP_COLS, how="left")
    mutations["Read support frac"] = mutations["Read support"] / mutations["Total read support"]
    mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4" if s.startswith("200") else "G2/3")

    # require at least 10 total spanning reads and at least 2 reads supporting a given diff
    mutations = mutations[
        (mutations["Read support"] > 1)
        & (mutations["Total read support"] >= 10)
        & (mutations["Total read support"] <= 100)
    ]

    # for every observed number of net CIGAR operations in a read,
    # figure out if that net CIGAR is closer to the expected net CIGAR
    # operation for a REF or ALT allele.
    mutations["allele_assignment"] = mutations.apply(lambda row: assign_allele(row), axis=1)

    print (mutations[mutations["diff"] == mutations["exp_allele_diff_non_denovo"]])

    # mutations["Element vs. TRGT AL difference"] = mutations.apply(lambda row: row["diff"] - row[f"exp_allele_diff_{row['allele_assignment']}"], axis=1)

    # allele_support = allele_support[allele_support["Total support"] >= 10]
    mutations_sorted = mutations[mutations["allele_assignment"].isin(["denovo"])]
    mutations_sorted[DIFF_COLUMN] = np.abs(mutations_sorted[DIFF_COLUMN])

    # mutations_sorted_totals = mutations_sorted.groupby(DIFF_COLUMN).agg({"Read support frac": lambda r: sum(r >= 0.9) / len(r)})
    # print(
    #     mutations_sorted.query("diff_to_allele < -10")[
    #         [
    #             "region",
    #             "sample_id",
    #             "diff",
    #             "Element vs. TRGT AL difference",
    #             "Read support",
    #             "allele_assignment",
    #             "exp_allele_diff_alt",
    #             "exp_allele_diff_ref",
    #             "Read support frac"
    #         ]
    #     ]
    # )

    g = sns.FacetGrid(data=mutations_sorted, row="generation", aspect=2.5)
    g.map(sns.boxplot, DIFF_COLUMN, "Read support frac", color="white", fliersize=0)
    g.map(sns.stripplot, DIFF_COLUMN, "Read support frac", alpha=0.25)
    g.figure.suptitle("Fraction of Element reads that\nlikely support the de novo allele")
    g.tight_layout()
    # sns.despine(ax=ax, top=True, right=True)
    g.savefig("frac.png", dpi=200)

    # print (mutations[["genotype", "sample_id", "region", "diff", DIFF_COLUMN]])
    g = sns.FacetGrid(
        data=mutations[mutations["allele_assignment"] == "denovo"].rename(columns={"diff": "Net CIGAR operations in Element read"}),
        row="generation",
        aspect=1.25,
    )
    g.map(
        sns.scatterplot,
        "Net CIGAR operations in Element read",
        DIFF_COLUMN,
        size=mutations["Read support frac"],
        alpha=0.25,
    )
    g.add_legend(title="Fraction of reads\nsupporting the AL")
    g.figure.suptitle("Element and TRGT concordance at small de novo STRs")
    g.tight_layout()
    # size=mutations["Read support frac"],
    # alpha=0.5,

    # g.set_ylabel("Difference between Element AL and TRGT AL")
    # g.set_xlabel("Expansion/contraction size")
    g.savefig("diff.png", dpi=200)

    # f, ax = plt.subplots()
    # sns.histplot(
    #     data=allele_support[allele_support["allele_assignment"] == "ref"],
    #     hue="one_bp_expansion",
    #     multiple="dodge",
    #     # cumulative=True,
    #     #stat="proportion",
    #     x="Fraction",
    #     bins=np.arange(0, 1.05, 0.05),
    #     ax=ax,
    # )
    # ax.set_ylabel("Count of DNMs")
    # ax.set_xlabel("Fraction of Element reads that\nlikely support the REF allele")
    # f.tight_layout()
    # f.savefig("frac.png", dpi=200)

    # print (kid_df.groupby("is_concordant").size())

    # f, ax = plt.subplots()

    # total = np.sum(kid_df["Read support"])
    # for col, lab in zip(("exp_allele_diff_ref", "exp_allele_diff_alt"), ("REF", "ALT")):
    #     ax.axvline(
    #         kid_df[col].values[0],
    #         ls=":",
    #         c="firebrick" if lab == "REF" else "dodgerblue",
    #         label=f"Expected {lab} AL",
    #     )

    # actual_diffs = kid_df["diff"]
    # read_support = kid_df["Read support"]
    # vals = []
    # for diff, support in zip(actual_diffs.to_list(), read_support.to_list()):
    #     vals.extend([diff] * support)
    # mean, std = (np.mean(vals), np.std(vals) * 2)
    # # ax.add_patch(
    # #     patches.Rectangle(
    # #         (mean - std, 0),
    # #         (mean + std) - (mean - std),
    # #         max(read_support / total),
    # #         facecolor="gainsboro",
    # #         alpha=0.5,
    # #         #label="estimated AL +/- 2 stdev",
    # #         # lw=1,
    # #         # ec="k",
    # #     ),

    # # )
    # ax.bar(actual_diffs, read_support, color="gainsboro", label="Observed ALs", ec="k", lw=1)

    # # ax.axvline(mean, c="gainsboro", zorder=-1)
    # # ax.set_title(ref_status)
    # # ax.set_ylabel(f"Fraction of {ref_status}\nreads supporting AL")
    # ax.legend(frameon=False)
    # ax.set_xlabel("Difference between Element and TRGT AL (bp)")

    # # f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # # for ax, ref_status in zip((ax1, ax2), ("REF", "ALT")):
    # #     kid_df_sub = kid_df[kid_df["is_ref"] == ref_status]
    # #     total = np.sum(kid_df_sub["Read support"])
    # #     if kid_df_sub.shape[0] == 0: continue
    # #     ax.axvline(kid_df_sub["Expected allele length"].values[0], ls=":", c="firebrick", label="Expected AL")

    # #     actual_diffs = kid_df_sub["diff"]
    # #     read_support = kid_df_sub["Read support"]
    # #     vals = []
    # #     for diff, support in zip(actual_diffs.to_list(), read_support.to_list()):
    # #         vals.extend([diff] * support)
    # #     mean, std = (np.mean(vals), np.std(vals) * 2)
    # #     ax.add_patch(
    # #         patches.Rectangle(
    # #             (mean - std, 0),
    # #             (mean + std) - (mean - std),
    # #             max(read_support / total),
    # #             facecolor="gainsboro",
    # #             alpha=0.5,
    # #             label="estimated AL +/- 2 stdev",
    # #             # lw=1,
    # #             # ec="k",
    # #         ),

    # #     )
    # #     ax.bar(actual_diffs, read_support / total, color="dodgerblue")

    # #     ax.axvline(mean, c="gainsboro", zorder=-1)
    # #     ax.set_title(ref_status)
    # #     ax.set_ylabel(f"Fraction of {ref_status}\nreads supporting AL")
    # #     ax.legend(frameon=False)
    # #     ax.set_xlabel("Difference between Element and TRGT AL (bp)")
    # f.tight_layout()
    # sns.despine(ax=ax, top=True, right=True)
    # f.savefig("test.png", dpi=200)


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
    main()
