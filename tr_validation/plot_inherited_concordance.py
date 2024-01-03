import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import itertools
import numpy as np

plt.rc("font", size=14)


def classify_motif(motif: str):
    motif_len = len(motif)
    if motif_len < 6:
        return motif_len
    else:
        return 6


def main():
    genotype = "0_1"
    fh = f"tr_validation/csv/2216.{genotype}.inherited.read_support.csv"
    res_df = pd.read_csv(fh)
    
    region_sums = (
        res_df.groupby(["sample", "trid"])
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    res_df = res_df.merge(region_sums, on=["sample", "trid"])

    

    if "0_1" in fh:
        res_df = res_df.query("diff != 0")

    res_df = res_df[res_df["Total read support"] >= 10]
    res_df["Motif size"] = res_df["motif"].apply(lambda m: len(m))
    res_df = res_df[res_df["Motif size"] <= 10]

    res_df["Matches ALT?"] = res_df["diff"] == res_df["exp_allele_diff_alt"]
    res_df["Matches REF?"] = res_df["diff"] == res_df["exp_allele_diff_ref"]

    

    res_df["Fraction of reads"] = (
        res_df["Read support"] / res_df["Total read support"]
    )
    # res_df = res_df[res_df["Fraction of reads"] >= 0.05]


    order = sorted(res_df["Motif size"].unique())
    # order = list(map(str, range(1, 10)))
    # order.append("10+")

    # just plot the fractions of reads that exactly support the ALT allele
    g = sns.FacetGrid(
        data=res_df[res_df["Matches ALT?"] == True], col="sample", aspect=1
    )
    g.map(
        sns.boxplot,
        "Motif size",
        "Fraction of reads",
        color="white",
        order=order,
        fliersize=0,
    )

    g.map(
        sns.stripplot,
        "Motif size",
        "Fraction of reads",
        "Motif size",
        hue_order=order,
        palette="colorblind",
        order=order,
        alpha=0.25,
    )

    gt = "homozygous" if genotype in ("0_0", "1_1") else "heterozygous"
    g.fig.suptitle(
        f"Allele balance at {gt} STRs present in all members of a trio (n = 1,000 randomly sampled loci)",
        fontsize=10,
    )
    g.tight_layout()
    g.savefig("test.png", dpi=200)

    res_df = res_df[res_df["sample"] == "kid"]


    print (res_df[res_df["exp_allele_diff_alt"] > 40][["trid", "diff", "exp_allele_diff_ref", "exp_allele_diff_alt"]])

    res_df["Is homopolymer?"] = res_df["motif"].apply(lambda m: len(m) == 1)

    x, y = "diff", "exp_allele_diff_alt"

    g = sns.relplot(
        data=res_df,
        x=x,
        y=y,
        col="Motif size",
        hue="Fraction of reads",
        alpha=1,
        s=100,
        col_wrap=5,
    )
    g.map(plt.axline, xy1=(0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)

    titles = [
        r"$\bf{a)}$" + " candidate STR DNMs (not in homopolymers)",
        r"$\bf{b)}$" + " candidate STR DNMs in homopolymers",
    ]
    for i, ax in enumerate(g.axes.flat):
        # ax.set_title(titles[i], fontsize=10)
        ax.set_xlabel(
            "Measured allele length difference\n(read vs. REF, using Element CIGAR)",
        )
        ax.set_ylabel(
            "Expected allele length difference\n(ALT vs. REF, using TRGT)",
        )
    # g.fig.suptitle(
    #     "Concordance between Element and TRGT allele lengths at homozygous STRs present in all members of a trio",
    #     fontsize=12,
    # )
    g.tight_layout()
    g.savefig("conc.png", dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv")
    p.add_argument("--out")
    args = p.parse_args()
    main()
