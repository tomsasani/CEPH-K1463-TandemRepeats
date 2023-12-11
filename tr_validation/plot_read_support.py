import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import itertools
import numpy as np


def main(args):
    res_df = pd.read_csv(args.csv)

    # res_df = res_df.query("n_with_denovo_allele_strict > 0")
    # res_df = res_df[res_df["Read support"] >= 5]
    res_df_totals = (
        res_df.groupby("region")
        .size()
        .reset_index()
        .rename(columns={0: "Total read support"})
    )
    res_df = res_df.merge(res_df_totals, on="region")

    res_df["diff_to_ref"] = res_df["diff"] - res_df["exp_allele_diff_ref"]
    res_df["exp_allele_diff"] = res_df["alt_al"] - res_df["ref_al"]
    # res_df["Read support"] = res_df["Read support"].astype(np.int32)

    res_df = res_df.query("diff_to_ref != 0")
    res_df = res_df[res_df["Read support"] > 1]

    # res_df["concordance"] = abs(res_df["diff"] - res_df["exp_allele_diff"])
    # res_df["concordance_pct"] = 1 - (res_df["concordance"] / res_df["alt_al"])
    res_df["Is transmitted?"] = res_df["n_with_denovo_allele_strict"] > 0
    print(res_df.groupby("Is transmitted?").size())

    x, y = "exp_allele_diff_alt", "diff"

    g = sns.relplot(
        data=res_df,
        x=x,
        y=y,
        col="In homopolymer?",
        hue="Read support",
        style="Is transmitted?",
        s=100,
    )
    g.map(plt.axline, xy1=(0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)

    titles = [
        r"$\bf{a)}$" + " candidate STR DNMs (not in homopolymers)",
        r"$\bf{b)}$" + " candidate STR DNMs in homopolymers",
    ]
    for i, ax in enumerate(g.axes.flat):
        ax.set_title(titles[i], fontsize=10)
        ax.set_xlabel(
            "Measured allele length difference\n(read vs. REF, using Element CIGAR)",
        )
        ax.set_ylabel(
            "Expected allele length difference\n(ALT vs. REF, using TRGT)",
        )
    g.tight_layout()
    g.savefig(args.out, dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
