import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import itertools
import numpy as np


def jitter(values, j):
    return values + np.random.normal(
        j,
        0.15,
        values.shape,
    )


def main(args):
    res_df = pd.read_csv(args.csv)

    # res_df = res_df.query("n_with_denovo_allele_strict > 0")
    #res_df = res_df[res_df["Read support"] >= 5]
    res_df_totals = res_df.groupby("region").size().reset_index().rename(columns={0: "Total read support"})
    res_df = res_df.merge(res_df_totals, on="region")

    res_df["diff_to_ref"] = res_df["diff"] - res_df["exp_allele_diff_ref"]
    res_df["exp_allele_diff"] = res_df["alt_al"] - res_df["ref_al"]
    

    res_df = res_df.query("diff_to_ref != 0")
    res_df = res_df[res_df["Read support"] > 1]
    print (res_df.query("diff_to_ref < -5"))

    # res_df["concordance"] = abs(res_df["diff"] - res_df["exp_allele_diff"])
    # res_df["concordance_pct"] = 1 - (res_df["concordance"] / res_df["alt_al"])
    res_df["Is transmitted?"] = res_df["n_with_denovo_allele_strict"] > 0

    f, axarr = plt.subplots(1, 2, figsize=(12, 6))

    x, y = "exp_allele_diff_alt", "diff"
    s = 100

    sns.scatterplot(
        data=res_df,
        x=x,
        y=y,
        hue="Read support",
        style="Is transmitted?",
        s=s,
        ax=axarr[0],
    )
    axarr[0].axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
    axarr[0].set_title(r"$\bf{a)}$" + " all candidate STR DNMs", fontsize=12)

    # sns.scatterplot(
    #     data=res_df,
    #     x=jitter(res_df["diff_to_ref"], 0),
    #     y="exp_allele_diff",
    #     hue="Read support",
    #     style="transmitted",
    #     s=s,
    #     ax=axarr[0, 1],
    # )
    # axarr[0, 1].axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
    # axarr[0, 1].set_xlim(-7, 10)
    # axarr[0, 1].set_title(
    #     r"$\bf{b)}$" + " zoomed to show subset of shorter STRs", fontsize=12
    # )
    # axarr[0, 1].set_ylim(-7, 10)

    res_df["motif_length"] = res_df["motif"].apply(lambda m: len(m))

    # sns.barplot(data=res_df, x="motif_length", y="concordance_pct", ax=axarr[1, 1])
    # axarr[1, 1].set_title(
    #     "Concordance between TRGT and Element AL\n(as pct. of total allele length)"
    # )

    sns.scatterplot(
        data=res_df[res_df["in_homopolymer"] == True],
        x=x,
        y=y,
        hue="Read support",
        style="Is transmitted?",
        s=s,
        ax=axarr[1],
    )
    axarr[1].axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
    axarr[1].set_title(
        r"$\bf{b)}$" + " all candidate STR DNMs in homopolymers", fontsize=12
    )

    for x in (0, 1):#itertools.product([0, 1], repeat=2):
        sns.despine(ax=axarr[x], top=True, right=True)
        #if (x, y) != (1, 1):
        axarr[x].set_xlabel(
            "Measured allele length difference\n(read vs. REF, using CIGAR)"
        )
        axarr[x].set_ylabel(
            "Expected allele length difference\n(ALT vs. REF, using AL)"
        )

    # ax1.set_ylabel("Expected allele length difference\n(REF vs. ALT, using AL field)")
    f.suptitle("How well do Element reads corroborate TRGT DNM calls?", fontsize=14)

    f.tight_layout()
    f.savefig(args.out, dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
