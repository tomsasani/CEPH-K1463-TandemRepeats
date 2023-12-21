import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import itertools
import numpy as np


def main(args):
    res_df = pd.read_csv(args.csv)

    res_df_totals = (
        res_df.groupby("region")
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    res_df = res_df.merge(res_df_totals, on="region")

    # require at least 10 total spanning reads
    res_df = res_df[res_df["Total read support"] >= 10]
    # require at least 2 reads for a given observation of allele length
    res_df["frac"] = res_df["Read support"] / res_df["Total read support"]

    # require an observation to be supported by at least 5% of reads
    # res_df = res_df[res_df["frac"] >= 0.05]

    # add a homopolymer column if it's not in there
    if "is_homopolymer" not in res_df:
        res_df["is_homopolymer"] = res_df["motif"].apply(lambda m: len(m) == 1)

    
    # ignore sites where the ref and alt ALs are identical
    res_df = res_df[res_df["exp_allele_diff_ref"] != res_df["exp_allele_diff_alt"]]

    res_df["Matches REF"] = res_df["diff"] == res_df["exp_allele_diff_ref"] 
    res_df["Matches ALT"] = res_df["diff"] == res_df["exp_allele_diff_alt"] 

    res_df_parental_ev = res_df[res_df["sample"] != "kid"].groupby(["region", "sample", "Matches ALT"]).size().reset_index().rename(columns={0: "count"})
    res_df_parental_ev = res_df_parental_ev[res_df_parental_ev["Matches ALT"] == True]
    
    to_filter = res_df_parental_ev["region"].to_list()


    res_df["Evidence for DNM allele in parents?"] = res_df["region"].apply(lambda r: r in to_filter)

    res_df["diff_to_ref"] = res_df["diff"] - res_df["exp_allele_diff_ref"]

    res_df = res_df.query("diff_to_ref != 0")

    if "n_with_denovo_allele_strict" in res_df.columns:
        res_df["Is transmitted?"] = res_df["n_with_denovo_allele_strict"] > 0
    else:
        res_df["Is transmitted?"] = "unknown"

    x, y = "exp_allele_diff_alt", "diff"

    g = sns.relplot(
        data=res_df[(res_df["sample"] == "kid") & (res_df["Matches REF"] == False)],
        x=x,
        y=y,
        col="is_homopolymer",
        # row="Evidence for DNM allele in parents?",
        hue="frac",
        style="Is transmitted?",
        s=100,
    )
    g.map(plt.axline, xy1=(0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)


    # titles = [
    #     r"$\bf{a)}$" + " candidate STR DNMs (not in homopolymers)",
    #     r"$\bf{b)}$" + " candidate STR DNMs in homopolymers",
    # ]
    # for i, ax in enumerate(g.axes.flat):
    #     ax.set_title(titles[i], fontsize=10)
    #     ax.set_xlabel(
    #         "Measured allele length difference\n(read vs. REF, using Element CIGAR)",
    #     )
    #     ax.set_ylabel(
    #         "Expected allele length difference\n(ALT vs. REF, using TRGT)",
    #     )
    g.tight_layout()
    g.savefig(args.out, dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
