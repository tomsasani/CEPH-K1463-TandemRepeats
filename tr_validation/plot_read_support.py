import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

def main(args):

    res_df = pd.read_csv(args.csv)

    res_df["qc"] = abs(res_df["diff"] - res_df["exp_allele_diff"])

    x = res_df["exp_allele_diff"].values
    y = res_df["diff"].values
    s = (res_df["count"].values * 1.5) ** 2

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    ax1.scatter(x, y, alpha=0.5, c="cornflowerblue", s=s, ec="w", lw=1)
    ax1.axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
    ax1.set_title(r"$\bf{a)}$" + " all candidate STR DNMs", fontsize=12)

    ax2.scatter(x, y, alpha=0.5, c="cornflowerblue", s=s, ec="w", lw=1)
    ax2.axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
    ax2.set_xlim(-7, 10)
    ax2.set_title(r"$\bf{b)}$" + " zoomed to show subset of shorter STRs", fontsize=12)
    ax2.set_ylim(-7, 10)

    for ax in (ax1, ax2):
        sns.despine(ax=ax, top=True, right=True)
        ax.set_xlabel("Measured allele length difference\n(read vs. REF)")

    ax1.set_ylabel("Expected allele length difference\n(REF vs. ALT, using AL field)")
    f.suptitle("How well do Element reads corroborate TRGT DNM calls?", fontsize=14)

    f.tight_layout()
    f.savefig(args.out, dpi=200)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv")
    p.add_argument("--out")
    args = p.parse_args()