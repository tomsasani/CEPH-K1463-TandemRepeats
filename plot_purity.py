import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
import matplotlib.patches as patches

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

rng = np.random.default_rng(42)

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"sample_id": str})

for c in mutations.columns:
   if np.any(mutations[c].isna()): print (c)

mutations = mutations[mutations["haplotype_in_parent_consensus"] != "UNK"]
# mutations = mutations[mutations["phase"] != "unknown"]

# remove non-complex
# mutations = mutations[mutations["n_motifs"] == 1]


# mutations = mutations[(mutations["precursor_AP"] > 0) & (mutations["untransmitted_AP"] > 0)]

mutations["precursor_AL"] = mutations["precursor_sequence_in_parent"].apply(lambda s: len(s))
mutations["untransmitted_AL"] = mutations["untransmitted_sequence_in_parent"].apply(lambda s: len(s))


mutations["direction"] = mutations["likely_denovo_size"].apply(
    lambda d: "expansion" if d > 0 else "contraction" if d < 0 else "interruption"
)


f, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharey=True)

palette = sns.color_palette("colorblind", 2)

vals = []

for i, (val, ax) in enumerate(zip(("AL", "AP"), (ax1, ax2))):

    tidy = mutations.melt(
        id_vars=["sample_id", "trid", "genotype", "likely_denovo_size"],
        value_vars=[f"precursor_{val}", f"untransmitted_{val}"],
        var_name="allele_status",
        value_name="size",
    )

    tidy["is_precursor"] = tidy["allele_status"].apply(lambda s: 1 if "precursor" in s else 0)

    # wilcoxon signed rank test on allele purities. alternative hypothesis
    # of "left" indicates that the precursor APs are on the whoel lower than the untransmitted
    stat, p = ss.wilcoxon(
        mutations[f"precursor_{val}"].values,
        y=mutations[f"untransmitted_{val}"].values,
        alternative="greater",
    )

    # diffs = tidy.groupby(["trid", "sample_id", "genotype"]).agg(diff=("value", lambda ap: np.diff(ap)[0]))["diff"]
    a, b = mutations[f"precursor_{val}"].values, mutations[f"untransmitted_{val}"].values
    diffs = a - b


    bin_edges = np.linspace(min(diffs), max(diffs), 50)

    v_hist, v_edges = np.histogram(diffs, bins=bin_edges)
    v_hist_frac = v_hist / np.sum(v_hist)
    v_hist_cumsum = np.cumsum(v_hist_frac)

    ind = np.arange(v_edges[1:].shape[0])

    ax.bar(v_edges[1:], v_hist, v_edges[1:] - v_edges[:-1], ec="w", color="k", lw=1,  zorder=1)
    rect1 = patches.Rectangle((0, 0), max(v_edges), max(v_hist), facecolor='green', zorder=-1, alpha=0.1)
    rect2 = patches.Rectangle((0, 0), min(v_edges), max(v_hist), facecolor='red', zorder=-1, alpha=0.1)

    # Add the patch to the axes
    ax.add_patch(rect1)
    ax.add_patch(rect2)

    # ax.set_title("p = {:0.2e}".format(p))
    ax.set_xlabel("Difference between allele lengths..." if val == "AL" else "...and allele purities")
    # ax.set_xticks(ind[::20])
    # ax.set_xticklabels(v_edges[1::20])
    # print (min(diffs))
    sns.despine(ax=ax)
    # ax.set_xscale("log")
    ax.set_yscale("log")

ax1.set_ylabel(f"Number of TR DNMs")

f.tight_layout()
f.savefig(f"length_purity.png", dpi=200)
