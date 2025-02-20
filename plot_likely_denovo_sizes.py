import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss
pd.set_option("display.precision", 8)


plt.rc("font", size=12)

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "CHM13v2"

# define the minimum size of events we want to consider
# (i.e., if MIN_SIZE = 1, don't include interruptions)
MIN_SIZE = 1
MAX_SIZE = 15

mutations = pd.read_csv(f"tr_validation/csv/mutations.{ASSEMBLY}.filtered.csv", dtype={"sample_id": str})

orig = []
# for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.tsv"):
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    orig.append(df)
mutations = pd.concat(orig)

mutations = mutations[mutations["paternal_id"] == "2209"]
mutations = mutations[np.abs(mutations["likely_denovo_size"]).between(MIN_SIZE, MAX_SIZE)]
mutations["Motif type"] = mutations["motif_size"].apply(lambda m: "homopolymer" if m == 1 else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR" if m > 6 else "?")
mutations["Likely DNM size"] = mutations["likely_denovo_size"]
mutations["Parent-of-origin"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else "unknown")

mutations["adj_size"] = mutations["likely_denovo_size"] // mutations["motif_size"]

mutations["Count"] = 1

f, ax = plt.subplots()
sns.pointplot(
    data=mutations[mutations["Motif type"] != "VNTR"],
    x="Likely DNM size",
    y="Count",
    hue="Motif type",
    estimator=lambda x: np.mean(x),
    #errorbar=lambda x: ss.sem(x),
    native_scale=True,
)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("test.png", dpi=200)


# test for expansions vs contractions
n_expansions = mutations.query("likely_denovo_size > 0").shape[0]
n_contractions = mutations.query("likely_denovo_size < 0").shape[0]

print (ss.binomtest(n_expansions, n_expansions + n_contractions, alternative="two-sided"))


val = "likely_denovo_size"

bins = np.arange(mutations[val].min(), mutations[val].max() + 1, 1)

sample_ids = mutations["sample_id"].unique()

f, ax = plt.subplots(figsize=(8, 4))

counts = np.zeros((2, MAX_SIZE * 2, len(sample_ids)))

cmap = sns.color_palette("colorblind", len(mutations["sample_id"].unique()))
mutations = mutations[mutations["Motif type"] != "VNTR"]
for mi, (motif_type, motif_type_df) in enumerate(mutations.groupby("Motif type")):
    for si, (sample, sample_df) in enumerate(motif_type_df.groupby("sample_id")):
        sizes = sample_df[val].values
        hist, edges = np.histogram(sizes, bins=bins)
        counts[mi, :, si] = hist

cmap = sns.color_palette("colorblind", 3)
# cmap = ["cornflowerblue", "goldenrod"]
for mi in (0, 1):

    motif_type_counts = counts[mi]

    fracs = motif_type_counts / np.sum(motif_type_counts, axis=0)[None, :]

    fracs[MAX_SIZE, :] = np.nan
    if mi == 1:
        # don't plot exp/con of 1bp at non-homopolymers
        fracs[MAX_SIZE - 1, :] = np.nan
        fracs[MAX_SIZE + 1, :] = np.nan

    mean_count = np.mean(fracs, axis=1)
    std_count = np.std(fracs, axis=1)
    sem_count = 1.96 * ss.sem(fracs, axis=1)

    x_adj = -0.33 if mi == 0 else 0.33

    # ax.bar(
    #     bins[:-1] + x_adj,
    #     mean_count,
    #     0.33,
    #     lw=1,
    #     ec="w",
    #     zorder=-1,
    #     yerr=sem_count,
    #     color=cmap[mi],
    #     label="homopolymer" if mi == 0 else "non-homopolymer STR",
    # )
    ax.plot(
        bins[:-1],
        mean_count,
        lw=1,
        zorder=-1,
        color=cmap[mi],
    )
    ax.scatter(
        bins[:-1],
        mean_count,
        s=30,
        # ec="w",
        c=cmap[mi],
        label="homopolymer" if mi == 0 else "non-homopolymer STR",
    )
    # ax.errorbar(
    #     bins[:-1],
    #     mean_count,
    #     lw=1,
    #     fmt="none",
    #     capsize=3,
    #     yerr=sem_count, c=cmap[mi]
    # )
    ax.fill_between(
        bins[:-1],
        mean_count - sem_count,
        mean_count + sem_count,
        color=cmap[mi],
        alpha=0.25,
    )
# ax.set_xticks(np.arange(-MAX_SIZE, MAX_SIZE))
xticklabs = list(map(str, np.arange(-MAX_SIZE, MAX_SIZE)))
xticklabs[MAX_SIZE] = ""
# ax.set_xticklabels(xticklabs)
ax.set_ylabel("Fraction of TR DNMs\n(mean +/- 95% CI)")
ax.set_xlabel("TR DNM size (bp)")
ax.set_xlim(-MAX_SIZE - 1, MAX_SIZE)
sns.despine(ax=ax)
ax.legend(shadow=True, fancybox=True, fontsize=10)
f.tight_layout()
f.savefig("sizes.png", dpi=200)
f.savefig("sizes.eps")

mutations["likely_motif_size"] = mutations["likely_denovo_size"] // mutations["motif_size"]
mutations["likely_motif_size_mod"] = mutations["likely_denovo_size"] % mutations["motif_size"]

print (mutations.shape)
mutations = mutations[mutations["likely_motif_size_mod"] == 0]
print (mutations.shape)

motif_sizes = mutations.groupby(["sample_id", "likely_motif_size"]).size().reset_index().rename(columns={0: "count"})
motif_totals = motif_sizes.groupby("sample_id").agg(total=("count", "sum")).reset_index()
motif_sizes = motif_sizes.merge(motif_totals)
motif_sizes["frac"] = motif_sizes["count"] / motif_sizes["total"]
print (motif_sizes)
f, ax = plt.subplots()
sns.pointplot(data=motif_sizes, x="likely_motif_size", y="frac", ax=ax, ls="none")
f.savefig("sizes.motif.png", dpi=200)