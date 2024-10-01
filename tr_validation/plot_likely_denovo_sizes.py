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
ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

# define the minimum size of events we want to consider
# (i.e., if MIN_SIZE = 1, don't include interruptions)
MIN_SIZE = 1
MAX_SIZE = 25

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == 2209]
samples = mutations["sample_id"].unique()

mutations = mutations[np.abs(mutations["likely_denovo_size"]).between(MIN_SIZE, MAX_SIZE)]
mutations["Motif type"] = mutations["simple_motif_size"]
mutations["Likely DNM size"] = mutations["likely_denovo_size"]
mutations["Parent-of-origin"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else "unknown")

bins = np.arange(mutations["likely_denovo_size"].min(), mutations["likely_denovo_size"].max(), 1)

print (mutations.groupby(["motif_size", "likely_denovo_size"]).size())

f, ax = plt.subplots(figsize=(5, 12))
#sns.stripplot(data=mutations[mutations["motif_size"].between(1, 10)], x="motif_size", y="likely_denovo_size", alpha=0.5, ax=ax)
sns.stripplot(data=mutations[mutations["likely_denovo_size"].between(-5, 5)], x="likely_denovo_size", y="motif_size", alpha=0.5, ax=ax)
f.savefig("test.png")

mutations = mutations.groupby(["sample_id",  "likely_denovo_size"]).size().reset_index().rename(columns={0: "count"})
totals = mutations.groupby("sample_id").size().reset_index().rename(columns={0: "total"})
mutations = mutations.merge(totals)
mutations["frac"] = mutations["count"] / mutations["total"]

f, ax = plt.subplots()

sns.pointplot(data=mutations, x = "likely_denovo_size", y="frac", ax=ax)
# for i, (motif_class, motif_df) in enumerate(mutations.groupby("simple_motif_size")):
#     bottom = np.zeros(bins.shape[0])

#     for sample, sample_df in motif_df.groupby("sample_id"):
#         sizes = sample_df["likely_denovo_size"].values
#         hist, edges = np.histogram(sizes, bins=bins)
#         axarr[i].bar(bins[:-1], hist, 1, bottom=bottom[:-1], ec="w", lw=0.5, label=sample)
#         bottom[:-1] += hist
#     axarr[i].set_title(motif_class)
#     axarr[i].set_ylabel("# of TR DNMs")
#     if i == 0:
#         axarr[i].legend(shadow=True, fancybox=True, title="Sample ID")
#     if i == 2:
#         axarr[i].set_xlabel("TR DNM size (bp)")
#     sns.despine(ax=axarr[i])
# # g = sns.FacetGrid(data=mutations, col="Motif type")
# # g.map(sns.histplot, "Likely DNM size", bins=bins, hue="sample_id", lw=1, edgecolor="w", stat="proportion")
# #ax.set_xlabel("Likely DNM size (bp)")
# f.tight_layout()
f.savefig("sizes.png", dpi=200)



#### PLOT MOTIF SIZE DISTRO ###

# read in all per-sample DNM files
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "motif_size": int})
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(samples)]
denoms = denoms.groupby(["sample_id", "motif_size"]).agg(denominator=("denominator", "sum"))

f, ax = plt.subplots()
sns.pointplot(data=denoms.query("motif_size > 0 and motif_size < 20"), x="motif_size", y="denominator", ax=ax)
f.savefig("motif.denominator.png", dpi=200)

# for i, (motif_class, motif_df) in enumerate(mutations.groupby("motif_size")):
#     bottom = np.zeros(bins.shape[0])

#     for sample, sample_df in motif_df.groupby("sample_id"):
#         sizes = sample_df["likely_denovo_size"].values
#         hist, edges = np.histogram(sizes, bins=bins)
#         axarr[i].bar(bins[:-1], hist, 1, bottom=bottom[:-1], ec="w", lw=0.5, label=sample)
#         bottom[:-1] += hist
#     axarr[i].set_title(motif_class)
#     axarr[i].set_ylabel("# of TR DNMs")
#     if i == 0:
#         axarr[i].legend(shadow=True, fancybox=True, title="Sample ID")
#     if i == 2:
#         axarr[i].set_xlabel("TR DNM size (bp)")
#     sns.despine(ax=axarr[i])
# # g = sns.FacetGrid(data=mutations, col="Motif type")
# # g.map(sns.histplot, "Likely DNM size", bins=bins, hue="sample_id", lw=1, edgecolor="w", stat="proportion")
# #ax.set_xlabel("Likely DNM size (bp)")
# f.tight_layout()
# f.savefig("sizes.png", dpi=200)