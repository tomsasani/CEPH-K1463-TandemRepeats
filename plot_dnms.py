import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import numpy as np

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


mutations = pd.read_csv("CHM13v2.filtered.tsv", sep="\t")

# mutations["is_complex"] = mutations["n_motifs"] > 1
# mutations = mutations[mutations["is_complex"] == False]
mutations["min_motif_size"] = mutations["motif_size"].apply(lambda m: str(m) if int(m) <= 6 else "7+")

counts = mutations.groupby(["sample_id", "min_motif_size"]).size().reset_index().rename(columns={0: "count"})

counts["generation"] = counts["sample_id"].apply(
    lambda s: (
        "G4A"
        if str(s).startswith("2000")
        else "G4B" if str(s).startswith("2001") else "G3"
    ),
)

print (counts)

smps = counts["sample_id"].unique()
smp2idx = dict(zip(smps, range(len(smps))))

ind = np.arange(counts["sample_id"].nunique())

palette = sns.color_palette("deep", 7)

f, ax = plt.subplots(figsize=(7, 5))
bottom = np.zeros_like(ind)
for motif_i, (motif, motif_df) in enumerate(counts.groupby("min_motif_size")):
    motif_counts = np.zeros_like(ind)
    for s, i in smp2idx.items():
        sc = motif_df[motif_df["sample_id"] == s]["count"]
        if sc.shape[0] == 0: continue
        sc = sc.values[0]
        motif_counts[i] = sc
    ax.bar(
        ind,
        motif_counts,
        0.85,
        ec="w",
        lw=1,
        color=palette[motif_i],
        bottom=bottom,
        label=motif,
    )
    bottom += motif_counts
ax.set_xticks(ind)
ax.set_xticklabels(counts["sample_id"].unique(), rotation=45)
ax.set_xlabel("Sample ID")
ax.set_ylabel("# of DNMs")
ax.legend(title="Minimum motif size in locus (bp)", shadow=True)
# g = sns.FacetGrid(data=counts, col="generation", sharex=False)
# g.map(sns.barplot, "sample_id", "count", "min_motif_size")
# for ax in g.axes.flat:
#     for label in ax.get_xticklabels():
#         label.set_rotation(90)
# g.add_legend()
sns.despine(ax=ax)
f.tight_layout()
f.savefig("counts.png", dpi=200)
