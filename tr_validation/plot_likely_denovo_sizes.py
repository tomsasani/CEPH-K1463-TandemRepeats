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
MAX_SIZE = 20

mutations = pd.read_csv(f"tr_validation/csv/mutations.{ASSEMBLY}.filtered.csv", dtype={"sample_id": str})

mutations = mutations[np.abs(mutations["likely_denovo_size"]).between(MIN_SIZE, MAX_SIZE)]
mutations["Motif type"] = mutations["motif_size"].apply(lambda m: "STR" if m <= 6 else "VNTR")
mutations["Likely DNM size"] = mutations["likely_denovo_size"]
mutations["Parent-of-origin"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else "unknown")

mutations["DNM event size (# of motifs)"] = mutations["likely_denovo_size"] // mutations["motif_size"]

val = "likely_denovo_size"

bins = np.arange(mutations[val].min(), mutations[val].max(), 1)

sample_ids = mutations["sample_id"].unique()
f, ax = plt.subplots(figsize=(8, 4))

counts = np.zeros((MAX_SIZE * 2 - 1, len(sample_ids)))

cmap = sns.color_palette("colorblind", len(mutations["sample_id"].unique()))

for si, (sample, sample_df) in enumerate(mutations.groupby("sample_id")):
    sizes = sample_df["likely_denovo_size"].values
    hist, edges = np.histogram(sizes, bins=bins)
    counts[:, si] = hist

fracs = counts / np.sum(counts, axis=0)[None, :]

fracs[MAX_SIZE, :] = np.nan

mean_count = np.mean(fracs, axis=1)
std_count = np.std(fracs, axis=1)
sem_count = 1.96 * ss.sem(fracs, axis=1)

ax.plot(bins[:-1], mean_count, lw=1, c="gainsboro", zorder=-1)
ax.scatter(bins[:-1], mean_count, s=50, ec="w", c="dodgerblue")
ax.errorbar(
    bins[:-1],
    mean_count,
    lw=1,
    fmt="none",
    c="dodgerblue",
    capsize=3,
    yerr=sem_count,
)
ax.set_ylabel("Fraction of TR DNMs\n(mean +/- 95% CI)")
ax.set_xlabel("TR DNM size (bp)")
ax.set_xlim(-MAX_SIZE - 1, MAX_SIZE)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("sizes.png", dpi=200)
f.savefig("sizes.svg")
