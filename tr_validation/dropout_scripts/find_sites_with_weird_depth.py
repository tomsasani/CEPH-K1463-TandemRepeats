import pandas as pd
import glob
import matplotlib.pyplot as plt
import utils
import seaborn as sns
import numpy as np

plt.rc("font", size=13)

mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.CHM13v2.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
mutations["assembly"] = "T2T"

# mutations = mutations[mutations["paternal_id"] == 200080]

print (mutations.groupby("sample_id").size())

# remove sites with weird depth imbalance between parents
mutations["dp_ratio"] = mutations["mom_dp"] / mutations["dad_dp"]

f, ax = plt.subplots()
bins=np.arange(0, 3, 0.05)
for smp, smp_df in mutations.groupby("sample_id"):
    ratios = smp_df["dp_ratio"].values
    hist, edges = np.histogram(ratios, bins=bins)
    hist_fracs = hist / np.sum(hist)
    ax.plot(np.cumsum(hist_fracs), c="grey")
ax.axvline(20, ls=":", c="k", zorder=-1)
ax.set_xticks(np.arange(bins.shape[0])[::10])
ax.set_xticklabels(bins[::10])
ax.set_xlabel("Ratio of maternal to paternal depth at DNMs")
ax.set_ylabel("Cumulative fraction of DNMs")
sns.despine(ax=ax)
f.tight_layout()
f.savefig('o.png', dpi=200)