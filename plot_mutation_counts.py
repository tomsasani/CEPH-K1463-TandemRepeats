import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List
import glob

plt.rc("font", size=16)

def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 10:
            return "unknown"
        else:
            return p.split(":")[0]

def calculate_phase_counts(df: pd.DataFrame, group_cols: List[str] = ["sample_id", "phase"]):
    phase_counts = df.groupby(group_cols).size().reset_index().rename(columns={0: "count"})
    totals = (
        phase_counts.groupby("alt_sample_id")
        .agg({"count": "sum"})
        .reset_index()
        .rename(columns={"count": "total"})
    )

    phase_counts = phase_counts.merge(totals)
    phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]
    return phase_counts

ASSEMBLY = "CHM13v2"


orig = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    orig.append(df)
orig = pd.concat(orig)

orig["Minimum motif size"] = orig["motif_size"].apply(lambda m: str(m) if m <= 6 else "7+")


f, ax = plt.subplots()

phase_counts = calculate_phase_counts(
    orig,
    group_cols=["alt_sample_id", "paternal_id", "Minimum motif size"],
)

phase_counts["generation"] = phase_counts["paternal_id"].apply(
    lambda s: (
        "G4A"
        if s == "200080"
        else (
            "G4B"
            if s == "2189"
            else "G3" if s == "2209" else "G2A" if s == "2281" else "G2B"
        )
    )
)

phase_counts = phase_counts[phase_counts["generation"].isin(["G3", "G4A", "G4B"])]

order = phase_counts.sort_values(["generation", "alt_sample_id"])["alt_sample_id"].unique()
order2idx = dict(zip(list(map(str, order)), range(len(order))))
phase_counts["order"] = phase_counts["alt_sample_id"].apply(lambda s: order2idx[str(s)])
print (phase_counts)

f, ax = plt.subplots(figsize=(7, 10))
bottom = np.zeros(len(order2idx))
ind = np.arange(len(order2idx))

cmap = sns.color_palette("colorblind", 7)
for i, (motif, motif_df) in enumerate(phase_counts.groupby("Minimum motif size")):
    
    motif_df_sorted = motif_df.sort_values("order")

    count_arr = np.zeros(phase_counts["alt_sample_id"].nunique())
    for o, oi in order2idx.items():
        o_count = motif_df_sorted[motif_df_sorted["alt_sample_id"].astype(str) == o]
        if o_count.shape[0] == 0:
            continue
        else:
            count_arr[oi] = o_count["count"].values[0]

    color = cmap[i]
    ax.barh(ind, count_arr, 0.8, left=bottom, label=motif, color=color, lw=1, ec="w")
    bottom += count_arr


ax.legend(title="Motif size (bp)", frameon=True, shadow=True, fancybox=True)
ax.set_yticks(ind)
ax.set_yticklabels([o for o,i in order2idx.items()])

ax.set_xlabel("Number of TR DNMs")
ax.set_ylabel("Sample ID")
sns.despine(ax=ax)
sns.set_style("ticks")

f.tight_layout()
f.savefig("o.png", dpi=200)
