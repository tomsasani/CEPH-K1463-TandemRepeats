import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List
import glob

plt.rc("font", size=10)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", dtype={"sample_id": str, "paternal_id": str}, sep="\t")

mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4A" if s.startswith("2000") else "G4B" if s.startswith("2001") else "G3")

mutations = mutations[mutations["phase"] != "unknown"]
mutations["TR_type"] = mutations["max_motiflen"].apply(lambda m: "homopolymer" if m == 1 else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR" if m > 6 else "?")

print (mutations[mutations["sample_id"].isin(["200101", "200084"])].groupby(["sample_id", "phase"]).size())


bs = []
for trial in range(100):
    _mutations = mutations.sample(frac=1, replace=True)

    counts = _mutations.groupby(["sample_id", "generation", "phase"]).size().reset_index().rename(columns={0: "count"})
    totals = counts.groupby(["sample_id"]).agg(total=("count", "sum")).reset_index()

    counts = counts.merge(totals)
    counts = counts[counts["phase"] == "dad"]
    counts["frac"] = counts["count"] / counts["total"]
    counts["trial"] = trial

    bs.append(counts)
bs = pd.concat(bs)

print (bs.groupby("generation").agg(avg=("frac", "mean")))
f, ax = plt.subplots(figsize=(12, 4))
sns.barplot(data=bs, x="sample_id", y="frac", hue="generation", ax=ax)
ax.set_ylabel("Paternal DNM fraction")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("phase.png", dpi=200)
