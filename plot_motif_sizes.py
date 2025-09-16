import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

CI = 95
N_BOOTS = 1_000

plt.rc("font", size=12)

def check_phase(p):
    if p == "unknown": return p
    else:
        support = float(p.split(":")[1])
        if support < 0.75:
            return "unknown"
        else:
            return p.split(":")[0]


def bootstrap(vals, n: int, bins):
    boots = np.zeros(shape=(n, bins.shape[0] - 1))
    for boot in range(n):
        # random sample
        boot_sizes = np.random.choice(vals, size=vals.shape[0], replace=True)
        hist, edges = np.histogram(boot_sizes, bins=bins)
        hist_fracs = hist / np.sum(hist)
        boots[boot, :] = hist_fracs

    mean = np.mean(boots, axis=0)
    lo_bound = (100 - CI) / 2
    lo = np.percentile(boots, lo_bound, axis=0)
    hi = np.percentile(boots, 100 - lo_bound, axis=0)

    return edges, mean, lo, hi

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv("CHM13v2.filtered.tsv", sep="\t", dtype={"sample_id": str, "paternal_id": str})

denoms = []
for fh in glob.glob(f"csv/denominators/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    denoms.append(df)
denoms = pd.concat(denoms)
denoms["min_motif_size"] = denoms["motif_size"].apply(lambda m: str(m) if int(m) <= 6 else "7+")

denom_counts = denoms.groupby(["sample_id", "min_motif_size"]).agg(count=("denominator", "sum")).reset_index()
denom_totals = denom_counts.groupby(["sample_id"]).agg(total=("count", "sum")).reset_index()
denom_counts = denom_counts.merge(denom_totals)
denom_counts["frac"] = denom_counts["count"] / denom_counts["total"]

mutations["min_motif_size"] = mutations["motif_size"].apply(lambda m: str(m) if int(m) <= 6 else "7+")

mutations["phase"] = mutations["phase_consensus"].apply(lambda p: check_phase(p))

mutations["generation"] = mutations["paternal_id"].apply(
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

f, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()
ax1.sharey(ax2)

counts = (
    mutations.groupby(["sample_id", "min_motif_size", "generation"])
    .size()
    .reset_index()
    .rename(
        columns={0: "count"},
    )
)
totals = counts.groupby(["sample_id", "generation"]).agg(total=("count", "sum")).reset_index()
counts = counts.merge(totals)
counts["frac"] = counts["count"] / counts["total"]

sns.boxplot(
    data=counts,
    x="min_motif_size",
    y="frac",
    dodge=True,
    color="white",
    ax=ax1,
    fliersize=0,
)
sns.stripplot(
    data=counts,
    x="min_motif_size",
    y="frac",
    palette="deep",
    ax=ax1,
)

ax1.set_xlabel("Minimum motif size in locus (bp)")
ax1.set_ylabel("Fraction of TR DNMs")
ax2.set_ylabel("Fraction of all TR loci")
sns.despine(ax=ax1, right=False, top=True)
sns.despine(ax=ax2, top=True)
f.tight_layout()
f.savefig("motifs.png", dpi=200)
