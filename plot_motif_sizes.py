import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

CI = 95
N_BOOTS = 1_000

plt.rc("font", size=12)


def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 10:
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

orig = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    orig.append(df)
orig = pd.concat(orig)

orig["phase"] = orig["phase_summary"].apply(lambda p: check_phase(p))

orig["generation"] = orig["paternal_id"].apply(
    lambda s: (
        "G4A"
        if s == "200080"
        else (
            "G4B"
            if s == "2189"
            else "G3" if s == "2209" else "G2"
        )
    )
)


# calculate total number of DNMs per individual
totals = orig.groupby("sample_id").size().reset_index().rename(columns={0: "total"})

orig = orig.merge(totals)
orig["count"] = 1
orig["frac"] = 1 / orig["total"]

g = sns.FacetGrid(data=orig, row="phase")
g.map(sns.barplot, "motif_size", "count", "generation", estimator="sum")
g.savefig("p.png")

f, ax = plt.subplots(figsize=(8, 4))
# ax2 = ax.twinx()


orig["motif_size"] = orig["motif_size"].apply(lambda m: m if m <= 6 else 7)

gen2color = dict(zip(orig["generation"].unique(), sns.color_palette("colorblind", orig["generation"].nunique())))

bins = np.arange(1, 9)
n_cats = orig["generation"].nunique()


MAX_MOTIF_SIZE = 7

ind = np.arange(0, MAX_MOTIF_SIZE) - 0.33
cur_x_adj = 0
labeled = []
for generation, generation_df in orig.groupby("generation"):
    motif_counts = generation_df["motif_size"].values
    print (motif_counts)
    edges, mean, lo, hi = bootstrap(motif_counts, n=N_BOOTS, bins=bins)
    print (edges, mean)
    offsets = np.array([mean - lo, hi - mean])
    ax.bar(
        ind + cur_x_adj,
        mean,
        0.8 / n_cats,
        yerr=offsets,
        label=generation if generation not in labeled else None,
        color=gen2color[generation],
        ec="w",
        lw=1,
    )

    motif_size, count = np.unique(motif_counts[motif_counts <= MAX_MOTIF_SIZE], return_counts=True)
    # ax2.scatter((motif_size - 1 - 0.33) + cur_x_adj, count, c="gainsboro", ec="k", lw=1)
    cur_x_adj += 0.8 / n_cats
    labeled.append(generation)


ax.set_xticks(np.arange(MAX_MOTIF_SIZE))
ax.set_xticklabels(list(map(lambda m: str(m) if m <= 6 else "7+", np.arange(MAX_MOTIF_SIZE) + 1)))
ax.legend(title="Generation", frameon=True, fancybox=True, shadow=True)
ax.set_ylabel("Fraction of DNMs (+/- 95% CI)")
ax.set_xlabel("Motif size (bp)")
# ax2.set_ylabel("Count of DNMs")
# ax2.set_yscale("log")
sns.despine(ax=ax)
f.tight_layout()
# sns.despine(ax=ax2, right=False)

f.savefig(f"motif_counts.png", dpi=200)
