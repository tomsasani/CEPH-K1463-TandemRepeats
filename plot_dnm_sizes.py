import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import scipy.stats as ss
import tqdm


CI = 95
N_BOOTS = 1_000

plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def get_motif_types(row):
    if row["n_motifs"] > 1:
        if row["max_motiflen"] > 6 and row["min_motiflen"] <= 6:
            return "complex"
        else:
            if row["min_motiflen"] == 1:
                return "homopolymer"
            elif row["min_motiflen"] > 1 and row["max_motiflen"] <= 6:
                return "non-homopolymer STR"
            else:
                return "VNTR"
    else:
        if row["min_motiflen"] == 1:
            return "homopolymer"
        elif 1 < row["min_motiflen"] <= 6:
            return "non-homopolymer STR"
        else:
            return "VNTR"

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

def get_size_in_motifs(row):
    size = row["likely_denovo_size"]
    motif_size = row["min_motiflen"]

    if size % motif_size != 0:
        return np.nan
    else:
        return size // motif_size

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv("CHM13v2.filtered.tsv", sep="\t", dtype={"paternal_id": str})

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

mutations["is_expansion"] = mutations["likely_denovo_size"].apply(lambda s: s > 0)


n_exp = mutations.groupby("is_expansion").size().to_dict()
print (ss.binomtest(n_exp[True], n_exp[True] + n_exp[False]))

mutations["TR type"] = mutations.apply(
    lambda row: get_motif_types(row),
    axis=1,
)

f, ax = plt.subplots(figsize=(9, 6), sharex=True)

MAX_DNM_SIZE = 10

cmap = sns.color_palette("colorblind", mutations["generation"].nunique())

n_bootstraps = 100

colors = sns.color_palette("deep", 2)

for ki, (kind, kind_df) in enumerate(
    mutations[
        mutations["TR type"].isin(("homopolymer", "non-homopolymer STR"))
    ].groupby("TR type")
):

    bootstrapped = np.zeros((n_bootstraps, MAX_DNM_SIZE * 2 + 1))
    subset = kind_df[(kind_df["likely_denovo_size"] >= -MAX_DNM_SIZE) & (kind_df["likely_denovo_size"] <= MAX_DNM_SIZE)]

    if kind == "non-homopolymer STR": 
        assert subset.query("likely_denovo_size == 1").shape[0] == 0
        bootstrapped[1 + MAX_DNM_SIZE] = np.nan
        bootstrapped[-1 + MAX_DNM_SIZE] = np.nan

    subset["is_expansion"] = subset["likely_denovo_size"].apply(lambda s: s > 0)

    n_exp = subset.groupby("is_expansion").size().to_dict()
    print (kind, ss.binomtest(n_exp[True], n_exp[True] + n_exp[False]))

    for bi in tqdm.tqdm(range(n_bootstraps)):
        resampled = subset.sample(frac=1, replace=True)
        counts = resampled.groupby("likely_denovo_size").size().reset_index().rename(columns={0: "count"})
        counts["frac"] = counts["count"] / counts["count"].sum()
        for motif_size in np.arange(-MAX_DNM_SIZE, MAX_DNM_SIZE + 1):
            motif_frac = counts[counts["likely_denovo_size"] == motif_size]
            if motif_frac.shape[0] == 0: continue
            bootstrapped[bi, motif_size + MAX_DNM_SIZE] = motif_frac["frac"].values[0]
    ind = np.arange(-MAX_DNM_SIZE, MAX_DNM_SIZE + 1)
    mean = np.mean(bootstrapped, axis=0)
    sem = ss.sem(bootstrapped, axis=0)
    mean[MAX_DNM_SIZE] = np.nan
    sem[MAX_DNM_SIZE] = np.nan
    # ax.plot(ind, mean, color=colors[ki], lw=1, label=kind)
    ax.errorbar(ind, mean, yerr=1.96 * sem, linestyle="none", color=colors[ki])
    ax.scatter(ind, mean, c=colors[ki],  label=kind)
    ax.plot(ind, mean, c=colors[ki])
    # ax.fill_between(ind, mean - 1.96 * sem, mean + 1.96 * sem, color=colors[ki], alpha=0.25)
    # sns.lineplot(data=fracs, x="likely_denovo_size", y="frac")
ax.set_xlabel(
    r"$\longleftarrow$"
    + "contraction"
    + "\t\t"
    + "expansion"
    + r"$\longrightarrow$"
    + "\n"
    + r"$\bf{Inferred\ size\ of\ DNM\ (bp)}$",
    size=12,
)
ax.set_xticks(np.arange(-MAX_DNM_SIZE, MAX_DNM_SIZE + 1))
ax.axvline(0, ls=":", c="gainsboro", lw=2)
ax.set_ylabel("Fraction of DNMs")
ax.legend(title="TR type")
sns.despine(ax=ax)

f.tight_layout()
f.savefig("sizes.png", dpi=200)


# do the same but for number of motifs
mutations = mutations[mutations["min_motiflen"] == mutations["max_motiflen"]]
print (mutations.shape)
mutations["size_in_motifs"] = mutations.apply(lambda row: get_size_in_motifs(row), axis=1)
mutations.dropna(inplace=True)
print (mutations.shape)

res = []
for bs in range(100):
    mutations_resampled = mutations.sample(frac=1, replace=True)
    counts = mutations_resampled.groupby("size_in_motifs").size().reset_index().rename(columns={0: "count"})
    counts["t"] = bs
    counts["frac"] = counts["count"] / mutations_resampled.shape[0]
    res.append(counts)
res = pd.concat(res)
res["size_in_motifs"] = res["size_in_motifs"].astype(int)


f, ax = plt.subplots()
sns.pointplot(data=res, x="size_in_motifs", y="frac", ax=ax,  linestyle="none")
f.savefig("sizes.motifs.png", dpi=200)