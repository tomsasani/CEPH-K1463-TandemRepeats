import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import upsetplot

CI = 95
N_BOOTS = 1_000

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
            else "G3" if s == "2209" else "NA12877" if s == "2281" else "NA12878"
        )
    )
)
print (orig)

f, ax = plt.subplots(figsize=(10, 5))
# ax2 = ax.twinx()

MAX_MOTIF_SIZE = 6

gen2color = dict(zip(orig["generation"].unique(), sns.color_palette("colorblind", orig["generation"].nunique())))

bins = np.arange(orig["motif_size"].min(), orig["motif_size"].max(), 1)
n_cats = orig["generation"].nunique()


ind = np.arange(MAX_MOTIF_SIZE) - 0.33
cur_x_adj = 0
labeled = []
for downsample, downsample_df in orig.groupby("generation"):
    # if downsample in ("60X", "70X"): continue
    gen = downsample_df["generation"].unique()[0]
    if gen in ("G4A", "G4B"): continue
    motif_counts = downsample_df["motif_size"].values
    edges, mean, lo, hi = bootstrap(motif_counts, n=N_BOOTS, bins=bins)
    offsets = np.array([mean - lo, hi - mean])
    ax.bar(
        ind + cur_x_adj,
        mean[:MAX_MOTIF_SIZE],
        0.8 / n_cats,
        yerr=offsets[:, :MAX_MOTIF_SIZE],
        label=gen if gen not in labeled else None,
        color=gen2color[gen],
    )

    motif_size, count = np.unique(motif_counts[motif_counts <= MAX_MOTIF_SIZE], return_counts=True)
    # ax2.scatter((motif_size - 1 - 0.33) + cur_x_adj, count, c="gainsboro", ec="k", lw=1)
    cur_x_adj += 0.8 / n_cats
    labeled.append(gen)


ax.set_xticks(np.arange(MAX_MOTIF_SIZE))
ax.set_xticklabels(np.arange(MAX_MOTIF_SIZE) + 1)
ax.legend(title="Generation", frameon=False, fancybox=False, shadow=False)
ax.set_title("Homopolymer DNMs are enriched in G4")
ax.set_ylabel("Fraction of DNMs (+/- 95% CI)")
ax.set_xlabel("Motif size (bp)")
# ax2.set_ylabel("Count of DNMs")
# ax2.set_yscale("log")
sns.despine(ax=ax)
# sns.despine(ax=ax2, right=False)

pref = "/scratch/ucgd/lustre-labs/quinlan/u1006375/sasani-lab-notebook/img/dropout"

f.savefig(f"motif_counts.parent.before.all.png", dpi=200)
f.savefig(f"{pref}/motif_counts.parent.before.all.png", dpi=200)
