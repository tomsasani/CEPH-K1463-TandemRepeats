import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import upsetplot

CI = 80
N_BOOTS = 1_000

def bootstrap(vals, n: int, bins):
    boots = np.zeros(shape=(n, bins.shape[0] - 1))
    for boot in range(n):
        # random sample
        boot_sizes = np.random.choice(vals, size=vals.shape[0], replace=True)
        hist, edges = np.histogram(boot_sizes, bins=bins)
        hist_fracs = hist / np.sum(hist)
        boots[boot, :] = hist

    mean = np.mean(boots, axis=0)
    lo_bound = (100 - CI) / 2
    lo = np.percentile(boots, lo_bound, axis=0)
    hi = np.percentile(boots, 100 - lo_bound, axis=0)

    return edges, mean, lo, hi

ASSEMBLY = "GRCh38"
WHO_TO_DOWNSAMPLE = "all"

orig = []
for fh in glob.glob(f"tr_validation/downsampled/csv/2187.{ASSEMBLY}.*.phased.2gen.tsv"):
    downsampling = fh.split("/")[-1].split(".")[2:5]
    kid, mom, dad = downsampling

    if WHO_TO_DOWNSAMPLE == "kid":
        if not (mom == "50" and dad == "50"): continue
    elif WHO_TO_DOWNSAMPLE == "mom":
        if not (dad == "50" and kid == "50"): continue
    elif WHO_TO_DOWNSAMPLE == "dad":
        if not (mom == "50" and kid == "50"): continue
    elif WHO_TO_DOWNSAMPLE == "all":
        if not (mom == kid == dad): continue
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    # df["downsampled_to"] = "_".join(fh.split(".")[-5:-2]) #+ "X"
    df["downsampled_to"] = ":".join(downsampling)
    
    orig.append(df)
orig = pd.concat(orig)

orig["phase"] = orig["phase_summary"].apply(lambda p: p.split(":")[0])
# orig = orig[orig["phase"] == "dad"]

f, ax = plt.subplots(figsize=(10, 5))
# ax2 = ax.twinx()

MAX_MOTIF_SIZE = 6

bins = np.arange(orig["motif_size"].min(), orig["motif_size"].max(), 1)
n_cats = orig["downsampled_to"].nunique()
ind = np.arange(MAX_MOTIF_SIZE) - 0.33
cur_x_adj = 0
for downsample, downsample_df in orig.groupby("downsampled_to"):
    # if downsample in ("60X", "70X"): continue
    motif_counts = downsample_df["motif_size"].values
    edges, mean, lo, hi = bootstrap(motif_counts, n=N_BOOTS, bins=bins)
    offsets = np.array([mean - lo, hi - mean])
    ax.bar(ind + cur_x_adj, mean[:MAX_MOTIF_SIZE], 0.8 / n_cats, yerr=offsets[:, :MAX_MOTIF_SIZE], label=f"{downsample} (n = {motif_counts.shape[0]})")

    motif_size, count = np.unique(motif_counts[motif_counts <= MAX_MOTIF_SIZE], return_counts=True)
    # ax2.scatter((motif_size - 1 - 0.33) + cur_x_adj, count, c="gainsboro", ec="k", lw=1)
    cur_x_adj += 0.8 / n_cats


ax.set_xticks(np.arange(MAX_MOTIF_SIZE))
ax.set_xticklabels(np.arange(MAX_MOTIF_SIZE) + 1)
ax.legend(title="Downsampled depth (kid:mom:dad)", frameon=False, fancybox=False, shadow=False)
# ax.set_title("Homopolymer DNMs are enriched when sequencing depth is low")
ax.set_ylabel("Fraction of DNMs (+/- 95% CI)")
ax.set_xlabel("Motif size (bp)")
# ax2.set_ylabel("Count of DNMs")
# ax2.set_yscale("log")
sns.despine(ax=ax)
# sns.despine(ax=ax2, right=False)

pref = "/scratch/ucgd/lustre-labs/quinlan/u1006375/sasani-lab-notebook/img/dropout"
f.savefig(f"{pref}/motif_counts.after.{WHO_TO_DOWNSAMPLE}.counts.png", dpi=200)
f.savefig(f"motif_counts.after.{WHO_TO_DOWNSAMPLE}.png", dpi=200)
