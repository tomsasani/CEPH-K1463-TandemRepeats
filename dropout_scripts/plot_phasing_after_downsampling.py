import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import upsetplot

CI = 95
N_BOOTS = 1_000

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

ASSEMBLY = "GRCh38"

orig = []
for fh in glob.glob(f"tr_validation/downsampled/csv/2187.{ASSEMBLY}.*.phased.2gen.tsv"):
    #if len(fh.split("/")[-1].split(".")) > 5: continue
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    mom_depth, dad_depth = list(map(int, fh.split(".")[-5:-3]))
    who_is_downsampled = "Mom" if mom_depth < dad_depth else "Dad"
    min_d, max_d = min([mom_depth, dad_depth]), max([mom_depth, dad_depth])
    frac = min_d / max_d
    df["who_is_downsampled"] = who_is_downsampled
    df["Downsampling relative to other parent"] = frac#f"{frac}% of other parent"
    # df["downsampled_to"] = fh.split(".")[-3] + "X"
    
    # df = df[df["downsampled_to"].isin(["2187", "200087", "200106"])]
    orig.append(df)
orig = pd.concat(orig)


orig["phase_int"] = orig["phase_summary"].apply(lambda p: 1 if p.split(":")[0] == "dad" else 0 if p.split(":")[0] == "mom" else -1)
orig = orig[orig["phase_int"] != -1]

f, ax = plt.subplots()
sns.barplot(data=orig, x="who_is_downsampled", hue="Downsampling relative to other parent", y="phase_int", ax=ax)
ax.set_xlabel("Which parent was downsampled?")
ax.set_ylabel("Fraction of DNMs assigned to dad")
ax.set_title("Imbalanced parental depths have a major impact\non parent-of-origin inference at TRs")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("o.png", dpi=200)




f, ax = plt.subplots(figsize=(10, 5))
ax2 = ax.twinx()

MAX_MOTIF_SIZE = 6

bins = np.arange(orig["motif_size"].min(), orig["motif_size"].max(), 1)
n_cats = orig["downsampled_to"].nunique()
ind = np.arange(MAX_MOTIF_SIZE) - 0.33
cur_x_adj = 0
for downsample, downsample_df in orig.groupby("downsampled_to"):
    motif_counts = downsample_df["motif_size"].values
    edges, mean, lo, hi = bootstrap(motif_counts, n=N_BOOTS, bins=bins)
    offsets = np.array([mean - lo, hi - mean])
    ax.bar(ind + cur_x_adj, mean[:MAX_MOTIF_SIZE], 0.8 / n_cats, yerr=offsets[:, :MAX_MOTIF_SIZE], label=f"{downsample} (n = {motif_counts.shape[0]})")

    motif_size, count = np.unique(motif_counts[motif_counts <= MAX_MOTIF_SIZE], return_counts=True)
    ax2.scatter((motif_size - 1 - 0.33) + cur_x_adj, count, ec="k", lw=1)
    cur_x_adj += 0.8 / n_cats


ax.set_xticks(np.arange(MAX_MOTIF_SIZE))
ax.set_xticklabels(np.arange(MAX_MOTIF_SIZE) + 1)
ax.legend(title="Downsampled depth in trio (child is 2187)", frameon=False, fancybox=False, shadow=False)
ax.set_title("Homopolymer DNMs are enriched when sequencing depth is low")
ax.set_ylabel("Fraction of DNMs (+/- 95% CI)")
ax.set_xlabel("Motif size (bp)")
ax2.set_ylabel("Count of DNMs")
# ax2.set_yscale("log")
sns.despine(ax=ax)
sns.despine(ax=ax2, right=False)

f.savefig("motif_counts.parent.png", dpi=200)


# what fraction of DNMs are removed as we sequence to higher and higher depths
f, ax = plt.subplots()

lenient_dnms = orig[(orig["downsampled_to"] == "10X") & (orig["motif_size"] <= MAX_MOTIF_SIZE)]["trid"].to_list()

not_lenient = orig[(orig["downsampled_to"] != "10X") & (orig["motif_size"] <= MAX_MOTIF_SIZE)]
not_lenient["overlap_type"] = not_lenient["trid"].apply(lambda t: "overlap" if t in lenient_dnms else "private")

grouped = not_lenient.groupby(["downsampled_to", "motif_size", "overlap_type"]).size().reset_index().rename(columns={0: "count"})
print (grouped)
g = sns.FacetGrid(data=grouped, row="downsampled_to", aspect=3, sharey=False)
g.map(sns.barplot, "motif_size", "count", "overlap_type", hue_order=["overlap", "private"])
g.add_legend()
g.savefig("overlap.png")
