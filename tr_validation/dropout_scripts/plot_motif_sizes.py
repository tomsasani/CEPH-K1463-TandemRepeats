import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob


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
for fh in glob.glob(f"tr_validation/downsampled/csv/2187.{ASSEMBLY}.50.50.*.prefiltered.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    df["downsampled_to"] = "_".join(fh.split(".")[-5:-2]) #+ "X"
    # df = df[df["downsampled_to"].isin(["2187", "200087", "200106"])]
    orig.append(df)
orig = pd.concat(orig)
print (orig["downsampled_to"].unique())
# orig = orig.query("likely_denovo_size != 0")

f, ax = plt.subplots(figsize=(10, 5))

bins = np.arange(orig["motif_size"].min(), orig["motif_size"].max(), 1)
n_cats = orig["downsampled_to"].nunique()
ind = np.arange(10) - 0.33
cur_x_adj = 0
for downsample, downsample_df in orig.groupby("downsampled_to"):
    motif_counts = downsample_df["motif_size"].values
    edges, mean, lo, hi = bootstrap(motif_counts, n=N_BOOTS, bins=bins)
    offsets = np.array([mean - lo, hi - mean])
    ax.bar(ind + cur_x_adj, mean[:10], 0.8 / n_cats, yerr=offsets[:, :10], label=downsample)
    cur_x_adj += 0.8 / n_cats
    # ax.plot(edges[:-1], mean, label=downsample)
    # ax.fill_between(edges[:-1], lo, hi, alpha=0.5)


# motif_counts = orig.groupby(["downsampled_to", "motif_size"]).size().reset_index().rename(columns={0: "count"})
# generation_counts = orig.groupby("downsampled_to").size().reset_index().rename(columns={0: "total"})
# motif_counts = motif_counts.merge(generation_counts)
# motif_counts["frac"] = motif_counts["count"] / motif_counts["total"]
# sns.barplot(data=motif_counts.query("motif_size <= 10"), x="motif_size", y="frac", hue="downsampled_to", ax=ax)
ax.set_xticks(np.arange(10))
ax.set_xticklabels(np.arange(10) + 1)
ax.legend(title="Downsampled depth in trio (child is 2187)", frameon=False, fancybox=False, shadow=False)
ax.set_title("Homopolymer DNMs are enriched when sequencing depth is low")
ax.set_ylabel("Fraction of DNMs (+/- 95% CI)")
ax.set_xlabel("Motif size (bp)")
sns.despine(ax=ax)
f.savefig("motif_counts.parent.png", dpi=200)

# orig["phase"] = orig["phase_summary"].apply(lambda p: check_phase(p))

# size_counts = orig.groupby(["generation", "likely_denovo_size"]).size().reset_index().rename(columns={0: "count"})
# size_totals = orig.groupby("generation").size().reset_index().rename(columns={0: "total"})
# size_merged = size_counts.merge(size_totals)
# size_merged["frac"] = size_merged["count"] / size_merged["total"]

# ratio_col = "Maternal:paternal depth ratio"
# val_col = "likely_denovo_size"

# bins = np.arange(1, 10_000, 1)
# f, ax = plt.subplots()
# for gen, gen_df in orig.groupby("generation"):
#     sizes = np.abs(gen_df[val_col].values)
#     edges, mean, lo, hi = bootstrap(sizes, N_BOOTS, bins)
#     ax.plot(edges[:-1], mean, label=gen)
#     ax.fill_between(edges[:-1], lo, hi, alpha=0.5)
# ax.set_xscale("log")
# ax.set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
# ax.set_xlabel(val_col)
# ax.legend(title="generation")
# f.savefig("boop.png", dpi=200)