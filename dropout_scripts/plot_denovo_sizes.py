import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import tqdm


def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 2:
            return "unknown"
        else:
            return p.split(":")[0]

CI = 95
N_BOOTS = 1_000

def bootstrap(vals, n: int, bins):
    boots = np.zeros(shape=(n, bins.shape[0] - 1))
    for boot in tqdm.tqdm(range(n)):
        # random sample
        boot_sizes = np.random.choice(vals, size=vals.shape[0], replace=True)
        hist, edges = np.histogram(boot_sizes, bins=bins)
        hist_fracs = hist / np.sum(hist)
        boots[boot, :] = np.cumsum(hist_fracs)

    mean = np.mean(boots, axis=0)
    lo_bound = (100 - CI) / 2
    lo = np.percentile(boots, lo_bound, axis=0)
    hi = np.percentile(boots, 100 - lo_bound, axis=0)

    return edges, mean, lo, hi

ASSEMBLY = "CHM13v2"

orig = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    if len(fh.split("/")[-1].split(".")) > 5: continue
    print (fh)
    # df["downsampled_to"] = "_".join(fh.split(".")[-5:-2]) #+ "X"
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    df["downsampled_to"] = fh.split(".")[-3] + "X"
    orig.append(df)
orig = pd.concat(orig)

orig = orig[
        (orig["motif_size"].isin([-1, 1])) | (np.abs(orig["likely_denovo_size"]) >= 2)
    ]
print (orig.groupby("sample_id").size())

# orig = pd.read_csv(f"tr_validation/csv/mutations.GRCh")

orig["likely_denovo_size"] = np.abs(orig["likely_denovo_size"])
orig["reference_span"] = orig["end"] - orig["start"]

orig["generation"] = orig["paternal_id"].apply(
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

orig["phase"] = orig["phase_summary"].apply(lambda p: check_phase(p))

# orig = orig[orig["phase"] == "mom"]

size_counts = orig.groupby(["generation", "likely_denovo_size"]).size().reset_index().rename(columns={0: "count"})
size_totals = orig.groupby("generation").size().reset_index().rename(columns={0: "total"})
size_merged = size_counts.merge(size_totals)
size_merged["frac"] = size_merged["count"] / size_merged["total"]

ratio_col = "Maternal:paternal depth ratio"
val_col = "likely_denovo_size"



cmap = sns.color_palette("colorblind", orig["generation"].nunique())

bins = np.arange(1, 10_000, 1)
f, ax = plt.subplots()
for gi, (gen, gen_df) in enumerate(orig.groupby("generation")):
    # if gen == "G3": continue
    sizes = np.abs(gen_df[val_col].values)
    edges, mean, lo, hi = bootstrap(sizes, N_BOOTS, bins)
    ax.plot(edges[:-1], mean, label=gen, c=cmap[gi], lw=2)
    ax.fill_between(edges[:-1], lo, hi, alpha=0.5, color=cmap[gi])
ax.set_xscale("log")
ax.set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
ax.set_xlabel(val_col)
ax.legend(title="Generation", fancybox=True, shadow=True)
ax.set_xlabel("Likely de novo size (bp)")
sns.despine(ax=ax)

pref = "/scratch/ucgd/lustre-labs/quinlan/u1006375/sasani-lab-notebook/img/dropout"


f.savefig(f"likely_denovo_size_dist.before.png", dpi=200)
f.savefig(f"{pref}/likely_denovo_size_dist.before.png", dpi=200)

orig = orig[orig["phase"] != "unknown"]

orig = orig[orig["generation"].isin(["G2A","G2B", "G3"])]

row_col, line_col = "generation", "phase"
n_subplots = orig[row_col].nunique()

f, axarr = plt.subplots(n_subplots, figsize=(6, 9))
for i, (a, adf) in enumerate(orig.groupby(row_col)):
    for b, bdf in adf.groupby(line_col):
        sizes = np.abs(bdf[val_col].values)
        edges, mean, lo, hi = bootstrap(sizes, N_BOOTS, bins)
        axarr[i].plot(edges[:-1], mean, label=b)
        axarr[i].fill_between(edges[:-1], lo, hi, alpha=0.5)

    axarr[i].set_xscale("log")
    axarr[i].legend(title=line_col)
    axarr[i].set_title(f"POI = {a}")
    axarr[i].set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
    axarr[i].set_xlabel(val_col)
f.tight_layout()
f.savefig("likely_denovo_size_dist.phased.png", dpi=200)

