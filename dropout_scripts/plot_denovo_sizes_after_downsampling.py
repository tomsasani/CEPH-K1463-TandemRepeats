import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob

def add_likely_de_novo_size(row: pd.Series, use_phase: bool = True):

    child_als = list(map(int, row["child_AL"].split(",")))
    # denovo idx
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]
    non_denovo_al = None

    if not row["trid"].startswith("chrX"):
        non_denovo_idx = 1 - denovo_idx
        non_denovo_al = child_als[non_denovo_idx]

    # figure out the inferred phase of the site, if not unknown
    phase = "unknown"
    if "phase_summary" in row:
        phase = row["phase_summary"].split(":")[0]

    # if phase is unknown, we *could* just use the minimum difference between the
    # de novo AL and any of the four parental ALs, assuming that the smallest
    # difference is the most parsimonious expansion/contraction. however, this is too
    # naive. for example, if the kid's non-denovo AL perfectly matches one of the two parents,
    # then we can assume inheritance of that allele from taht parent. we can call this "fuzzy"
    # phasing.

    # what if the kid perfectly inherited an AL form one parent?
    if use_phase:
        if phase == "unknown":
            parent_als = []
            for parent in ("father", "mother"):
                parent_als.extend(list(map(int, row[f"{parent}_AL"].split(","))))

            abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
            denovo_diffs = np.array([denovo_al - al for al in parent_als])
            min_diff_idx = np.argmin(abs_denovo_diffs)
            # likely original AL
            return denovo_diffs[min_diff_idx]
        else:
            if phase == "dad":
                parent_als = list(map(int, row["father_AL"].split(",")))
            elif phase == "mom":
                parent_als = list(map(int, row["mother_AL"].split(",")))

            abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
            denovo_diffs = np.array([denovo_al - al for al in parent_als])
            min_diff_idx = np.argmin(abs_denovo_diffs)

            # likely original AL
            return denovo_diffs[min_diff_idx]
    else:
        parent_als = []
        for parent in ("father", "mother"):
            parent_als.extend(list(map(int, row[f"{parent}_AL"].split(","))))

        abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
        denovo_diffs = np.array([denovo_al - al for al in parent_als])
        min_diff_idx = np.argmin(abs_denovo_diffs)
        # likely original AL
        return denovo_diffs[min_diff_idx]

def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 10:
            return "unknown"
        else:
            return p.split(":")[0]

CI = 80
N_BOOTS = 1_000

def bootstrap(vals, n: int, bins):
    boots = np.zeros(shape=(n, bins.shape[0] - 1))
    for boot in range(n):
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

ASSEMBLY = "GRCh38"

# which parent to downsample
WHO_TO_DOWNSAMPLE = "kid"

orig = []
for fh in glob.glob(f"tr_validation_new/downsampled/csv/2187.{ASSEMBLY}.*.phased.2gen.tsv"):
# for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    #if len(fh.split("/")[-1].split(".")) > 5: continue
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
    df["downsampled_to"] = ":".join(downsampling) #+ "X"

    #df["downsampled_to"] = fh.split(".")[-3] + "X"
    orig.append(df)
orig = pd.concat(orig)

orig["likely_denovo_size"] = orig.apply(lambda row: add_likely_de_novo_size(row), axis=1)

orig = orig[
        (orig["motif_size"].isin([-1, 1])) | (np.abs(orig["likely_denovo_size"]) >= 2)
    ]
print (orig.groupby("downsampled_to").size())

# orig = pd.read_csv(f"tr_validation/csv/mutations.GRCh")

orig["likely_denovo_size"] = np.abs(orig["likely_denovo_size"])
orig["reference_span"] = orig["end"] - orig["start"]

orig["phase"] = orig["phase_summary"].apply(lambda p: check_phase(p))

size_counts = orig.groupby(["downsampled_to", "likely_denovo_size"]).size().reset_index().rename(columns={0: "count"})
size_totals = orig.groupby("downsampled_to").size().reset_index().rename(columns={0: "total"})
size_merged = size_counts.merge(size_totals)
size_merged["frac"] = size_merged["count"] / size_merged["total"]

ratio_col = "Maternal:paternal depth ratio"
val_col = "likely_denovo_size"

bins = np.arange(1, 10_000, 1)
f, ax = plt.subplots()
for gen, gen_df in orig.groupby("downsampled_to"):
    sizes = np.abs(gen_df[val_col].values)
    edges, mean, lo, hi = bootstrap(sizes, N_BOOTS, bins)
    ax.plot(edges[:-1], mean, label=gen)
    ax.fill_between(edges[:-1], lo, hi, alpha=0.5)
ax.set_xscale("log")
ax.set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
ax.set_xlabel(val_col)
ax.legend(title="Downsampled depth (kid:mom:dad)", fancybox=True, shadow=True)
ax.set_xlabel("Likely de novo size (bp)")
sns.despine(ax=ax)

pref = "/scratch/ucgd/lustre-labs/quinlan/u1006375/sasani-lab-notebook/img/dropout"


f.savefig(f"likely_denovo_size_dist.after.{WHO_TO_DOWNSAMPLE}.png", dpi=200)
f.savefig(f"{pref}/likely_denovo_size_dist.after.{WHO_TO_DOWNSAMPLE}.png", dpi=200)

orig = orig[orig["phase"] != "unknown"]

row_col, line_col = "phase", "downsampled_to"
n_subplots = orig[row_col].nunique()

f, axarr = plt.subplots(n_subplots, figsize=(6, 9))
for i, (a, adf) in enumerate(orig.groupby(row_col)):
    for b, bdf in adf.groupby(line_col):
        sizes = np.abs(bdf[val_col].values)
        edges, mean, lo, hi = bootstrap(sizes, N_BOOTS, bins)
        axarr[i].plot(edges[:-1], mean, label=b)
        axarr[i].fill_between(edges[:-1], lo, hi, alpha=0.5)

    axarr[i].set_xscale("log")
    axarr[i].legend(title="Downsampled depth (kid:mom:dad)", fancybox=True, shadow=True)
    axarr[i].set_title(f"POI = {a}")
    axarr[i].set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
    axarr[i].set_xlabel("Likely de novo size (bp)")


f.tight_layout()
f.savefig(f"likely_denovo_size_dist.phased.{WHO_TO_DOWNSAMPLE}.png", dpi=200)
f.savefig(f"{pref}/likely_denovo_size_dist.phased.{WHO_TO_DOWNSAMPLE}.png", dpi=200)

