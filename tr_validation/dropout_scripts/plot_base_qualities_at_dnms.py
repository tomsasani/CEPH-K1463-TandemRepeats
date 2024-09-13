import pysam
import pandas as pd
import glob
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns


CI = 80

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

def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 10:
            return "unknown"
        else:
            return p.split(":")[0]

ASSEMBLY = "GRCh38"
BAM_PREF = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/hifi-bams/{ASSEMBLY}"

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str,})

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

orig = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str, "maternal_id": str})
    orig.append(df)
orig = pd.concat(orig)
orig = orig.query("likely_denovo_size != 0")
orig["likely_denovo_size"] = np.abs(orig["likely_denovo_size"])
orig["phase"] = orig["phase_summary"].apply(lambda p: check_phase(p))


res = []
for sample, sample_df in orig.groupby("sample_id"):
    sample = SMP2ALT[sample]
    mom, dad = sample_df["maternal_id"].unique()[0], sample_df["paternal_id"].unique()[0]

    for parent, parent_name in zip((mom, dad), ("mom", "dad")):
        parent_alt = SMP2ALT[parent]
        # read in bam
        bamfile = pysam.AlignmentFile(f"{BAM_PREF}/{parent_alt}.GRCh38.haplotagged.bam", "rb")

        counted = 0
        for i,row in tqdm.tqdm(sample_df.iterrows()):
            # if counted > 5: break
            chrom, start, end = row["#chrom"], row["start"], row["end"]
            row_dict = row.to_dict()
            bqs, lens, mqs, nps, rqs = [], [], [], [], []
            for read in bamfile.fetch(chrom, start, end):
                if read.is_unmapped: continue
                bqs.append(np.mean(read.query_qualities))
                lens.append(len(read.query_qualities))
                mqs.append(read.mapping_quality)
                nps.append(read.get_tag("np"))
                rqs.append(read.get_tag("rq"))
            res.append(
                {
                    "sample": sample,
                    "parent": parent_alt,
                    "trid": row["trid"],
                    "read": read.query_name,
                    "phase": row["phase"],
                    "mean_bq_per_read": np.mean(bqs),
                    "mean_read_len": np.mean(lens),
                    "mean_mq": np.mean(mqs),
                    "mean_np": np.mean(nps),
                    "mean_rq": np.mean(rqs),
                }
            )
            counted += 1


res_df = pd.DataFrame(res)

res_df_tidy = res_df.melt(
    id_vars=["sample", "parent", "trid", "read", "phase"],
    var_name="metric",
    value_name="value",
)

res_df_tidy["generation"] = res_df_tidy["parent"].apply(
    lambda s: (
        "G4A"
        if s in ("NA12879", "200080")
        else (
            "G4B"
            if s in ("NA12886", "200100")
            else "G3" if s in ("NA12877", "NA12878") else "G2"
        )
    )
)

res_df_tidy = res_df_tidy[res_df_tidy["generation"].isin(["G3", "G4A", "G4B"])]

metrics = res_df_tidy["metric"].unique()
phases = res_df_tidy["phase"].unique()

metric2i = dict(zip(metrics, range(len(metrics))))
phase2i = dict(zip(phases, range(len(phases))))

# g = sns.FacetGrid(data=res_df_tidy, row="metric", sharex=False, aspect=2)
# g.map(sns.histplot, "value", None, "phase", multiple="layer")

# g = sns.displot(data=res_df_tidy[~res_df_tidy["generation"].str.startswith("G2")], row="metric", col="phase", hue="generation", x="value", kind="ecdf", facet_kws={"sharex": False})
# g.savefig("t.png")

print (res_df_tidy.groupby(["metric", "phase", "parent"]).size())

f, axarr = plt.subplots(len(metrics), len(phases), figsize=(14, 20), sharey=True)

for metric, metric_df in res_df_tidy.groupby("metric"):
    vals = metric_df["value"].values
    bins = np.linspace(min(vals), max(vals), num=100)
    for phase, phase_df in metric_df.groupby("phase"):
        mi, pi = metric2i[metric], phase2i[phase]
        print (metric, phase, phase_df.shape)
        for gen, gen_df in phase_df.groupby("parent"):
            if gen.startswith("G2"): continue
            vals = gen_df["value"].values
            edges, mean, lo, hi = bootstrap(vals, 1_000, bins)
            axarr[mi, pi].plot(edges[:-1], mean, label=gen)
            axarr[mi, pi].fill_between(edges[:-1], lo, hi, alpha=0.5)
        ylim = (0, 1)
        if metric == "mean_mq":
            ylim = (0, 0.25)
        axarr[mi, pi].set_ylabel("Cumulative fraction of DNMs\n(mean +/- bootstrap 95% CI)")
        axarr[mi, pi].set_xlabel(metric)
        axarr[mi, pi].legend(title="parent")
        axarr[mi, pi].set_title(f"POI = {phase}")
        # if metric == "mean_mq":
        #     axarr[mi, pi].set_yscale("log")
        sns.despine(ax=axarr[mi, pi])
f.tight_layout()
f.savefig("test.png", dpi=200)
