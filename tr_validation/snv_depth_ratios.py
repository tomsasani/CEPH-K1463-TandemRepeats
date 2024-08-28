import pandas as pd
import glob
import matplotlib.pyplot as plt
import utils
import seaborn as sns
import numpy as np
from cyvcf2 import VCF
import tqdm
from collections import Counter
import glob

plt.rc("font", size=10)

MIN_DEPTH = 10
ASSEMBLY = "GRCh38"

# read in G3 and G2 DNMs
dnms = pd.read_csv("tr_validation/data/Supplementary_Table11 - SNVs.tsv", sep="\t", skiprows=1)
dnms = dnms[dnms["variant origin"] == "germline"]
dnms = dnms[dnms["inheritance"] != "cannot_determine"]
dnms = dnms.dropna(subset=["GRCh38 position"])

mapping = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str, "paternal_id": str, "maternal_id": str},
)

smp2dad, smp2mom = (
    dict(zip(mapping["alt_sample_id"], mapping["paternal_id"])),
    dict(zip(mapping["alt_sample_id"], mapping["maternal_id"])),
)
orig2alt = dict(zip(mapping["sample_id"], mapping["alt_sample_id"]))

vcf_fh = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.{ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz"
vcf = VCF(vcf_fh, gts012=True)

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

snvs = []
snv_pos_col = "T2T-CHM13v2.0 position" if ASSEMBLY == "CHM13v2" else "GRCh38 position"

for i, row in tqdm.tqdm(dnms.iterrows()):
    sample = row["sample"]
    dad, mom = orig2alt[smp2dad[sample]], orig2alt[smp2mom[sample]]
    si, di, mi = smp2idx[sample], smp2idx[dad], smp2idx[mom]

    chrom, start = row["chromosome"], int(row[snv_pos_col])
    end = start + 1
    region = f"{chrom}:{start}-{end}"
    for v in vcf(region):

        rd = v.gt_ref_depths
        ad = v.gt_alt_depths
        td = rd + ad
        if td[si] < MIN_DEPTH or td[di] < MIN_DEPTH or td[mi] < MIN_DEPTH: continue
        gts = v.gt_types
        if not (gts[si] == 1 and (gts[di] + gts[mi]) == 0): continue

        snvs.append(
            {
                "region": region,
                "sample": sample,
                "kid_depth": td[si],
                "mom_depth": td[mi],
                "dad_depth": td[di],
                "poi": row["inheritance"],
            }
        )


snv_df = pd.DataFrame(snvs)
snv_df["ratio"] = snv_df["mom_depth"] / snv_df["dad_depth"]
snv_df["mutation_type"] = "SNV"

trs = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    trs.append(df)
tr_df = pd.concat(trs).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
tr_df = tr_df[np.abs(tr_df["likely_denovo_size"]) > 0]


tr_df = tr_df[(tr_df["mom_dp"] >= MIN_DEPTH) & (tr_df["dad_dp"] >= MIN_DEPTH) & (tr_df["child_coverage"] >= MIN_DEPTH)]
tr_df["ratio"] = tr_df["mom_dp"] / tr_df["dad_dp"]
tr_df["mutation_type"] = "TR"

tr_df["region"] = tr_df["trid"].apply(lambda t: t.split("_")[0] + ":" + t.split("_")[1] + "-" + t.split("_")[2])
tr_df["chrom"] = tr_df["trid"].apply(lambda t: t.split("_")[0])

print (tr_df[tr_df["region"] == "chr1:2217040-2217960"][["sample_id", "trid", "index", "child_AL", "mother_AL", "father_AL", "ratio", "mom_dp", "dad_dp", "likely_denovo_size"]])

tr_df["poi"] = tr_df["phase_summary"].apply(
    lambda p: (
        "maternal"
        if p.split(":")[0] == "mom"
        else "paternal" if p.split(":")[0] == "dad" else "cannot_determine"
    )
)
tr_df = tr_df[tr_df["poi"] != "cannot_determine"]
tr_df["sample"] = tr_df["alt_sample_id"]

tr_df["reference_len"] = tr_df.apply(lambda row: row["end"] - row["start"], axis=1)


combined = pd.concat([snv_df[["sample", "region", "ratio", "poi", "mutation_type"]], tr_df[["sample", "region", "ratio", "poi", "mutation_type"]]])
combined["generation"] = combined["sample"].apply(
    lambda s: (
        "G2"
        if str(s) in ("NA12877", "NA12878")
        else (
            "G4A"
            if str(s).startswith("20008")
            else "G4B" if str(s).startswith("2001") else "G3"
        )
    )
)

meds = []

f, axarr = plt.subplots(2, figsize=(8, 8))
bins=np.arange(0, 3, 0.05)

gen2color = {"G2": "k", "G3": "dodgerblue", "G4A": "firebrick", "G4B": "olive"}
mut2idx = {"TR": 0, "SNV": 1}
seen = []
for (smp, generation, mut_type), smp_df in combined.groupby(["sample", "generation", "mutation_type"]):

    if smp_df.shape[0] == 0: continue

    bad = smp_df.query("ratio < 0.5")
    if generation == "G4" and bad.shape[0] > 0:
        print (bad)

    pois = smp_df["poi"].values
    poi_dict = dict(Counter(pois).most_common())
    if "maternal" not in poi_dict:
        maternal_frac = 0
    else:
        maternal_frac = poi_dict["maternal"] / pois.shape[0]
    ratios = smp_df["ratio"].values
    hist, edges = np.histogram(ratios, bins=bins)
    hist_fracs = hist / np.sum(hist)

    axarr[mut2idx[mut_type]].plot(np.cumsum(hist_fracs), c=gen2color[generation], label=generation if generation not in seen else None)
    if mut_type == "TR":
        seen.append(generation)

    meds.append(
        {
            "sample": smp,
            "generation": generation,
            "mutation_type": mut_type,
            "Median maternal:paternal depth ratio across DNMs": np.median(ratios),
            "Maternal POI fraction": maternal_frac,
        }
    )

    axarr[mut2idx[mut_type]].axvline(20, ls=":", c="gainsboro", zorder=-1)
    axarr[mut2idx[mut_type]].set_xticks(np.arange(bins.shape[0])[::10])
    axarr[mut2idx[mut_type]].set_xticklabels(bins[::10])

for lab, idx in mut2idx.items():
    axarr[idx].set_title(f"{lab} DNMs")
    sns.despine(ax=axarr[idx])
    axarr[idx].legend(fancybox=True)
axarr[1].set_xlabel("Ratio of maternal to paternal depth at DNMs")
axarr[1].set_ylabel("Cumulative fraction of DNMs")
f.tight_layout()
f.savefig('o.png', dpi=200)


med_df = pd.DataFrame(meds)

f, ax = plt.subplots()
g = sns.relplot(
    data=med_df,
    x="Median maternal:paternal depth ratio across DNMs",
    y="Maternal POI fraction",
    hue="generation",
    col="mutation_type",
    # facet_kws={"sharex": False, "sharey": False},
    s=75,
)
g.savefig("frac.png", dpi=200)


# remove X chromosome
combined = combined[~(combined["region"].str.startswith("chrX")) | (combined["region"].str.startswith("chr19"))]
# combined = combined[combined["region"].str.startswith("chr1")]
combined = combined[combined["generation"] == "G3"]

# genomic distribution
genome = pd.read_csv("tr_validation/data/hg38.genome", sep="\t")
chrom2len = dict(zip(genome["chrom"], genome["size"]))
combined["chrom"] = combined["region"].apply(lambda r: r.split(":")[0])

genome_dist = []

for (chrom, mutation_type), chrom_df in combined.groupby(["chrom", "mutation_type"]):
    genome_bins = np.arange(0, chrom2len[chrom], 1_000_000)
    start = chrom_df["region"].apply(lambda r: int(r.split(":")[1].split("-")[0]))
    end = chrom_df["region"].apply(lambda r: int(r.split(":")[1].split("-")[1]))
    mid = (end + start) / 2
    bin_i = np.digitize(mid, genome_bins)
    
    for b in bin_i:
        if b - 1 == -1: print ("AH")
        genome_dist.append({"chromosome": chrom, "mutation_type": mutation_type, "binned_start": genome_bins[b - 1]})

genome_dist_df = pd.DataFrame(genome_dist)#.groupby(["chrom", "mutation_type"]).
g = sns.displot(
    genome_dist_df,
    x="binned_start",
    row="chromosome",
    col="mutation_type",
    aspect=4,
    kind="hist",
    # stat="proportion",
)
g.savefig("dist.png")
