import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tqdm
import glob

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))


ASSEMBLY = "CHM13v2"

VCF_FH = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.chm13.deepvariant.glnexus.phased.vcf.gz"
AWS_PREF = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25/"
vcf = VCF(VCF_FH, gts012=True)
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

KID = "2216"
MOM, DAD = "2209", "2188"
KID_VCF, MOM_VCF, DAD_VCF = [AWS_PREF + "hiphase/" + sample + "_" + SMP2SUFF[sample] + f"_{ASSEMBLY.lower()}" + "_50bp_merge.sorted.phased.vcf.gz" for sample in (KID, MOM, DAD)]
KID_VCF, MOM_VCF, DAD_VCF = [VCF(vcf_fh, gts012=True) for vcf_fh in (KID_VCF, MOM_VCF, DAD_VCF)]


SLOP = 10_000

mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"}).drop_duplicates("trid")
mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else "unknown")
mutations = mutations[mutations["maternal_id"] == 2216]

mutations["ratio"] = mutations["mom_dp"] / mutations["dad_dp"]

mutations = mutations.query("ratio < 2")

res = []

for row_i, (_, row) in tqdm.tqdm(enumerate(mutations.iterrows())):
    trid = row["trid"]
    chrom, start, end, method = trid.split("_")
    start, end = int(start) - SLOP, int(end) + SLOP
    if start < 0: start = 1
    region = f"{chrom}:{start}-{end}"
    trid_mid = (int(end) + int(start)) / 2

    genotypes = []

    for li, (vcf_, label) in enumerate(zip((KID_VCF, MOM_VCF, DAD_VCF), ("kid", "mom", "dad"))):
        for v in vcf_(region):
            gts = v.gt_types
            genotypes.append({"pos": v.POS, "gt": gts[0], "sample": label})
    genotypes_df = pd.DataFrame(genotypes)

    
    for pos, pos_df in genotypes_df.groupby("pos"):
        if pos_df["gt"].nunique() != 3: continue
        mom_gt, dad_gt, kid_gt = [pos_df[pos_df["sample"] == s]["gt"].values[0] for s in ("mom", "dad", "kid")]
        if abs(mom_gt - dad_gt) != 2: continue
        obligate = kid_gt == 1

        row_dict = row.to_dict()

        row_dict.update({"pos": pos, "obligate": obligate})
        
        res.append(row_dict)
    

res_df = pd.DataFrame(res)


group_cols = ["trid", "ratio"]

totals = res_df.groupby(group_cols).size().reset_index().rename(columns={0: "total"})
obligates = res_df.groupby(group_cols).agg(obligate=("obligate", "sum")).reset_index()

new = totals.merge(obligates)
new["frac"] = new["obligate"] / new["total"]


f, ax = plt.subplots()
ax.scatter(new["ratio"], new["frac"], alpha=0.5)
f.savefig("obligate.png", dpi=200)
            
