from cyvcf2 import VCF
import cyvcf2
import pandas as pd
import csv
import tqdm
from typing import List, Dict
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import gzip
from schema import DeNovoSchema
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def catalog_informative_sites(
    vcf: VCF,
    region: str,
):
    res = []
    for v in vcf(region):
        if v.var_type != "snp": continue
        if v.gt_quals[0] < 20: continue
        rd, ad = v.gt_ref_depths, v.gt_alt_depths
        if rd[0] + ad[0] < 10: continue

        gts = v.genotypes
        hap_0, hap_1, is_phased = gts[0]
        if not is_phased: continue
        if hap_0 + hap_1 != 1: continue

        out_dict = {
            "chrom": v.CHROM,
            "pos": v.POS,
            "gt": "|".join((str(hap_0), str(hap_1))),
            "haplotype_A_allele": hap_0,
            "haplotype_B_allele": hap_1,
        }

        res.append(out_dict)

    return pd.DataFrame(res)

def catalog_informative_sites_parents(
    vcf: VCF,
    region: str,
    parents: List[str],
    smp2idx: Dict[str, int],
):
    
    idxs = np.array([smp2idx[p] for p in parents])
    res = []
    for v in vcf(region):
        if v.var_type != "snp": continue
        if np.any(v.gt_quals[idxs] < 20): continue
        rd, ad = v.gt_ref_depths, v.gt_alt_depths
        td = ad + rd
        if np.any(td[idxs] < 10): continue

        gts = v.gt_types[idxs]
        if np.any(gts == 3): continue
        if gts[0] == gts[1]: continue

        # we consider the informative parent to be the one
        # that donated the ALT allele
        if gts[0] > gts[1]:
            inf_parent, oth_parent = parents

        else:
            inf_parent, oth_parent = parents[::-1]

        out_dict = {
                "chrom": v.CHROM,
                "pos": v.POS,
                "alt_parent": inf_parent,
                "ref_parent": oth_parent,
            }
    
        res.append(out_dict)

    return pd.DataFrame(res)

JOINT_SNP_VCF = VCF("/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/Illumina-dragen_v4.2.4/dragen_v4.2.4_joint_genotyped/palladium.v4.2.4.grc38.multisample.joint_genotyped.hard-filtered.vcf.gz", gts012=True)
SMP2IDX = dict(zip(JOINT_SNP_VCF.samples, range(len(JOINT_SNP_VCF.samples))))

MOM_VCF = VCF("tr_validation/data/hiphase/NA12878.GRCh38.deepvariant.glnexus.phased.vcf.gz", gts012=True)
DAD_VCF = VCF("tr_validation/data/hiphase/NA12877.GRCh38.deepvariant.glnexus.phased.vcf.gz", gts012=True)
KID_VCF = VCF("tr_validation/data/hiphase/NA12886.GRCh38.deepvariant.glnexus.phased.vcf.gz", gts012=True)
KID_STR_VCF = VCF("tr_validation/data/hiphase/2189_SF_GRCh38_50bp_merge.sorted.phased.vcf.gz", gts012=True)

KID_DNMS = "tr_validation/data/df_transmission_C0_T0.parquet.gzip"
DNM_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/"
KID_DNMS = f"{DNM_PREF}/2189_SF_GRCh38_50bp_merge_trgtdn.csv.gz"

KID_PHASE_BLOCKS = "tr_validation/data/hiphase/NA12886.GRCh38.hiphase.blocks.tsv"

CHROMS = [f"chr{c}" for c in range(1, 23)]
#CHROMS = ["chr21"]

informative_sites = []

with open(KID_PHASE_BLOCKS, "r") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    header = None
    for i,l in tqdm.tqdm(enumerate(csvf)):
        if header is None:
            header = l
            continue
        d = dict(zip(header, l))
        chrom, start, end = d["chrom"], d["start"], d["end"]
        if chrom not in CHROMS:
            continue
        region = f"{chrom}:{start}-{end}"

        # inf_sites_mom = catalog_informative_sites(
        #     MOM_VCF,
        #     region,
        #     individual="mom",
        #     expected_gts=[1, 2],
        # )
        # inf_sites_dad = catalog_informative_sites(
        #     DAD_VCF,
        #     region,
        #     individual="dad",
        #     expected_gts=[1, 2],
        # )
        inf_sites = catalog_informative_sites_parents(JOINT_SNP_VCF, region, ["2209", "2188"], SMP2IDX)
        het_sites_kid = catalog_informative_sites(
            KID_VCF,
            region,
        )

        if any([idf.shape[0] == 0 for idf in (inf_sites, het_sites_kid)]): continue

        #inf_sites = inf_sites_mom.merge(inf_sites_dad, on=["chrom", "pos"])

        # figure out the parental origin of the REF or ALT alleles in the child.
        # since we only consider sites at which one parent is HOM_ALT and the other is HET,
        # we can always know for sure which parent donated either allele.
        #inf_sites = inf_sites[inf_sites["mom_gt"] != inf_sites["dad_gt"]]
        # if dad's genotype is 1, we know that the origin of the REF allele is dad. if 
        # dad's genotype is 2, we know the REF allele is from mom. similarly, if dad's 
        # genotype is 1, we know the ALT came from mom.
        #inf_sites["ref_origin"] = inf_sites["dad_gt"].apply(lambda g: "F" if g == 1 else "M")
        #inf_sites["alt_origin"] = inf_sites["dad_gt"].apply(lambda g: "M" if g == 1 else "F")
        inf_sites = inf_sites.merge(het_sites_kid, on=["chrom", "pos"])
        if inf_sites.shape[0] == 0: continue
        #print (inf_sites)

        # now that we know the parental origins of the REF/ALT alleles, we can now determine
        # the parent-of-origin for each of the two haplotypes in the child. e.g., if the child's
        # phased genotype at this informative site is 0|1, and we know that mom was 1/1 and dad was 0/1
        # at this site, we know that mom donated the ALT allele and dad donated the REF allele. thus,
        # the A haplotype is derived from dad and the B haplotype is derived from mom.
        inf_sites["haplotype_A_parent"] = inf_sites.apply(lambda row: row["alt_parent"] if row["haplotype_A_allele"] == 1 else row["ref_parent"], axis=1)
        inf_sites["haplotype_B_parent"] = inf_sites.apply(lambda row: row["alt_parent"] if row["haplotype_B_allele"] == 1 else row["ref_parent"], axis=1)
        if inf_sites.shape[0] == 0: continue
        inf_sites["start"] = start
        inf_sites["end"] = end
        informative_sites.append(inf_sites)

informative_sites = pd.concat(informative_sites)

informative_sites.to_csv("inf_sites.csv", index=False)
