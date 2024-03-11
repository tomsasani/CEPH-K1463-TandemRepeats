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
import argparse

# def catalog_informative_sites(
#     vcf: VCF,
#     region: str,
#     parent_dict: Dict,
# ):
#     res = []
#     for v in vcf(region):
#         if not var_pass(v): continue
#         if f"{v.CHROM}:{v.POS}" not in parent_dict: continue

#         gts = v.genotypes
#         hap_0, hap_1, is_phased = gts[0]
#         if not is_phased: continue
#         if hap_0 + hap_1 != 1: continue

#         (
#             dad_hap_0,
#             dad_hap_1,
#             dad_phased,
#             mom_hap_0,
#             mom_hap_1,
#             mom_phased,
#         ) = parent_dict[f"{v.CHROM}:{v.POS}"].values()

#         dad_gt, mom_gt = dad_hap_0 + dad_hap_1, mom_hap_0 + mom_hap_1

#         if dad_gt == mom_gt: continue

#         # figure out origin of haplotype A and B
#         hap_0_origin, hap_1_origin = None, None
#         if dad_gt > mom_gt:
#             hap_0_origin = "dad" if hap_0 == 1 else "mom"
#             hap_1_origin = "dad" if hap_1 == 1 else "mom"
#         elif mom_gt > dad_gt:
#             hap_0_origin = "mom" if hap_0 == 1 else "dad"
#             hap_1_origin = "mom" if hap_1 == 1 else "dad"

#         dad_haplotype_origin, mom_haplotype_origin = "?", "?"
#         if dad_gt > mom_gt and dad_phased:
#             dad_haplotype_origin = "A" if dad_hap_0 == 1 else "B"
#         elif mom_gt > dad_gt and mom_phased:
#             mom_haplotype_origin = "A" if mom_hap_0 == 1 else "B"
#         elif dad_gt < mom_gt and dad_phased:
#             dad_haplotype_origin = "A" if dad_hap_0 == 0 else "B"
#         elif mom_gt < dad_gt and mom_phased:
#             mom_haplotype_origin = "A" if mom_hap_0 == 0 else "B"

#         out_dict = {
#             "chrom": v.CHROM,
#             "pos": v.POS,
#             "dad_gt": "|".join((str(dad_hap_0), str(dad_hap_1))) if dad_phased else "/".join((str(dad_hap_0), str(dad_hap_1))),
#             "mom_gt": "|".join((str(mom_hap_0), str(mom_hap_1))) if mom_phased else "/".join((str(mom_hap_0), str(mom_hap_1))),
#             "kid_gt": "|".join((str(hap_0), str(hap_1))),
#             "haplotype_A_origin": hap_0_origin,
#             "haplotype_B_origin": hap_1_origin,
#             "dad_haplotype_origin": dad_haplotype_origin,
#             "mom_haplotype_origin": mom_haplotype_origin,
#         }
#         #print (out_dict)
#         #break

#         res.append(out_dict)

#     return pd.DataFrame(res)

def var_pass(v: cyvcf2.Variant, idxs: np.ndarray):

    GT2ALT_AB = {0: (0., 0.05), 1: (0.2, 0.8), 2: (0.95, 1.)}

    if v.var_type != "snp": return False
    if v.call_rate < 1.: return False
    if np.any(v.gt_quals[idxs] < 20): return False
    if len(v.ALT) > 1: return False
    rd, ad = v.gt_ref_depths, v.gt_alt_depths
    td = ad + rd
    ab = ad / td
    if np.any(td[idxs] < 10): return False

    gts = v.gt_types
    for idx in idxs:
        min_ab, max_ab = GT2ALT_AB[gts[idx]]
        if ab[idx] < min_ab or ab[idx] > max_ab: 
            return False

    return True

def catalog_informative_sites(
    vcf: VCF,
    region: str,
    dad: str,
    mom: str,
    child: str,
    smp2idx: Dict[str, int],
):

    dad_idx, mom_idx = smp2idx[dad], smp2idx[mom]
    kid_idx = smp2idx[child]

    res = []
    for v in vcf(region):
        if not var_pass(
            v,
            np.array([dad_idx, mom_idx, kid_idx]),
        ):
            continue

        gts = v.genotypes

        dad_hap_0, dad_hap_1, dad_phased = gts[dad_idx]
        mom_hap_0, mom_hap_1, mom_phased = gts[mom_idx]
        kid_hap_0, kid_hap_1, kid_phased = gts[kid_idx]

        dad_gt, mom_gt = (
            dad_hap_0 + dad_hap_1,
            mom_hap_0 + mom_hap_1,
        )

        if dad_gt == mom_gt: continue

        if not kid_phased: continue
        if kid_hap_0 + kid_hap_1 != 1: continue

        if not dad_phased: 
            if dad_hap_0 != dad_hap_1: continue
        if not mom_phased: 
            if mom_hap_0 != mom_hap_1: continue

        hap_0_origin, hap_1_origin = None, None
        if dad_gt > mom_gt:
            hap_0_origin = "dad" if kid_hap_0 == 1 else "mom"
            hap_1_origin = "dad" if kid_hap_1 == 1 else "mom"
        elif mom_gt > dad_gt:
            hap_0_origin = "mom" if kid_hap_0 == 1 else "dad"
            hap_1_origin = "mom" if kid_hap_1 == 1 else "dad"

        dad_haplotype_origin, mom_haplotype_origin = "?", "?"
        if dad_gt > mom_gt and dad_phased:
            dad_haplotype_origin = "A" if dad_hap_0 == 1 else "B"
        elif mom_gt > dad_gt and mom_phased:
            mom_haplotype_origin = "A" if mom_hap_0 == 1 else "B"
        elif dad_gt < mom_gt and dad_phased:
            dad_haplotype_origin = "A" if dad_hap_0 == 0 else "B"
        elif mom_gt < dad_gt and mom_phased:
            mom_haplotype_origin = "A" if mom_hap_0 == 0 else "B"

        out_dict = {
            "chrom": v.CHROM,
            "pos": v.POS,
            "dad_gt": "|".join((str(dad_hap_0), str(dad_hap_1))) if dad_phased else "/".join((str(dad_hap_0), str(dad_hap_1))),
            "mom_gt": "|".join((str(mom_hap_0), str(mom_hap_1))) if mom_phased else "/".join((str(mom_hap_0), str(mom_hap_1))),
            "kid_gt": "|".join((str(kid_hap_0), str(kid_hap_1))),
            "haplotype_A_origin": hap_0_origin,
            "haplotype_B_origin": hap_1_origin,
            "dad_haplotype_origin": dad_haplotype_origin,
            "mom_haplotype_origin": mom_haplotype_origin,
        }
        res.append(out_dict)

    return pd.DataFrame(res)

def main(args):

    SNP_VCF = VCF(args.phased_snp_vcf, gts012=True)

    SMP2IDX = dict(zip(SNP_VCF.samples, range(len(SNP_VCF.samples))))

    CHROMS = [f"chr{c}" for c in range(1, 23)]

    informative_sites = []

    with open(args.phase_blocks, "r") as infh:
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

            
            informative_phases = catalog_informative_sites(
                SNP_VCF,
                region,
                args.dad_id,
                args.mom_id,
                args.focal,
                SMP2IDX,
            )

            if informative_phases.shape[0] == 0: continue
            informative_phases["phase_block_chrom"] = chrom
            informative_phases["phase_block_start"] = start
            informative_phases["phase_block_end"] = end
            informative_sites.append(informative_phases)

    informative_sites = pd.concat(informative_sites)

    informative_sites.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--phased_snp_vcf")
    p.add_argument("--phase_blocks")
    p.add_argument("--focal")
    p.add_argument("--out")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")
    args = p.parse_args()
    main(args)
