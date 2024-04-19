from cyvcf2 import VCF
import cyvcf2
import pandas as pd
import tqdm
from typing import List, Dict
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import argparse
from utils import filter_mutation_dataframe

def trid_to_region(trid: str):
    chrom, start, end, _ = trid.split("_")
    return f"{chrom}:{start}-{end}"

def variant_pass(v: cyvcf2.Variant, idxs: np.ndarray, min_gq: int = 20, min_dp: int = 10,):
    if v.var_type != "snp": 
        return False
    if v.call_rate < 1.: return False
    if np.any(v.gt_quals[idxs] < min_gq): 
        return False
    rd, ad = v.gt_ref_depths, v.gt_alt_depths
    td = rd + ad
    if np.any(td[idxs] < min_dp): 
        return False
    return True

def ab_ok(idx: int, rd: np.ndarray, ad: np.ndarray):
    ab = ad[idx] / (rd[idx] + ad[idx])
    return 0.3 <= ab <= 0.7


def catalog_informative_sites(
    *,
    vcf: VCF,
    region: str,
    mom: str,
    dad: str,
    focal: str,
    focal_spouse: str,
    kids: List[str],
    kids_with_str: List[str],
    smp2idx: Dict[str, int],
    min_gq: int = 20,
    min_dp: int = 10,
):
    all_smps = [mom, dad, focal, focal_spouse] + kids
    all_idxs = np.array([smp2idx[s] for s in all_smps])

    mom_idx, dad_idx = smp2idx[mom], smp2idx[dad]
    focal_idx, spouse_idx = smp2idx[focal], smp2idx[focal_spouse]

    res = []
    for v in vcf(region):
        # only look at SNPs that pass
        if not variant_pass(v, all_idxs, min_gq=min_gq, min_dp=min_dp):
            continue
        
        gts = v.gt_types
        rd, ad = v.gt_ref_depths, v.gt_alt_depths

        # figure out informative parent
        inf_parent = None
        if gts[mom_idx] == 0 and gts[dad_idx] == 1 and ab_ok(dad_idx, rd, ad):
            inf_parent = "dad"
        elif gts[mom_idx] == 1 and gts[dad_idx] == 0 and ab_ok(mom_idx, rd, ad):
            inf_parent = "mom"
        if inf_parent is None: continue

        # make sure kid is HET
        if not (gts[focal_idx] == 1 and ab_ok(focal_idx, rd, ad)): continue
        # make sure spouse is HOM_REF
        if gts[spouse_idx] != 0: continue

        # loop over kids, figure out who has informative allele
        inf_string = []
        for k in kids:
            ki = smp2idx[k]
            # if the kid that inherited the informative site also has the STR, 
            # make a note. if more kids inherit the informative allele than inherit the DNM,
            # we may be looking at a post-zygotic.
            has_str = "Y" if k in kids_with_str else "N"
            if gts[ki] == 1 and ab_ok(ki, rd, ad):
                inf_string.append(f"{k}-{has_str}")
        if len(inf_string) == 0: continue
        res.append(f"{v.CHROM}:{v.POS}:{inf_parent}:{'|'.join(inf_string)}")

    return res


def check_for_dnm_inheritance(inf: str):
    
    children = inf.split(":")[-1]
    has_dnm = [k.split("-")[-1] for k in children.split("|")]
    # only interested in the inheritance patterns for which at least one kid inherited
    # the DNM allele
    return all([h == "N" for h in has_dnm])


def main(args):
    # assume we have a VCF with SNP genotypes for all individuals in the pedigree
    SNP_VCF = VCF(args.joint_snp_vcf, gts012=True)
    SMP2IDX = dict(zip(SNP_VCF.samples, range(len(SNP_VCF.samples))))

    PED_FILE = args.ped
    ped = pd.read_csv(
        PED_FILE,
        sep=",",
        dtype={"paternal_id": str, "maternal_id": str, "sample_id": str},
    )

    SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

    CHILDREN = ped[(ped["paternal_id"] == args.focal) | (ped["maternal_id"] == args.focal)]["sample_id"].to_list()
    EX_CHILD = ped[ped["sample_id"] == CHILDREN[0]]
    SPOUSE = EX_CHILD["paternal_id"].values[0] if EX_CHILD["paternal_id"].values[0] != args.focal else EX_CHILD["maternal_id"].values[0]

    SPOUSE = SMP2ALT[SPOUSE]
    CHILDREN = [SMP2ALT[c] for c in CHILDREN]

    # define the amount of slop (in bp) around a de novo to look for informative sites
    SLOP = 250_000

    INH_COL = "children_with_denovo_allele"

    transmission = pd.read_csv(
        args.transmission,
        sep=";",
        dtype={
            "sample_id": str,
            INH_COL: str,
        },
    ).dropna(subset=[INH_COL])

    # filter to the samples of interest
    transmission = transmission[transmission["sample_id"] == args.focal]

    mutations = pd.read_csv(args.mutations, sep="\t", dtype={"sample_id": str})

    mutations = mutations.merge(transmission)#, on=["trid", "index", "denovo_coverage"])

    # remove pathogenics
    mutations["region"] = mutations["trid"].apply(lambda t: trid_to_region(t))
    mutations["chrom"] = mutations["region"].apply(lambda r: r.split(":")[0])
    mutations = mutations[mutations["region"] != "UNK"]

    mutations["n_with_denovo_allele"] = mutations[INH_COL].apply(lambda c: len(c.split(",")))
    mutations = mutations[mutations["n_with_denovo_allele"] >= 1]

    # loop over mutations in the dataframe
    res = []
    for chrom, chrom_df in mutations.groupby("chrom"):

        for _, row in tqdm.tqdm(chrom_df.iterrows(), total=chrom_df.shape[0]):
            # extract informative sites from the parents
            region = row["region"]
            chrom = region.split(":")[0]
            start, end = list(map(int, region.split(":")[1].split("-")))
            adj_start, adj_end = start - SLOP, end + SLOP
            if adj_start < 0:
                adj_start = 1

            slop_region = f"{chrom}:{adj_start}-{adj_end}"

            phase_info = catalog_informative_sites(
                vcf=SNP_VCF,
                region=slop_region,
                mom=args.mom_id,
                dad=args.dad_id,
                focal=SMP2ALT[args.focal],
                focal_spouse=SPOUSE,
                kids=CHILDREN,
                kids_with_str=[SMP2ALT[k] for k in row[INH_COL].split(",")],
                smp2idx=SMP2IDX,
                min_gq=20,
                min_dp=10,
            )
            row_dict = row.to_dict()
            if len(phase_info) == 0:
                most_common_hap, most_common_freq, is_pz = "UNK", 0., False
            else:
                # count up most common inheritance pattern, ignoring cases where *none* of the
                # children inherited the DNM
                relevant_phase_infos = []
                for p in phase_info:
                    if check_for_dnm_inheritance(p): continue
                    relevant_phase_infos.append(p)
                if len(relevant_phase_infos) == 0:
                    most_common_hap, most_common_freq, is_pz = "UNK", 0., False
                else:
                    inheritance_combos = Counter([":".join(p.split(":")[2:]) for p in relevant_phase_infos]).most_common()
                    total_combos = len(relevant_phase_infos)

                    most_common_hap = inheritance_combos[0][0]
                    most_common_freq = inheritance_combos[0][1] / total_combos
                    most_common_children = [c.split("-")[-1] for c in most_common_hap.split(":")[-1].split("|")]
                    is_pz = len(set(most_common_children)) == 2

            row_dict["most_common_hap"] = most_common_hap
            row_dict["most_common_freq"] = most_common_freq
            row_dict["candidate_postzygotic"] = is_pz
            row_dict["n_inf"] = len(phase_info)

            res.append(row_dict)

    res_df = pd.DataFrame(res)

    res_df["phase"] = res_df["most_common_hap"].apply(lambda c: c.split(":")[0] if c != "UNK" and len(c.split(":")[1]) > 0 else "unknown")
    res_df["consistent"] = res_df["most_common_freq"] >= 0.8
    #res_df["true_phase"] = res_df.apply(lambda row: row["true_phase"] if row["consistent"] else "UNK", axis=1)
    # res_df.drop(
    #     columns=[
    #         INH_COL,
    #         "region",
    #         "n_with_denovo_allele",
    #     ]
    # )res_df.to_csv(args.out, sep="\t", index=False)
    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("-transmission")
    p.add_argument("--ped")
    p.add_argument("--focal")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")

    p.add_argument("--out")
    p.add_argument("--joint_snp_vcf")
    p.add_argument("-mutation_type")

    args = p.parse_args()
    main(args)
