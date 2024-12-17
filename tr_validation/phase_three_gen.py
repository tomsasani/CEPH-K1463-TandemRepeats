from cyvcf2 import VCF
import cyvcf2
import pandas as pd
import tqdm
from typing import List, Dict
from collections import Counter
import numpy as np
import argparse

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

def ab_ok(idx: int, rd: np.ndarray, ad: np.ndarray, gt: int):
    ab = ad[idx] / (rd[idx] + ad[idx])
    if gt == 1:
        return 0.3 <= ab <= 0.7
    elif gt == 2:
        return ab >= 0.9
    else: return False

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

        # figure out informative parent. as long as parental genotypes
        # don't match AND the kid is HET, we have an informative site.
        # we know for a fact that if the kid is HET, the parent with more ALTs
        # donated the ALT allele (e.g., if kid is 0/1, dad is 1/1, and mom is 0/1,
        # dad donated the 1 and mom donated the 0).
        inf_parent = None
        dad_gt, mom_gt = gts[dad_idx], gts[mom_idx]
        if dad_gt > mom_gt and ab_ok(dad_idx, rd, ad, dad_gt):
            inf_parent = "dad"
        elif mom_gt > dad_gt and ab_ok(mom_idx, rd, ad, mom_gt):
            inf_parent = "mom"
        if inf_parent is None: continue

        # make sure kid is HET
        if not (gts[focal_idx] == 1 and ab_ok(focal_idx, rd, ad, gts[focal_idx])): continue
        
        # make sure this site is informative w/r/t spouse, as well.
        # kid and spouse must have different genotypes -- in this code, 
        # we require that kid is HET and spouse in HOM_REF.
        if gts[spouse_idx] != 0: continue

        # loop over kids, figure out who has informative allele
        inf_string = []
        for k in kids:
            ki = smp2idx[k]
            # if the kid that inherited the informative site also has the STR, 
            # make a note. if some kids inherit the informative allele but don't inherit the DNM,
            # we may be looking at a post-zygotic.
            has_str = "Y" if k in kids_with_str else "N"
            if gts[ki] == 1 and ab_ok(ki, rd, ad, gts[ki]):
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

    INH_COL = "children_with_denovo_allele_strict"

    transmission = pd.read_csv(
        args.transmission,
        sep=";",
        dtype={
            "sample_id": str,
            INH_COL: str,
        },
    ).dropna(subset=[INH_COL, "children_consistent"])

    # filter to the samples of interest
    transmission = transmission[transmission["sample_id"] == args.focal]

    # we'll only look at sites at which *no G4 samples had an additional DNM at this locus*. in other
    # words, every G4 sample has to be mendelian consistent with the G3 parents.
    # transmission = transmission[(transmission["n_children"] == transmission["n_children_consistent"])]

    mutations = pd.read_csv(args.mutations, sep="\t", dtype={"sample_id": str})

    mutations = mutations.merge(transmission)

    mutations["region"] = mutations["trid"].apply(lambda t: trid_to_region(t))
    mutations = mutations[mutations["region"] != "UNK"]

    mutations["n_with_denovo_allele"] = mutations[INH_COL].apply(lambda c: len(c.split(",")))
    mutations = mutations[mutations["n_with_denovo_allele"] >= 1]

    # loop over mutations in the dataframe
    res = []
    for chrom, chrom_df in mutations.groupby("#chrom"):

        for _, row in tqdm.tqdm(chrom_df.iterrows(), total=chrom_df.shape[0]):
            # extract informative sites from the parents
            region = row["region"]
            chrom = region.split(":")[0]
            start, end = list(map(int, region.split(":")[1].split("-")))
            adj_start, adj_end = start - SLOP, end + SLOP
            if adj_start < 0:
                adj_start = 1

            slop_region = f"{chrom}:{adj_start}-{adj_end}"

            # get the list of kids that inherited the STR from the focal individual.
            kids_with_str = [SMP2ALT[k] for k in row[INH_COL].split(",")]
            # get the list of kids that are mendelian consistent with the focal individual.
            kids_mendel_consist = [SMP2ALT[k] for k in row["children_consistent"].split(",")]
            # subset to kids that both inherited the STR AND are mendelian consistent (i.e., 
            # don't have evidence of an additional DNM at this locus)
            kids_with_str_consist = list(set(kids_with_str).intersection(set(kids_mendel_consist)))
        
            phase_info = catalog_informative_sites(
                vcf=SNP_VCF,
                region=slop_region,
                mom=args.mom_id,
                dad=args.dad_id,
                focal=SMP2ALT[args.focal],
                focal_spouse=SPOUSE,
                kids=CHILDREN,
                kids_with_str=kids_with_str_consist,
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
                    # keep track of the "inheritance combos" -- that is, the number of instances where
                    # a particular informative allele (known to come from dad or mom) was inherited by 
                    # a list of children in the family.
                    inheritance_combos = Counter([":".join(p.split(":")[2:]) for p in relevant_phase_infos]).most_common()
                    total_combos = len(relevant_phase_infos)

                    most_common_hap = inheritance_combos[0][0]
                    most_common_freq = inheritance_combos[0][1] / total_combos
                    # most common hap is formatted like
                    # dad:NA12886-Y|NA12887-Y|NA12882-N|NA12879-Y|NA12883-Y
                    # keep track of whether all of the children that inherited the informative
                    # allele at this site also inherited the DNM
                    children_inherited_both = [c.split("-")[-1] for c in most_common_hap.split(":")[-1].split("|")]
                    is_pz = len(set(children_inherited_both)) == 2

            row_dict["most_common_hap"] = most_common_hap
            row_dict["most_common_freq"] = most_common_freq
            row_dict["candidate_postzygotic"] = is_pz
            row_dict["n_inf"] = len(phase_info)

            res.append(row_dict)

    res_df = pd.DataFrame(res)

    res_df["phase"] = res_df["most_common_hap"].apply(lambda c: c.split(":")[0] if c != "UNK" and len(c.split(":")[1]) > 0 else "unknown")
    
    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--transmission")
    p.add_argument("--ped")
    p.add_argument("--focal")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")

    p.add_argument("--out")
    p.add_argument("--joint_snp_vcf")
    p.add_argument("-mutation_type")

    args = p.parse_args()
    main(args)
