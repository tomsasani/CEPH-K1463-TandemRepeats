from cyvcf2 import VCF
import pandas as pd
import tqdm
import numpy as np 
import argparse
from typing import List, Dict
import cyvcf2
from collections import namedtuple


def var_pass(
    v: cyvcf2.Variant,
    idxs: np.ndarray,
    min_depth: int = 10,
    min_gq: int = 20,
) -> bool:
    """method to ensure that informative
    variants meet basic filtering criteria.

    Args:
        v (cyvcf2.Variant): cyvcf2.Variant object
        idxs (np.ndarray): indices of kid, mom, and dad in the VCF FORMAT field
        min_depth (int): minimum depth required in all members of trio. Defaults to 10.
        min_gq: minimum genotype quality required in all members of trio. Defaults to 20.

    Returns:
        bool: boolean indicating whether or not variant passes filters.
    """

    # map genotypes to expected allele balance range
    GT2ALT_AB = {0: (0., 0.), 1: (0.2, 0.8), 2: (1, 1.)}

    if v.var_type != "snp": return False
    # if the PF tag is set to TR_OVERLAP, then the SNV
    # overlaps a TR and shouldn't be used
    if v.format("PF") is not None:
        if "TR_OVERLAP" in v.format("PF"): return False
    
    # remove multi-allelic variants
    if len(v.ALT) > 1: return False
    # require minimum genotype quality and depth
    gq = v.gt_quals[idxs]
    if np.any(gq < min_gq): return False

    ad, rd = v.gt_alt_depths, v.gt_ref_depths
    td = ad + rd
    # calculate allele balance, and require reasonable
    # allele balance for HOM_REF/HET/HOM_ALT
    ab = ad / td
    gts = v.gt_types
    # filter UNK genotypes
    if np.any(gts[idxs] == 3): return False
    for idx in idxs:
        min_ab, max_ab = GT2ALT_AB[gts[idx]]
        if not (min_ab <= ab[idx] <= max_ab):
            return False

    return True


def find_informative_sites_in_parents(
    *,
    vcf: VCF,
    region: str,
    dad: str,
    mom: str,
    child: str,
    smp2idx: Dict[str, int],
    is_male_x: bool,
    is_male_y: bool,
) -> pd.DataFrame:
    """given a VCF file and a region of interest, catalog all sites that
    are 'informative' for determining the phase of a candidate variant.

    Args:
        vcf (VCF): cyvcf2.VCF object, should be a joint-genotyped VCF with all samples.
        region (str): region to query, formatted as chr:start-end.
        dad (str): sample ID of dad.
        mom (str): sample ID of mom.
        child (str): sample ID of child.
        smp2idx (Dict[str, int]): dictionary mapping sample IDs to integer indices in VCF FORMAT field
        is_male_x (bool): whether the region is on the X chromosome AND the child is male.
        is_male_y (bool): whether the region is on the Y chromosome AND the child is male.

    Returns:
        pd.DataFrame: pandas DataFrame containing relevant information about informative sites in the region.
    """

    dad_idx, mom_idx = smp2idx[dad], smp2idx[mom]
    kid_idx = smp2idx[child]

    res = []
    for v in vcf(region):
        if not var_pass(v, np.array([kid_idx, dad_idx, mom_idx])): 
            continue
        # access genotypes in kid, mom, and dad. don't need to be
        # phased, as long as they are informative.
        gts = v.gt_types
        dad_gt, mom_gt, kid_gt = (
            gts[dad_idx],
            gts[mom_idx],
            gts[kid_idx],
        )
        # make sure parents don't have the same genotype.
        if dad_gt == mom_gt: continue
        # make sure kid is HET (or HOM if a male sex chromosome)
        if is_male_x or is_male_y:
            if kid_gt != 2: continue
        else:
            if kid_gt != 1: continue

        # create namedtuple object with information about this site
        InfSite = namedtuple("InfSite", "chrom pos dad_gt mom_gt is_male_x is_male_y")
        inf_site = InfSite(v.CHROM, v.POS, dad_gt, mom_gt, is_male_x, is_male_y)

        res.append(inf_site)

    return res

def adjudicate_str_and_snv():
    """given a phased STR genotype and a phased SNV genotype that are present
    in the same HiPhase phase block, 
    """
    pass

def phase_informative_sites(
    *,
    inf_sites: List[namedtuple],
    vcf: cyvcf2.VCF,
    smp2idx: Dict[str, int],
    child: str,
) -> pd.DataFrame:
    """given a VCF file and a region of interest, catalog all sites that
    are 'informative' for determining the phase of a candidate variant.

    Args:
        vcf (VCF): cyvcf2.VCF object, representing a phased SNP VCF in the kid.
        

    Returns:
        pd.DataFrame: pandas DataFrame containing relevant information about informative sites in the region.
    """

    kid_idx = smp2idx[child]
    res = []

    # loop over the list of informative sites surrounding the TR DNM
    for site in inf_sites:
        
        chrom, start, end = site.chrom, site.pos - 1, site.pos
        is_male_x, is_male_y = site.is_male_x, site.is_male_y

        # we require the kid to be HET, so we don't need to access their genotype
        dad_gt, mom_gt = site.dad_gt, site.mom_gt

        region = f"{chrom}:{start}-{end}"
        for v in vcf(region):
            # access phase block PS tags. for a given sample, if two
            # variants in two different VCFs share a PS tag, their phases
            # can be inferred to be the same. e.g., a 0|1 SNP and 0|1 STR
            # both have the non-reference allele on haplotype B.
            kid_ps = -1
            if is_male_x or is_male_y:
                pass
            else:
                ps_tags = v.format("PS")
                # only need a PS tag if we're not on a male X
                if ps_tags is None: continue
                ps_tags = ps_tags[:, 0]
                assert ps_tags.shape[0] == len(smp2idx)
                # get kid's PS tag
                kid_ps = ps_tags[kid_idx]

            # get phased genotype information in the child.
            phased_genotypes = v.genotypes
            kid_hap_0, kid_hap_1, kid_phased = phased_genotypes[kid_idx]
            # only need kid to be phased if we're *not* on a male sex chromosome.
            if is_male_x or is_male_y:
                pass
            else: 
                if not kid_phased: continue

            # if the kid is het, whichever parent has more
            # alt alleles is the one that donated the ALT.
            # even if one parent is UNK, we can assume they donated
            # the REF allele if we're on an autosome.
            hap_0_origin, hap_1_origin = None, None
            # determining origin is easy if we're on a male X or Y
            if is_male_x:
                hap_0_origin, hap_1_origin = "mom", "mom"
            elif is_male_y:
                hap_0_origin, hap_1_origin = "dad", "dad"
            else:
                if dad_gt > mom_gt:
                    hap_0_origin = "dad" if kid_hap_0 == 1 else "mom"
                    hap_1_origin = "dad" if kid_hap_1 == 1 else "mom"
                elif mom_gt > dad_gt:
                    hap_0_origin = "mom" if kid_hap_0 == 1 else "dad"
                    hap_1_origin = "mom" if kid_hap_1 == 1 else "dad"
                else:
                    continue
            # store results in dictionary
            out_dict = {
                "inf_chrom": v.CHROM,
                "inf_pos": v.POS,
                "dad_inf_gt": str(dad_gt), 
                "mom_inf_gt": str(mom_gt),
                "kid_inf_gt": "|".join((str(kid_hap_0), str(kid_hap_1))),
                "kid_inf_ps": kid_ps,
                "haplotype_A_origin": hap_0_origin,
                "haplotype_B_origin": hap_1_origin,
            }
            res.append(out_dict)

    return pd.DataFrame(res)


def measure_consistency(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """given a dataframe that contains information about informative sites
    surrounding a candidate DNM, find the longest stretch of 'consistent'
    informative sites that all support the same haplotype assignment.

    Args:
        df (pd.DataFrame): pandas DataFrame object containing information about informative sites.
        columns (List[str]): list of columns to use in consistency checks

    Returns:
        pd.DataFrame: pandas DataFrame containing a subset of only the informative sites in the
            longest continuous stretch of consistent sites.
    """

    res = []
    for (trid, genotype), trid_df in df.groupby(["trid", "genotype"]):
        # sort the informative sites by absolute distance to the STR
        trid_df_sorted = trid_df.sort_values(
            "abs_diff_to_str",
            ascending=True,
        ).reset_index(drop=True)

        sorted_values = trid_df_sorted[columns].values

        # figure out how many sites are consistent closest to the STR. we can do this
        # simply by figuring out the first index where they are *inconsistent.*
        inconsistent_phases = np.where(sorted_values[1:] != sorted_values[:-1])[0]

        # if we have more than 1 informative site and none of them are inconsistent,
        # then all of them are consistent
        if sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] == 0:
            trid_df_sorted_consistent = trid_df_sorted.iloc[:]
        # if we have more than 1 informative site and some of them are inconsistent,
        # return the consistent ones up until the first inconsistency.
        elif sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] > 0:
            trid_df_sorted_consistent = trid_df_sorted.iloc[:inconsistent_phases[0]]
        # if we only have one informative site
        elif sorted_values.shape[0] == 1:# and sorted_values[0] in ("mom", "dad"):
            trid_df_sorted_consistent = trid_df_sorted.iloc[:]
        # if we don't have any informative sites
        else:
            continue

        res.append(trid_df_sorted_consistent)
    return pd.concat(res, ignore_index=True)


def main(args):

    KID_STR_VCF = VCF(args.str_vcf, gts012=True)
    JOINT_SNP_VCF = VCF(args.joint_snp_vcf, gts012=True)
    KID_SNP_VCF = VCF(args.kid_snp_vcf, gts012=True)

    KID_STR_SMP2IDX = dict(zip(KID_STR_VCF.samples, range(len(KID_STR_VCF.samples))))
    kid_str_idx = KID_STR_SMP2IDX[args.focal_alt]

    SLOP = 50_000

    # CHROMS = [f"chr{c}" for c in range(1, 23)]
    # CHROMS.extend(["chrX", "chrY"])

    mutations = pd.read_csv(args.mutations, sep="\t")

    dnm_phases = []

    informative_sites: List[pd.DataFrame] = []

    # loop over STR de novos in the dataframe
    for _, row in tqdm.tqdm(mutations.iterrows()):
        row_dict = row.to_dict()
        # extract chrom, start and end
        trid = row["trid"]
        trid_chrom, trid_start, trid_end = row["#chrom"], row["start"], row["end"]
        trid_start, trid_end = int(trid_start) - 1, int(trid_end)

        is_male_x = row_dict["suffix"].startswith("S") and trid_chrom == "chrX"
        is_male_y = row_dict["suffix"].startswith("S") and trid_chrom == "chrY"

        als = list(map(int, row_dict["child_AL"].split(",")))
        idx = row_dict["index"]
        denovo_al = als[idx]

        phase_start, phase_end = trid_start - SLOP, trid_end + SLOP
        if phase_start < 1: phase_start = 1

        # collate informative SNV sites within a defined region surrounding the
        # STR. this is just a list of every site where the parent genotypes are
        # different, and the child is heterozygous.
        informative_parent_sites = find_informative_sites_in_parents(
                vcf=JOINT_SNP_VCF,
                region=f"{trid_chrom}:{phase_start}-{phase_end}",
                dad=args.dad_id,
                mom=args.mom_id,
                child=args.focal_alt,
                smp2idx=dict(zip(JOINT_SNP_VCF.samples, range(len(JOINT_SNP_VCF.samples)))),
                is_male_x=is_male_x,
                is_male_y=is_male_y,
            )

        if len(informative_parent_sites) == 0: 
            continue

        # using those informative sites, we can now search the kid's *phased*
        # SNV VCF to determine which of the kid's haplotypes (A or B) carries the informative
        # allele. that is, if mom is 0/1 and dad is 1/1, and the kid is 0|1, we know that
        # haplotype A is from mom and B is from dad.
        informative_phases = phase_informative_sites(
            inf_sites=informative_parent_sites,
            vcf=KID_SNP_VCF,
            smp2idx=dict(zip(KID_SNP_VCF.samples, range(len(KID_SNP_VCF.samples)))),
            child=args.focal_alt,
        )

        # add information about the region in which we searched for informative
        # sites to the output dataframe.
        informative_phases["trid"] = trid

        informative_sites.append(informative_phases)

        # now, we look in the kid's phased STR VCF to figure out which
        # haplotype the STR DNM occurred on
        for var in KID_STR_VCF(f"{trid_chrom}:{trid_start}-{trid_end}"):
            # ensure the VCF TRID matches the TR TRID
            if var.INFO.get("TRID") != trid: 
                continue
            # get phased genotype if we're not on a male sex chromosome
            if is_male_x or is_male_y:
                pass
            else:
                hap_a, hap_b, is_phased = var.genotypes[0]
                if not is_phased: continue

            allele_lengths = [len(var.REF)] + [len(a) for a in var.ALT]
            assert denovo_al == allele_lengths[row["genotype"]]
            
            # figure out the haplotype on which the de novo allele occurred
            if hap_a == row["genotype"]:
                denovo_hap_id = "A"
            elif hap_b == row["genotype"]:
                denovo_hap_id = "B"
            else:
                continue
            # access PS tag for child at the STR locus
            kid_ps = -1
            if is_male_x or is_male_y:
                pass
            else:
                ps_tags = var.format("PS")[:, 0]
                kid_ps = ps_tags[kid_str_idx]

            row_dict = row.to_dict()
            row_dict.update(
                {
                    "denovo_hap_id": denovo_hap_id,
                    "str_chrom": trid_chrom,
                    "kid_str_ps": kid_ps,
                    "kid_str_gt": "|".join((str(hap_a), str(hap_b))),
                    "kid_str_hap_a": hap_a,
                    "kid_str_hap_b": hap_b,
                }
            )
            dnm_phases.append(row_dict)

    informative_sites = pd.concat(informative_sites)

    dnm_phases = pd.DataFrame(dnm_phases)

    # merge de novo STR information with informative site information.
    # it's critical to ensure that we only merge STRs with informative sites
    # that share the same phase block (PS) tag, so that we can use phased
    # genotypes at informative sites to infer phases at STRs.
    merged_dnms_inf = dnm_phases.merge(
        informative_sites,
        left_on=["trid", "str_chrom", "kid_str_ps"],
        right_on=["trid", "inf_chrom", "kid_inf_ps"],
    )

    # measure the consistency of the haplotype phasing. subset the dataframe to only
    # include the N informative sites that share a consistent haplotype origin
    # immediately surrounding the STR locus.
    # first, figure out the distance between each informative site and the de novo STR
    merged_dnms_inf["str_midpoint"] = merged_dnms_inf.apply(lambda row: np.mean([int(row["start"]), int(row["end"])]), axis=1)
    merged_dnms_inf["diff_to_str"] = merged_dnms_inf["inf_pos"] - merged_dnms_inf["str_midpoint"]
    merged_dnms_inf["abs_diff_to_str"] = np.abs(merged_dnms_inf["diff_to_str"])

    merged_dnms_inf["str_parent_of_origin"] = merged_dnms_inf.apply(
        lambda row: row["haplotype_{}_origin".format(row["denovo_hap_id"])],
        axis=1,
    )

    merged_dnms_inf_consistent = measure_consistency(merged_dnms_inf, "str_parent_of_origin")

    # now that we have a consistent set of informative sites surrounding the STR,
    # infer parent of origin for the haplotype that carries the de novo allele.
    # by figuring out the haplotype ID on which the de novo STR occurred (A or B)
    # and then asking which parent consistently donated the A or B haplotype at
    # the surrounding informative sites.
    

    merged_dnms_inf_consistent.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--joint_snp_vcf")
    p.add_argument("--kid_snp_vcf")
    p.add_argument("--focal")
    p.add_argument("--focal_alt")
    p.add_argument("--out")
    p.add_argument("--str_vcf")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")
    args = p.parse_args()
    main(args)
