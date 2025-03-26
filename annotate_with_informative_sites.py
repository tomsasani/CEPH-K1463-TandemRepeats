from cyvcf2 import VCF
import pandas as pd
import tqdm
import numpy as np 
import argparse
from typing import List, Dict
import cyvcf2
from collections import namedtuple

from schema import TRGTDeNovoSchema

def var_pass(
    v: cyvcf2.Variant,
    idxs: np.ndarray,
    min_depth: int = 10,
    min_gq: int = 20,
    is_male_x: bool = False,
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
    if np.any(td < min_depth): return False
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
        # make sure kid is HET
        if kid_gt != 1: continue
        # create namedtuple object with information about this site
        InfSite = namedtuple("InfSite", "chrom pos dad_gt mom_gt")
        inf_site = InfSite(v.CHROM, v.POS, dad_gt, mom_gt)

        res.append(inf_site)

    return res


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
        # we require the kid to be HET, so we don't need to access their genotype
        dad_gt, mom_gt = site.dad_gt, site.mom_gt
        region = f"{chrom}:{start}-{end}"
        for v in vcf(region):
        
            # access phase block PS tags. for a given sample, if two
            # variants in two different VCFs share a PS tag, their phases
            # can be inferred to be the same. e.g., a 0|1 SNP and 0|1 STR
            # both have the non-reference allele on haplotype B.
            ps_tags = v.format("PS")
            if ps_tags is None: continue
            ps_tags = ps_tags[:, 0]
            assert ps_tags.shape[0] == len(smp2idx)
            # get kid's PS tag
            kid_ps = ps_tags[kid_idx]
            # get phased genotype information in the child. if the kid isn't
            # phased in the SNV VCF, move on.
            phased_genotypes = v.genotypes
            kid_hap_0, kid_hap_1, kid_phased = phased_genotypes[kid_idx]
            if not kid_phased: continue
            # if the kid is het, whichever parent has more
            # alt alleles is the one that donated the ALT.
            # even if one parent is UNK, we can assume they donated
            # the REF allele if we're on an autosome.
            hap_0_origin, hap_1_origin = None, None            
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

    # sort the informative sites by absolute distance to the STR
    df_sorted = df.sort_values(
        "abs_diff_to_str",
        ascending=True,
    ).reset_index(drop=True)

    sorted_values = df_sorted[columns].values

    # figure out how many sites are consistent closest to the STR. we can do this
    # simply by figuring out the first index where they are *inconsistent.*
    inconsistent_phases = np.where(sorted_values[1:] != sorted_values[:-1])[0]

    # if we have more than 1 informative site and none of them are inconsistent,
    # then all of them are consistent
    if sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] == 0:
        df_sorted_consistent = df_sorted.iloc[:]
    # if we have more than 1 informative site and some of them are inconsistent,
    # return the consistent ones up until the first inconsistency.
    elif sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] > 0:
        df_sorted_consistent = df_sorted.iloc[:inconsistent_phases[0]]
    # if we only have one informative site
    else:
        df_sorted_consistent = df_sorted.iloc[:]

    return df_sorted_consistent


def main(args):

    # read in child's TRGT STR VCF (must be phased with HiPhase)
    KID_STR_VCF = VCF(args.str_vcf, gts012=True)
    # read in an SNV VCF with genotypes in (at least) the child's trio.
    # this VCF doesn't need to be phased.
    JOINT_SNP_VCF = VCF(args.joint_snp_vcf, gts012=True)
    # read in a VCF with the child's phased SNV genotypes
    KID_SNP_VCF = VCF(args.kid_snp_vcf, gts012=True)

    KID_STR_SMP2IDX = dict(zip(KID_STR_VCF.samples, range(len(KID_STR_VCF.samples))))
    kid_str_idx = KID_STR_SMP2IDX[args.focal_alt]

    # read in mutations and validate schema
    mutations = pd.read_csv(args.mutations, sep="\t", dtype={"sample_id": str})
    TRGTDeNovoSchema.validate(mutations)

    res: List = []

    # loop over STR de novos in the dataframe
    for _, row in tqdm.tqdm(mutations.iterrows()):
        row_dict = row.to_dict()

        # extract chrom, start and end
        trid = row["trid"]
        trid_chrom, trid_start, trid_end = row["#chrom"], row["start"], row["end"]
        trid_start, trid_end = int(trid_start) - 1, int(trid_end)

        # if this de novo is on a male sex chromosome, we don't even need to
        # look at informative sites -- by default, we know which parent's
        # germline experienced the de novo mutation. so we'll skip it for now.
        is_male_x = row_dict["suffix"].startswith("S") and trid_chrom == "chrX"
        is_male_y = row_dict["suffix"].startswith("S") and trid_chrom == "chrY"
        if is_male_x or is_male_y:
            continue

        # figure out the allele length of the de novo allele in the kid
        child_als = list(map(int, row_dict["child_AL"].split(",")))
        denovo_idx = row_dict["index"]
        denovo_al = child_als[denovo_idx]

        # define the bounds of the region in which we'll search for informative sites
        phase_start, phase_end = trid_start - args.slop, trid_end + args.slop
        if phase_start < 1: phase_start = 1

        # collate informative SNV sites within a defined region surrounding the
        # STR. this is just a list of every site where the parent genotypes are
        # different and the child is heterozygous.
        informative_sites = find_informative_sites_in_parents(
            vcf=JOINT_SNP_VCF,
            region=f"{trid_chrom}:{phase_start}-{phase_end}",
            dad=args.dad_id,
            mom=args.mom_id,
            child=args.focal_alt,
            smp2idx=dict(
                zip(JOINT_SNP_VCF.samples, range(len(JOINT_SNP_VCF.samples))),
            ),
        )
        # if we didn't find any informative sites *at all*, we can't determine
        # a reliable phase, and we can move on.
        if len(informative_sites) == 0: 
            continue

        # using those informative sites, we can now search the kid's *phased*
        # SNV VCF to determine which of the kid's haplotypes (A or B) carries the informative
        # allele. that is, if mom is 0/1 and dad is 1/1, and the kid is 0|1, we know that
        # haplotype A is from mom and B is from dad.
        informative_sites_phased = phase_informative_sites(
            inf_sites=informative_sites,
            vcf=KID_SNP_VCF,
            smp2idx=dict(
                zip(KID_SNP_VCF.samples, range(len(KID_SNP_VCF.samples))),
            ),
            child=args.focal_alt,
        )
        if len(informative_sites_phased) == 0:
            continue

        # now, we look in the kid's phased STR VCF to figure out which
        # haplotype the STR DNM occurred on
        for var in KID_STR_VCF(f"{trid_chrom}:{trid_start}-{trid_end}"):
            # ensure the VCF TRID matches the TR TRID
            assert var.INFO.get("TRID") == trid
            # get phased genotype for the STR
            hap_a, hap_b, is_phased = var.genotypes[0]
            # if we can't get a phased genotype for the STR, we can't
            # reliably phase it
            if not is_phased:
                continue

            allele_lengths = [len(var.REF)] + [len(a) for a in var.ALT]
            assert denovo_al == allele_lengths[row["genotype"]]
            # figure out the haplotype on which the de novo allele occurred
            denovo_hap_id = None
            if hap_a == row["genotype"]:
                denovo_hap_id = "A"
            elif hap_b == row["genotype"]:
                denovo_hap_id = "B"
            if denovo_hap_id is None:
                raise ValueError

            # access PS tag for child at the STR locus
            ps_tags = var.format("PS")[:, 0]
            kid_str_ps = ps_tags[kid_str_idx]

            # access the informative sites at which the PS tag in the
            # child's SNV VCF matches that of their phased STR genotype
            matching_inf_sites = informative_sites_phased[informative_sites_phased["kid_inf_ps"] == kid_str_ps]

            # sort those matching sites by their absolute distance to the STR
            str_midpoint = (trid_end + trid_start) // 2
            matching_inf_sites["diff_to_str"] = matching_inf_sites["inf_pos"] - str_midpoint
            matching_inf_sites["abs_diff_to_str"] = matching_inf_sites["diff_to_str"].apply(lambda d: abs(d))

            # then, infer the likely parent-of-origin for the STR at each site
            matching_inf_sites["str_parent_of_origin"] = matching_inf_sites.apply(
                lambda row: row["haplotype_{}_origin".format(denovo_hap_id)],
                axis=1,
            )
            
            # for each informative site, add informatiion about the TR locus and append
            # to output file
            for _, inf_row in matching_inf_sites.iterrows():
                inf_row_dict = inf_row.to_dict()
                inf_row_dict.update(row_dict)
                res.append(inf_row_dict)

    # combine informative sites with original mutation dataframe
    combined_informative_sites = pd.DataFrame(res).drop(
        columns=[
            "kid_inf_ps",
            "haplotype_A_origin",
            "haplotype_B_origin",
            "diff_to_str",
        ]
    )
    merged = mutations.merge(combined_informative_sites, how="left")

    # fill NA values for DNMs without any informative sites
    merged = merged.fillna(
        value={
            "inf_chrom": "UNK",
            "inf_pos": -1,
            "dad_inf_gt": -1,
            "mom_inf_gt": -1,
            "str_parent_of_origin": "UNK",
        }
    )

    merged.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        required=True,
        help="""TSV file containing TRGT-denovo output for the child with candidate DNMs""",
    )
    p.add_argument(
        "--joint_snp_vcf",
        required=True,
        help="""VCF file containing SNV genotypes for (at least) the child's corresponding trio. does not need to be phased with HiPhase.""",
    )
    p.add_argument(
        "--kid_snp_vcf",
        required=True,
        help="""VCF file containing phased SNV genotypes in the child.""",
    )
    p.add_argument(
        "--focal",
        required=True,
        help="""sample ID for the child (CEPH LABID)""",
    )
    p.add_argument(
        "--focal_alt", required=True, help="""sample ID for the child (Coriell ID)"""
    )
    p.add_argument(
        "--out",
        required=True,
        help="""name of output TSV file.""",
    )
    p.add_argument(
        "--str_vcf",
        required=True,
        help="""VCF file containing STR genotypes in the child. must be phased with HiPhase.""",
    )
    p.add_argument(
        "--dad_id",
        required=True,
        help="""Coriell ID of father in trio.""",
    )
    p.add_argument(
        "--mom_id",
        required=True,
        help="""Coriell ID of mother in trio""",
    )
    p.add_argument(
        "-slop",
        default=250_000,
        required=False,
        type=int,
        help="""how many base pairs upstream and downstream of the candidate STR DNM to look for informative sites for parent-of-origin inference.""",
    )
    args = p.parse_args()
    main(args)
