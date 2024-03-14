from cyvcf2 import VCF
import pandas as pd
import csv
import tqdm
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict, Counter
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 
import argparse
from utils import filter_mutation_dataframe
from collect_informative_sites import catalog_informative_sites


def measure_consistency(df: pd.DataFrame, column: str = "str_parent_of_origin"):

    res = []
    for trid, trid_df in df.groupby("trid"):
        # sort the informative sites by absolute distance to the STR
        trid_df_sorted = trid_df.sort_values("abs_diff_to_str", ascending=True).reset_index()

        # access the inferred phases at every informative site
        sorted_values = trid_df[column].values

        # figure out how many sites are consistent closest to the STR. we can do this
        # simply by figuring out the first index where they are *inconsistent.*
        inconsistent_phases = np.where(sorted_values[1:] != sorted_values[:-1])[0]

        # if none are inconsistent, all are consistent!
        if inconsistent_phases.shape[0] == 0:
            final_consistent = -1
        else:
            final_consistent = inconsistent_phases[0]

        trid_df_sorted_consistent = trid_df_sorted.iloc[:final_consistent]
        res.append(trid_df_sorted_consistent)

    return pd.concat(res)

def main(args):

    KID_STR_VCF = VCF(args.str_vcf, gts012=True)
    SNP_VCF = VCF(args.joint_snp_vcf, gts012=True)

    SLOP = 500_000

    SMP2IDX = dict(zip(SNP_VCF.samples, range(len(SNP_VCF.samples))))

    CHROMS = [f"chr{c}" for c in range(1, 23)]

    mutations = pd.read_csv(args.mutations, sep="\t")
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=True,
        remove_duplicates=True,
    )
    print (f"Total of {mutations.shape[0]} DNMs")

    dnm_phases = []

    WRONG_TRID, NOT_PHASED, NON_AUTOSOMAL = 0, 0, 0

    informative_sites = []

    # loop over mutations in the dataframe
    for _, row in tqdm.tqdm(mutations.iterrows()):
        row_dict = row.to_dict()
        # extract chrom, start and end
        trid = row["trid"]
        try:
            trid_chrom, trid_start, trid_end, _ = trid.split("_")
            trid_start, trid_end = int(trid_start) - 1, int(trid_end)
        except ValueError:
            dnm_phases.append(row_dict)
            continue


        if trid_chrom not in CHROMS:
            NON_AUTOSOMAL += 1
            dnm_phases.append(row_dict)
            continue

        phase_chrom = trid_chrom
        phase_start, phase_end = trid_start - SLOP, trid_end + SLOP
        if phase_start < 1: phase_start = 1

        # collate informative sites within a defined region surrounding the
        # STR. we record the PS (phase block) tag in the mom, dad, and child
        # at every informative site.
        informative_phases = catalog_informative_sites(
                SNP_VCF,
                f"{phase_chrom}:{phase_start}-{phase_end}",
                args.dad_id,
                args.mom_id,
                args.focal_alt,
                SMP2IDX,
            )
        if informative_phases.shape[0] == 0: 
            dnm_phases.append(row_dict)
            continue

        # add information about the region in which we searched for informative
        # sites to the output dataframe.
        informative_phases["phase_chrom"] = phase_chrom
        informative_phases["phase_start"] = phase_start
        informative_phases["phase_end"] = phase_end
        informative_sites.append(informative_phases)

        denovo_gt = row["genotype"]

        # loop over the STR VCF to get the focal STR entry
        for var in KID_STR_VCF(f"{trid_chrom}:{trid_start}-{trid_end}"):
            # ensure the VCF TRID matches the TR TRID
            if var.INFO.get("TRID") != trid: 
                WRONG_TRID += 1
                continue

            # get phased genotype
            is_phased = False
            try:
                hap_a, hap_b, is_phased = var.genotypes[0]
            except ValueError: 
                print ("### WARNING ###", var.genotypes[0])
                continue
            if not is_phased: 
                NOT_PHASED += 1
                dnm_phases.append(row_dict)
                continue

            # concatenate the alleles at this locus
            alleles = [var.REF] + var.ALT
            allele_lengths = [len(al) for al in alleles]

            # figure out haplotype ID on which DNM resides
            denovo_al, orig_al = None, None
            if hap_a == denovo_gt:
                denovo_hap_id = "A"
                denovo_al, orig_al = allele_lengths[hap_a], allele_lengths[hap_b]
            elif hap_b == denovo_gt:
                denovo_hap_id = "B"
                denovo_al, orig_al = allele_lengths[hap_b], allele_lengths[hap_a]
            else:
                continue

            # access PS tag for child at the STR locus
            kid_ps = var.format("PS")[:, 0][0]

            row_dict = row.to_dict()
            row_dict.update(
                {
                    "denovo_hap_id": denovo_hap_id,
                    "phase_chrom": phase_chrom,
                    "phase_start": phase_start,
                    "phase_end": phase_end,
                    "denovo_al": denovo_al,
                    "non_denovo_al": orig_al,
                    "kid_str_ps": kid_ps,
                    "kid_str_gt": "|".join((str(hap_a), str(hap_b))),
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
        left_on=["phase_chrom", "phase_start", "phase_end", "kid_str_ps"],
        right_on=["phase_chrom", "phase_start", "phase_end", "kid_inf_ps"],
    )

    # infer the parent of origin at every informative site. we do this
    # by figuring out the haplotype ID on which the de novo STR occurred (A or B)
    # and then asking which parent donated the A or B haplotype at the informative site.
    merged_dnms_inf["str_parent_of_origin"] = merged_dnms_inf.apply(
        lambda row: row["haplotype_{}_origin".format(row["denovo_hap_id"])],
        axis=1,
    )

    # for a given informative site, we can now ask if we were able
    # to determine the actual *haplotype* the parent-of-origin donated
    # to the child.

    # NOTE: it's critical to bear in mind that these haplotype origins are
    # only informative if they occur at sites with PS tags that match the PS
    # tag of the STR locus in the relevant parent. for example, at an informative site,
    # if the dad's genotype is 0|1 and mom's genotype is 0/0 and kid's genotype is 1|0,
    # we know that haplotype A in the child is paternal in origin. at the de novo locus,
    # if the dad's allele lengths are 12,15 and his genotype is 1|0, we know that the 15 AL
    # is on the same haplotype that he donated to the kid, ONLY IF the PS tag is the same
    # at the STR locus.
    # merged_dnms_inf["parent_haplotype_of_origin"] = merged_dnms_inf.apply(
    #     lambda row: row[f"{row['str_parent_of_origin']}_haplotype_origin"],
    #     axis=1,
    # )

    # merged_dnms_inf = merged_dnms_inf[merged_dnms_inf["parent_haplotype_of_origin"] != "?"]

    merged_dnms_inf["str_midpoint"] = merged_dnms_inf["trid"].apply(lambda t: np.mean(list(map(int, t.split("_")[1:-1]))))
    merged_dnms_inf["diff_to_str"] = merged_dnms_inf["pos"] - merged_dnms_inf["str_midpoint"]
    # figure out the distance between each informative site and the de novo STR
    merged_dnms_inf["abs_diff_to_str"] = np.abs(merged_dnms_inf["diff_to_str"])
    merged_dnms_inf["sample_id"] = args.focal

    # get a subset of the informative sites, including only those that were consistent
    # closest to the STR locus
    merged_dnms_inf_consistent = measure_consistency(merged_dnms_inf, column="str_parent_of_origin")

    merged_dnms_inf_consistent.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--joint_snp_vcf")
    p.add_argument("--focal")
    p.add_argument("--focal_alt")
    p.add_argument("--out")
    p.add_argument("--str_vcf")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")
    args = p.parse_args()
    main(args)
