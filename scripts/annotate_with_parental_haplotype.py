import pandas as pd
from cyvcf2 import VCF
import numpy as np
import tqdm
import argparse
from collections import namedtuple

from schema import InformativeSiteSchema
import utils

def main(args):

    # read in dataframe with informative sites for every candidate DNM
    dtypes = utils.DTYPES.copy()
    dtypes.update({"sample_id": "string"})

    informative_sites = pd.read_csv(
        args.annotated_dnms, sep="\t", dtype=dtypes,
    )
    InformativeSiteSchema.validate(informative_sites)

    # file handles for (p/m)aternal SNV and STR VCFs (phased)
    dad_str_vcf = VCF(args.dad_phased_str_vcf)
    dad_snv_vcf = VCF(args.dad_phased_snv_vcf)
    mom_str_vcf = VCF(args.mom_phased_str_vcf)
    mom_snv_vcf = VCF(args.mom_phased_snv_vcf)

    res = []

    # loop over every informative site
    for i, row in tqdm.tqdm(informative_sites.iterrows()):
        # figure out the inferred parent of origin for this STR
        parent_of_origin = row["str_parent_of_origin"]
        if parent_of_origin not in ("mom", "dad"):
            continue

        # figure out the genotype in the informative and uninformative parent
        informative_parent_gt, other_parent_gt = None, None
        if parent_of_origin == "dad":
            informative_parent_gt = row["dad_inf_gt"]
            other_parent_gt = row["mom_inf_gt"]
        else:
            informative_parent_gt = row["mom_inf_gt"]
            other_parent_gt = row["dad_inf_gt"]

        # if the informative parent's genotype isn't HET, we can't infer
        # a haplotype of origin, since both haplotypes are identical and
        # we won't have phasing information
        if informative_parent_gt != 1:
            continue


        # get the location of this informative SNP
        chrom, start, end = row["inf_chrom"], row["inf_pos"] - 1, row["inf_pos"]
        snv_region = f"{chrom}:{start}-{end}"
        # get the location of the TR locus
        str_region = f"{row['#chrom']}:{row['start']}-{row['end']}"

        # define the parental SNV and STR VCFs that we'll need to query given
        # the parent that's informative here
        snv_vcf2query = dad_snv_vcf if parent_of_origin == "dad" else mom_snv_vcf
        str_vcf2query = dad_str_vcf if parent_of_origin == "dad" else mom_str_vcf

        # look in the STR VCF and ensure that the informative parent's STR 
        # is phased.
        str_ps = None
        for v in str_vcf2query(str_region):
            assert v.INFO.get("TRID") == row["trid"]
            
            str_ps_tag = v.format("PS")
            if str_ps_tag is None:
                continue
            str_ps = str_ps_tag[:, 0][0]

        if str_ps is None: continue

        # now, look in the SNV VCF and ensure that this informative site
        # shared a PS tag with the STR. also, figure out whether the A or B haplotype
        # in the informative parent is the one they donated to the kid.
        informative_haplotype = None
        for v in snv_vcf2query(snv_region):
            snv_ps_tag = v.format("PS")
            if snv_ps_tag is None:
                continue
            snv_ps = snv_ps_tag[:, 0][0]

            # ensure PS tags match. if the SNV and STR don't share
            # the same phase block, we can't be sure that e.g., at a 0|1
            # SNP and 1|2 STR, the 1 and 2 alleles were transmitted together
            if snv_ps != str_ps: continue

            hap_a, hap_b, is_phased = np.array(v.genotypes)[0]
            if not is_phased:
                continue
            
            # now there are four possible options
            # inf parent 0|1, other parent 0/0
            # inf parent 1|0, other parent 0/0
            # inf parent 0|1, other parent 1/1
            # inf parent 1|0, other parent 1/1
            # since the child is always 0/1 at these SNVs,
            # we know that...
            # ...at this informative site, if dad's genotype is larger
            # than mom's, then he must have donated the ALT allele. the informative
            # haplotype is therefore the haplotype with the ALT allele.
            if informative_parent_gt > other_parent_gt:
                informative_haplotype = "A" if hap_a > hap_b else "B"
            # if dad's genotype is less than mom's, he must have donated
            # the REF allele. the informative haplotype is therefore the
            # haplotype with the REF allele.
            elif informative_parent_gt < other_parent_gt:
                informative_haplotype = "A" if hap_a < hap_b else "B"

        if informative_haplotype is None: continue

        row_dict = row.to_dict()
        row_dict.update(
            {
                "haplotype_in_parent": informative_haplotype,
            }
        )
        res.append(row_dict)

    res_df = pd.DataFrame(res)

    res_df = informative_sites.merge(res_df, how="left").fillna(
        value={
            "haplotype_in_parent": "unknown",
        }
    )
    res_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--annotated_dnms",
        required=True,
        help="""TSV file containing DNM calls with one or more informative site annotations. \
            each DNM will be represented by as many rows as there were informative sites.""",
    )
    p.add_argument(
        "--dad_id",
        required=True,
        help="""LAB ID of father""",
    )
    p.add_argument("--dad_phased_str_vcf", required=True)
    p.add_argument("--dad_phased_snv_vcf", required=True)
    p.add_argument("--mom_phased_str_vcf", required=True)
    p.add_argument("--mom_phased_snv_vcf", required=True)
    p.add_argument(
        "--mom_id",
        required=True,
        help="""LAB ID of mother""",
    )
    p.add_argument(
        "--out",
        required=True,
        help="""name of TSV outfile""",
    )
    args = p.parse_args()
    main(args)
