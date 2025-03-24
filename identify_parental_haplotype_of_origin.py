import pandas as pd
from cyvcf2 import VCF
import numpy as np
import tqdm
import argparse
from annotate_with_informative_sites import measure_consistency

from schema import InformativeSiteSchema

def main(args):

    informative_sites = pd.read_csv(args.annotated_dnms, sep="\t", dtype={"sample_id": str})
    InformativeSiteSchema.validate(informative_sites)

    # file handles for (p/m)aternal SNV and STR VCFs (phased)
    dad_str_vcf = VCF(f"{args.hiphase_pref}/{args.dad_id_orig}_{args.assembly}_50bp_merge.sorted.phased.vcf.gz")
    dad_snv_vcf = VCF(f"{args.hiphase_pref}/{args.dad_id}.{args.assembly}.deepvariant.glnexus.phased.vcf.gz")
    mom_str_vcf = VCF(f"{args.hiphase_pref}/{args.mom_id_orig}_{args.assembly}_50bp_merge.sorted.phased.vcf.gz")
    mom_snv_vcf = VCF(f"{args.hiphase_pref}/{args.mom_id}.{args.assembly}.deepvariant.glnexus.phased.vcf.gz")

    res = []

    # loop over every informative site
    for i, row in tqdm.tqdm(informative_sites.iterrows()):
        # figure out the inferred parent of origin for this STR
        parent_of_origin = row["str_parent_of_origin"]
        if parent_of_origin not in ("mom", "dad"): continue
        # figure out the informative genotype in the parent
        informative_parent_gt, other_parent_gt = None, None
        if parent_of_origin == "dad":
            informative_parent_gt = row["dad_inf_gt"]
            other_parent_gt = row["mom_inf_gt"]
        else:
            informative_parent_gt = row["mom_inf_gt"]
            other_parent_gt = row["dad_inf_gt"]

        # get the location of this informative SNP
        chrom, start, end = row["inf_chrom"], row["inf_pos"] - 1, row["inf_pos"]
        snv_region = f"{chrom}:{start}-{end}"
        # get the location of the TR locus
        str_region = f"{row['#chrom']}:{row['start']}-{row['end']}"

        snv_vcf2query = dad_snv_vcf if parent_of_origin == "dad" else mom_snv_vcf
        str_vcf2query = dad_str_vcf if parent_of_origin == "dad" else mom_str_vcf

        str_ps_tag, allele_lengths = None, None

        str_hap_to_al = {}
        for v in str_vcf2query(str_region):
            assert v.INFO.get("TRID") == row["trid"]
            # get PS tag for STR in the parent
            try:
                str_ps_tag = v.format("PS")[:, 0][0]
            except TypeError:
                continue

            allele_lengths = [len(v.REF)] + [len(a) for a in v.ALT]

            # keep track of the STR allele lengths on either haplotype
            # in the parent
            str_hap_a, str_hap_b, is_phased = np.array(v.genotypes)[0]
            # use haplotype indices as keys
            str_hap_to_al[0] = allele_lengths[str_hap_a]
            str_hap_to_al[1] = allele_lengths[str_hap_b]

        if str_ps_tag is None: 
            continue

        # figure out the phase block corresponding to this STR in the parent
        for v in snv_vcf2query(snv_region):
            try:
                snv_ps_tag = v.format("PS")[:, 0][0]
            except TypeError:
                continue
            if snv_ps_tag != str_ps_tag: continue
            # if the SNV PS tag matches the STR PS tag, we can use
            # it to infer the haplotype that the parent donated to the kid
            # on which the DNM occurred. importantly, we need to know the
            # informative genotypes in mom and dad. if dad's genotype is
            # bigger than mom's, we should figure out the haplotype that has
            # more ALT alleles. if dad's genotype is smaller, we need to figure
            # out which haplotype has fewer ALT alleles
            hap_a, hap_b, is_phased = np.array(v.genotypes)[0]
            if informative_parent_gt > other_parent_gt:
                informative_haplotype = 0 if hap_a > hap_b else 1
            elif informative_parent_gt < other_parent_gt:
                informative_haplotype = 0 if hap_a < hap_b else 1

            row_dict = row.to_dict()

            row_dict.update(
                {
                    "haplotype_in_parent": "A" if informative_haplotype == 0 else "B",
                    "allele_length_in_parent": str_hap_to_al[informative_haplotype],
                }
            )
            res.append(row_dict)

    res_df = pd.DataFrame(res)

    res_df = informative_sites.merge(res_df, how="left").fillna(
        value={
            "haplotype_in_parent": "UNK",
            "allele_length_in_parent": -1,
        }
    )

    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--hiphase_pref")
    p.add_argument("--dad_id")
    p.add_argument("--dad_id_orig")
    p.add_argument("--mom_id")
    p.add_argument("--mom_id_orig")
    p.add_argument("--assembly")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
