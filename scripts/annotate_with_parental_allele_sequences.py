import pandas as pd
from cyvcf2 import VCF

mutations = pd.read_csv(snakemake.input.annotated_dnms, sep="\t")

# read in VCF files from the trio. we read in each VCF separately because we need
# each to be phased with HiPhase, and our merged VCF is no longer phased
kid_vcf = VCF(snakemake.input.kid_phased_str_vcf)
mom_vcf = VCF(snakemake.input.mom_phased_str_vcf)
dad_vcf = VCF(snakemake.input.dad_phased_str_vcf)

res = []
for i, row in mutations.iterrows():
    trid, chrom, start, end = row["trid"], row["#chrom"], row["start"], row["end"]
    region = f"{chrom}:{start}-{end}"

    out_dict = row.to_dict()

    is_sex_chrom = row["#chrom"] in ("chrX", "chrY")

    # get the allele sequence of the de novo allele in the
    # child's VCF
    for v in kid_vcf(region):
        assert v.INFO.get("TRID") == trid
        alleles = [v.REF] + v.ALT
        denovo_gt = row["genotype"]
        # update our dictionary with the sequence and AP of the DNM
        out_dict.update({"denovo_allele_sequence": alleles[denovo_gt]})

    # depending on the parent of origin at this DNM, decide which
    # parental VCF we'll use to query for allele sequence information
    parent_of_origin, poi_support = row["phase_consensus"].split(":")

    if parent_of_origin == "dad":
        vcf2query = dad_vcf
    elif parent_of_origin == "mom":
        vcf2query = mom_vcf
    else:
        out_dict.update(
            {
                "precursor_allele_length_in_parent": -1,
                "precursor_sequence_in_parent": "unknown",
                "untransmitted_sequence_in_parent": "unknown",
                "precursor_AP": -1,
                "untransmitted_AP": -1,
            }
        )
        res.append(out_dict)
        continue

    # check to see if we were able to confidently infer which parental
    # haplotype the de novo mutation occurred on
    informative_haplotype, hoi_support = row["haplotype_in_parent_consensus"].split(":")

    if informative_haplotype == "unknown":
        out_dict.update(
            {
                "precursor_allele_length_in_parent": -1,
                "precursor_sequence_in_parent": "unknown",
                "untransmitted_sequence_in_parent": "unknown",
                "precursor_AP": -1,
                "untransmitted_AP": -1,
            }
        )
    else:
        # if we could, iterate over the parent's VCF
        for v in vcf2query(region):
            assert v.INFO.get("TRID") == trid
            allele_seqs = [v.REF] + v.ALT

            allele_lengths = [len(a) for a in allele_seqs]
            str_hap_a, str_hap_b, is_phased = v.genotypes[0]
            
            # assert is_phased
            if not is_phased: continue

            als = v.format("AL")[0]
            aps = v.format("AP")[0]

            # NOTE: THIS IS DANGEROUS WITHOUT CONFIRMING
            # THAT TRGT ALWAYS OUTPUTS GENOTYPES IN ASCENDING
            # ORDER!!!
            # get index into other FORMAT fields for each genotype,
            # assuming that the order of FORMAT values *didn't* change
            # after phasing and the unphased VCF always reported genotypes
            # in ascending order.
            if informative_haplotype == "A":
                informative_genotype, other_genotype = str_hap_a, str_hap_b
                fmt_idx = 0 if str_hap_a < str_hap_b else 1
            elif informative_haplotype == "B":
                informative_genotype, other_genotype = str_hap_b, str_hap_a
                fmt_idx = 0 if str_hap_b < str_hap_a else 1
            else:
                continue

            denovo_allele_sequence = allele_seqs[informative_genotype]
            denovo_allele_length = als[fmt_idx]
            
            # assert denovo_allele_length == allele_lengths[informative_genotype]
            denovo_allele_purity = aps[fmt_idx]

            if is_sex_chrom and parent_of_origin == "dad":
                non_denovo_allele_sequence = "UNK"
                non_denovo_allele_purity = -1
                non_denovo_allele_length = -1
            else:
                non_denovo_allele_sequence = allele_seqs[other_genotype]
                non_denovo_allele_length = als[1 - fmt_idx]
                # assert non_denovo_allele_length == allele_lengths[other_genotype]
                non_denovo_allele_purity = aps[1 - fmt_idx]

            out_dict.update(
                {
                    "precursor_allele_length_in_parent": denovo_allele_length,
                    "precursor_sequence_in_parent": denovo_allele_sequence,
                    "untransmitted_sequence_in_parent": non_denovo_allele_sequence,
                    "precursor_AP": denovo_allele_purity,
                    "untransmitted_AP": non_denovo_allele_purity,
                }
            )
    res.append(out_dict)

df = pd.DataFrame(res)
df = mutations.merge(df, how="left").fillna(
    {
        "precursor_allele_length_in_parent": -1,
        "precursor_sequence_in_parent": "unknown",
        "untransmitted_sequence_in_parent": "unknown",
        "precursor_AP": -1,
        "untransmitted_AP": -1,
    }
)
print (df.shape)
df.to_csv(snakemake.output.out, sep="\t", index=False)

