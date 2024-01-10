import pandas as pd

SAMPLE_MANIFEST = pd.read_csv("tr_validation/data/file_mapping.tsv", sep="\t", dtype={"sample_id": str, "paternal_id": str, "maternal_id": str})

BAM_TECH = "element" # illumina, element, ont

SAMPLE2VCF = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["vcf_fh"]))
SAMPLE2BAM = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST[f"{BAM_TECH}_bam_fh"]))
SAMPLE2MOM = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["maternal_id"]))
SAMPLE2DAD = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["paternal_id"]))

print (SAMPLE2MOM)

# define samples of interest
SAMPLES = ["2188", "2216", "2189", "2209"]


rule all:
    input:
        "tr_validation/validate_dnms.py",
        expand("tr_validation/fig/{sample}.dnm.read_support.png", sample=SAMPLES),
        #expand("tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv", sample=["2216", "2189"], genotype=["0_0", "0_1", "1_1"]),


rule validate_dnms:
    input:
        py_script="tr_validation/validate_dnms.py",
        dnms="tr_validation/data/df_trgtdn_palladium_GRCh38_C0.parquet.gzip",
        #dnms = "tr_validation/data/df_transmission_C5.parquet.gzip",
    output:
        "tr_validation/csv/{sample}.dnm.read_support.csv",
    params:
        vcf=lambda wildcards: SAMPLE2VCF[wildcards.sample],
        kid_bam=lambda wildcards: SAMPLE2BAM[wildcards.sample],
        mom_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]] if SAMPLE2MOM[wildcards.sample] != "missing" else None,
        dad_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]] if SAMPLE2DAD[wildcards.sample] != "missing" else None,
    shell:
        """
        python {input.py_script} --dnms {input.dnms} \
                                 --vcf {params.vcf} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.sample} \
                                 --out {output} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                                 
        """


rule validate_inherited:
    input:
        py_script="tr_validation/validate_inherited.py",
        mutations="/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/GRCh38/dashnow/fd6baac/element_validation_consistent_GTs/consistent_{genotype}.element_trio.tsv",
    output:
        "tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv",
    params:
        vcf=lambda wildcards: SAMPLE2VCF[wildcards.sample],
        kid_bam=lambda wildcards: SAMPLE2BAM[wildcards.sample],
        mom_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]],
        dad_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]],
    shell:
        """
        python {input.py_script} --mutations {input.mutations} \
                                 --vcf {params.vcf} \
                                 --kid_bam {params.kid_bam} \
                                 --mom_bam {params.mom_bam} \
                                 --dad_bam {params.dad_bam} \
                                 --sample_id {wildcards.sample} \
                                 --out {output} 
        """


rule plot_support:
    input:
        py_script="tr_validation/plot_read_support.py",
        csv="tr_validation/csv/{sample}.dnm.read_support.csv",
    output:
        "tr_validation/fig/{sample}.dnm.read_support.png",
    shell:
        """
        python {input.py_script} --csv {input.csv} --out {output}
        """
