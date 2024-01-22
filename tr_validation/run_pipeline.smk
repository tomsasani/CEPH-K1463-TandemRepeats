import pandas as pd

SAMPLE_MANIFEST = pd.read_csv("tr_validation/data/file_mapping.tsv", sep="\t", dtype={"sample_id": str, "paternal_id": str, "maternal_id": str})

BAM_TECH = "element" # illumina, element, ont

SAMPLE2VCF = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["vcf_fh"]))
SAMPLE2BAM = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST[f"{BAM_TECH}_bam_fh"]))
SAMPLE2MOM = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["maternal_id"]))
SAMPLE2DAD = dict(zip(SAMPLE_MANIFEST["sample_id"], SAMPLE_MANIFEST["paternal_id"]))
INHERITANCE2FILE = {"dnm": "tr_validation/data/df_transmission_C0_T0.parquet.gzip", 
                    "inherited": "/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/GRCh38/dashnow/fd6baac/element_validation_consistent_GTs/consistent_1_1.element_trio.tsv"}


# define samples of interest
SAMPLES = ["2188", "2216", "2189", "2209"]
SAMPLES = ["2216", "2189"]
KINDS = ["dnm", "inherited"]

rule all:
    input:
        "tr_validation/validate_with_bams.py",
        expand("tr_validation/fig/{sample}.{tech}.{kind}.read_support.png", sample=SAMPLES, tech=[BAM_TECH], kind=KINDS),


rule validate_mutations:
    input:
        py_script="tr_validation/validate_with_bams.py",
    output:
        "tr_validation/csv/{sample}.{tech}.{kind}.read_support.csv",
    params:
        vcf=lambda wildcards: SAMPLE2VCF[wildcards.sample],
        mutations=lambda wildcards: INHERITANCE2FILE[wildcards.kind],
        kid_bam=lambda wildcards: SAMPLE2BAM[wildcards.sample],
        mom_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]] if SAMPLE2MOM[wildcards.sample] != "missing" else None,
        dad_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]] if SAMPLE2DAD[wildcards.sample] != "missing" else None,
        max_al=lambda wildcards: 100_000 if wildcards.tech == "ont" else 120,

    shell:
        """
        python {input.py_script} --mutations {params.mutations} \
                                 --vcf {params.vcf} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.sample} \
                                 --out {output} \
                                 --variant_type {wildcards.kind} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                                 -max_al {params.max_al} \                                 
        """


# rule validate_inherited:
#     input:
#         py_script="tr_validation/validate_inherited.py",
#         mutations="/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/GRCh38/dashnow/fd6baac/element_validation_consistent_GTs/consistent_1_1.element_trio.tsv",
#     output:
#         "tr_validation/csv/{sample}.{tech}.inherited.read_support.csv",
#     params:
#         vcf=lambda wildcards: SAMPLE2VCF[wildcards.sample],
#         kid_bam=lambda wildcards: SAMPLE2BAM[wildcards.sample],
#         mom_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]] if SAMPLE2MOM[wildcards.sample] != "missing" else None,
#         dad_bam=lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]] if SAMPLE2DAD[wildcards.sample] != "missing" else None,
#     shell:
#         """
#         python {input.py_script} --mutations {input.mutations} \
#                                  --vcf {params.vcf} \
#                                  --kid_bam {params.kid_bam} \
#                                  --mom_bam {params.mom_bam} \
#                                  --dad_bam {params.dad_bam} \
#                                  --sample_id {wildcards.sample} \
#                                  --out {output} 
#         """


rule plot_support:
    input:
        py_script="tr_validation/plot_read_support.py",
        csv="tr_validation/csv/{sample}.{tech}.{kind}.read_support.csv",
    output:
        "tr_validation/fig/{sample}.{tech}.{kind}.read_support.png",
    shell:
        """
        python {input.py_script} --csv {input.csv} --out {output}
        """
