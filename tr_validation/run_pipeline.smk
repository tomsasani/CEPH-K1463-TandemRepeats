VCF_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/GRCh38/adotto/"
# VCF_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/GRCh38/dashnow/fd6baac/"
BAM_PREF = "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/complete/"
ORIG_BAM_PREF = "/scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00065564/A1291-220303-VAR-Masked_CEPH_WGS/UCGD/GRCh38/Data/PolishedCrams/"

SAMPLES = ["2188", "2216", "2189", "2209"]
VCF_FHS = [VCF_PREF + c for c in ("2188_DM_adotto_v02.sorted.vcf.gz", "2216_D_adotto_v02.sorted.vcf.gz", "2189_SF_adotto_v02.sorted.vcf.gz", "2209_SF_adotto_v02.sorted.vcf.gz",)]
# VCF_FHS = [VCF_PREF + c for c in ("2188_DM_GRCh38.sorted.vcf.gz", "2216_DM_GRCh38.sorted.vcf.gz", "2189_SF_GRCh38.sorted.vcf.gz", "2209_SF_GRCh38.sorted.vcf.gz",)]

BAM_FHS = [BAM_PREF + c for c in ("2188_1_merged_sort.bam", "2216_1_merged_sort.bam", "2189_1_merged_sort.bam", "2209_2_merged_sort.bam",)]
ORIG_BAM_FHS = [ORIG_BAM_PREF + c for c in ("2188.cram", "2216.cram", "2189.cram", "2209.cram",)]

SAMPLE2VCF = dict(zip(SAMPLES, VCF_FHS))
SAMPLE2BAM = dict(zip(SAMPLES, BAM_FHS))
SAMPLE2BAM["missing"] = None

SAMPLE2MOM = {"2189": "2188", "2216": "2188", "2188": "missing", "2209": "missing"}
SAMPLE2DAD = {"2189": "2209", "2216": "2209", "2188": "missing", "2209": "missing"}

rule all:
    input: 
        expand("tr_validation/csv/{sample}.dnm.read_support.csv", sample=["2189", "2216", "2188", "2209"]),
        expand("tr_validation/fig/{sample}.dnm.read_support.png", sample=["2189", "2216", "2188", "2209"]),
        # expand("tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv", sample=["2216", "2189"], genotype=["0_0", "0_1", "1_1"]),

        # expand("tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv", sample=["2216", "2189"], genotype=["0_0", "0_1", "1_1"])
        #expand("tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv", sample=["2216"], genotype=["1_1"])



rule validate_dnms:
    input:
        py_script = "tr_validation/validate_dnms.py",
        # dnms = "tr_validation/data/df_trgtdn_palladium_GRCh38_C0.parquet.gzip",
        dnms = "tr_validation/data/df_transmission_C5.parquet.gzip",

    output: "tr_validation/csv/{sample}.dnm.read_support.csv"
    params:
        vcf = lambda wildcards: SAMPLE2VCF[wildcards.sample],
        kid_bam = lambda wildcards: SAMPLE2BAM[wildcards.sample],
        mom_bam = lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]],
        dad_bam = lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]]
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
        py_script = "tr_validation/validate_inherited.py",
        mutations = "/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/GRCh38/dashnow/fd6baac/element_validation_consistent_GTs/consistent_{genotype}.element_trio.tsv"
    output: "tr_validation/csv/{sample}.{genotype}.inherited.read_support.csv"
    params:
        vcf = lambda wildcards: SAMPLE2VCF[wildcards.sample],
        kid_bam = lambda wildcards: SAMPLE2BAM[wildcards.sample],
        mom_bam = lambda wildcards: SAMPLE2BAM[SAMPLE2MOM[wildcards.sample]],
        dad_bam = lambda wildcards: SAMPLE2BAM[SAMPLE2DAD[wildcards.sample]]

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
        py_script = "tr_validation/plot_read_support.py",
        csv = "tr_validation/csv/{sample}.dnm.read_support.csv"
    output: "tr_validation/fig/{sample}.dnm.read_support.png"
    shell:
        """
        python {input.py_script} --csv {input.csv} --out {output}
        """