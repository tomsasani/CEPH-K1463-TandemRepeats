VCF_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/GRCh38/adotto/"
BAM_PREF = "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/complete/"

SAMPLES = ["2188", "2216", "2189", "2209"]
VCF_FHS = [VCF_PREF + c for c in ("2188_DM_adotto_v02.sorted.vcf.gz", "2216_D_adotto_v02.sorted.vcf.gz", "2189_SF_adotto_v02.sorted.vcf.gz", "2209_SF_adotto_v02.sorted.vcf.gz",)]
BAM_FHS = [BAM_PREF + c for c in ("2188_1_merged_sort.bam", "2216_1_merged_sort.bam", "2189_1_merged_sort.bam", "2209_2_merged_sort.bam",)]

SAMPLE2VCF = dict(zip(SAMPLES, VCF_FHS))
SAMPLE2BAM = dict(zip(SAMPLES, BAM_FHS))

rule all:
    input: expand("tr_validation/fig/{sample}.read_support.png", sample=SAMPLES)


rule run_iterator:
    input:
        py_script = "tr_validation/vcf_iterator.py",
        dnms = "tr_validation/data/df_transmission_C5.parquet.gzip"
    output: "tr_validation/csv/{sample}.read_support.csv"
    params:
        vcf = lambda wildcards: SAMPLE2VCF[wildcards.sample],
        bam = lambda wildcards: SAMPLE2BAM[wildcards.sample]
    shell:
        """
        python {input.py_script} --dnms {input.dnms} \
                                 --vcf {params.vcf} \
                                 --bam {params.bam} \
                                 --sample_id {wildcards.sample} \
                                 --out {output} \
        """

rule plot_support:
    input:
        py_script = "tr_validation/plot_read_support.py",
        csv = "tr_validation/csv/{sample}.read_support.csv"
    output: "tr_validation/fig/{sample}.read_support.png"
    shell:
        """
        python {input.py_script} --csv {input.csv} --out {output}
        """