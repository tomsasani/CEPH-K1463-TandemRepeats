import glob
from collections import defaultdict
import pandas as pd


ASSEMBLY2CATALOG = {"GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz"}

REF_FH = "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/CHM13v2.0/primary_assembly_decoy_phix.fa"

# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))

CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()

# map sample names to ALT (Coriell) IDs, suffixes, and parents
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

ORIG_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/CHM13"
NEW_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025"

# map samples to list of sub-bams
samples, seqruns, bcids = [], [], []
SMP2SEQRUNS, SMP2BCIDS = defaultdict(list), defaultdict(list)
for fh in glob.glob(f"{NEW_PATH}/K**/*.bam"):
    sample = fh.split("/")[-2].split("_")[1]
    seqrun = fh.split("/")[-1].split(".")[0]
    bcid = fh.split("/")[-1].split(".")[-2]
    samples.append(sample)
    seqruns.append(seqrun)
    bcids.append(bcid)
    SMP2SEQRUNS[sample].append(seqrun)
    SMP2BCIDS[sample].append(bcid)

# add orig
smp2orig = {}
for sample in ("200080", "200100", "2189", "2216"):
    smp2orig[sample] = f"{ORIG_PATH}/{SMP2ALT[sample]}.chm13.haplotagged.bam"

ALL_SAMPLES = ped["sample_id"].unique()
CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

def get_raw_ccs_bams(wildcards):
    seqruns = SMP2SEQRUNS[wildcards.SAMPLE]
    bcids = SMP2BCIDS[wildcards.SAMPLE]
    o = []
    for s,b in zip(seqruns, bcids):
        o.append(f"data/bam/raw/{wildcards.SAMPLE}.{s}.{b}.sorted.bam")
    return o

def get_complete_bams(wildcards):
    if wildcards.SAMPLE in ("200100", "2189", "2216", "200080"):
        return NEW_PATH + f"/merged/{wildcards.SAMPLE}.merged.bam"
    else:
        assembly_adj = wildcards.ASSEMBLY.split('v2')[0]
        sample_id = SMP2ALT[wildcards.SAMPLE]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in wildcards.ASSEMBLY else assembly_adj,)

wildcard_constraints:
    SAMPLE = "[0-9]+",
    ASSEMBLY = "CHM13v2",
    CHROM = "chr[0-9]{1,2}|chrX|chrY"

rule all:
    input:
        dnms = expand("csv/{SAMPLE}.{ASSEMBLY}.trgt-denovo.csv", SAMPLE=CHILDREN[:2], ASSEMBLY=["CHM13v2"]),
        vcf = expand("data/trgt/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz", SAMPLE=CHILDREN[:2], ASSEMBLY=["CHM13v2"])

rule align_ccs_bam:
    input:
        ref = REF_FH,
        bam = NEW_PATH + "/K1463_{SAMPLE}/{SEQRUN}.hifi_reads.{BCID}.bam"
    output:
        bam = "data/bam/raw/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam",
        bai = "data/bam/raw/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam.bai",
    threads: 8
    shell:
        """
        ~/bin/pbmm2 align \
                    --num-threads {threads} \
                    --sort-memory 4G \
                    --preset CCS \
                    --sample {wildcards.SAMPLE} \
                    --sort \
                    --bam-index BAI \
                    --unmapped \
                    {input.ref} \
                    {input.bam} \
                    {output.bam} \
    """


rule combine_sample_bams:
    input:
        new_bams = lambda wildcards: get_raw_ccs_bams(wildcards),
        old_bam = lambda wildcards: smp2orig[wildcards.SAMPLE]
    output:
        bam_fh = NEW_PATH + "/merged/{SAMPLE}.merged.bam",
    threads: 4
    shell:
        """
        module load samtools
        
        samtools merge --threads {threads} -O BAM -o {output.bam_fh} {input.new_bams} {input.old_bam}
        """

rule index_sample_bams:
    input:
        bam_fh = NEW_PATH + "/merged/{SAMPLE}.merged.bam"
    output:
        bam_idx = NEW_PATH + "/merged/{SAMPLE}.merged.bam.bai",
    threads: 4
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam_fh}
        """

rule run_trgt:
    input:
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = REF_FH,
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v3.0.0-x86_64-unknown-linux-gnu/trgt",
        repeats = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY]
    output:
        vcf = temp("data/trgt/{SAMPLE}.{ASSEMBLY}.vcf.gz"),
        bam = temp("data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.bam")
    threads: 8
    shell:
        """
        {input.trgt_binary} genotype --threads 8 \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            --output-prefix "data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}"
        """



rule sort_bam:
    input: "data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.bam"
    output: "data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.sorted.bam"    
    shell:
        """
        module load samtools
        
        samtools sort -o {output} {input}
        samtools index {output}
        """

rule sort_vcf:
    input: "data/trgt/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    output:
        vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz",
        idx = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz.csi"
    shell:
        """
        module load bcftools
        
        bcftools sort -Ob -o {output.vcf} {input}
        bcftools index {output.vcf}
        """


rule call_snvs:
    input:
        ref = REF_FH,
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf",
        gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.gvcf",
    threads: 4
    resources:
        mem_mb = 64_000
    shell:
        """
        module load singularity

        export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

        singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
                                            {input.sif} \
                                            run_deepvariant \
                                                --model_type PACBIO \
                                                --num_shards 4 \
                                                --output_vcf {output.vcf} \
                                                --output_gvcf {output.gvcf} \
                                                --reads {input.bam} \
                                                --ref {input.ref} \
                                                --regions {wildcards.CHROM} \
                                                --sample_name {wildcards.SAMPLE}
        """


rule merge_chrom_vcfs:
    input:
        vcfs = expand("data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{CHROM}.vcf", CHROM=CHROMS)
    output: "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz"    
    shell:
        """
        module load bcftools
        
        bcftools concat -Oz -o {output} {input.vcfs}
        """

rule index_merged_vcf:
    input:
        "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    output:
        "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


rule run_hiphase:
    input:
        snv_vcf = "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz",
        snv_vcf_idx = "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi",
        str_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz",
        str_vcf_idx = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz.csi",
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = REF_FH
    output:
        snv_vcf = "data/vcf/snv/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz",
        str_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz",
    threads: 4
    shell:
        """
        ~/bin/hiphase-v1.4.4-x86_64-unknown-linux-gnu/hiphase --threads 4 \
                      --bam {input.bam} \
                      --vcf {input.snv_vcf} \
                      --vcf {input.str_vcf} \
                      --reference {input.reference} \
                      --output-vcf {output.snv_vcf} \
                      --output-vcf {output.str_vcf}

        """


rule run_trgt_denovo:
    input:
        reference = REF_FH,
        repeats = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY],
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz",
        mom_vcf = lambda wildcards: "data/trgt/" + SMP2MOM[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".sorted.vcf.gz",
        dad_vcf = lambda wildcards: "data/trgt/" + SMP2DAD[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".sorted.vcf.gz",
        kid_bam = "data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.bam",
        mom_bam = lambda wildcards: "data/trgt/" + SMP2MOM[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".spanning.sorted.bam",
        dad_bam = lambda wildcards: "data/trgt/" + SMP2DAD[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".spanning.sorted.bam",
    output:
        output = "csv/{SAMPLE}.{ASSEMBLY}.trgt-denovo.csv"
    threads: 8
    params:
        kid_pref = lambda wildcards: f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}",
        mom_pref = lambda wildcards: f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}",
        dad_pref = lambda wildcards: f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}",
    shell:
        """
        {input.trgt_denovo_binary} trio --reference {input.reference} \
                                        --bed {input.repeats} \
                                        --father {params.dad_pref} \
                                        --mother {params.mom_pref} \
                                        --child {params.kid_pref} \
                                        --out {output} \
                                        -@ 8
        """
