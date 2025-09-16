import glob
from collections import defaultdict
import pandas as pd


ASSEMBLY2CATALOG = {"GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
                    "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz"}

ASSEMBLY2REF = {"GRCh38": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix.fa", 
                "CHM13v2": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/CHM13v2.0/primary_assembly_decoy_phix.fa"}

# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))
ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )
SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))

# map sample names to ALT (Coriell) IDs, suffixes, and parents
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

ORIG_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams"
NEW_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025"

# map samples to list of sub-bams
# samples, seqruns, bcids = [], [], []
SMP2SEQRUNS, SMP2BCIDS = defaultdict(list), defaultdict(list)
for fh in glob.glob(f"{NEW_PATH}/K**/*.bam"):
    sample = fh.split("/")[-2].split("_")[1]
    seqrun = fh.split("/")[-1].split(".")[0]
    bcid = fh.split("/")[-1].split(".")[-2]
    SMP2SEQRUNS[sample].append(seqrun)
    SMP2BCIDS[sample].append(bcid)


def get_orig_sample_bam(wildcards):
    assembly_adj = wildcards.ASSEMBLY.lower().rstrip("v2") if wildcards.ASSEMBLY == "CHM13v2" else wildcards.ASSEMBLY
    return f"{ORIG_PATH}/{wildcards.ASSEMBLY.rstrip("v2")}/{SMP2ALT[wildcards.SAMPLE]}.{assembly_adj}.haplotagged.bam"
    

def get_raw_ccs_bams(wildcards):
    seqruns = SMP2SEQRUNS[wildcards.SAMPLE]
    bcids = SMP2BCIDS[wildcards.SAMPLE]
    o = []
    for s,b in zip(seqruns, bcids):
        o.append(f"data/bam/raw/{wildcards.ASSEMBLY}/{wildcards.SAMPLE}.{s}.{b}.sorted.bam")
    return o


def get_complete_bams(wildcards):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """
    if wildcards.USE_NEW_BAM == "TOPUP" and wildcards.SAMPLE in ("200100", "2189", "2216", "200080"):
        return f"{NEW_PATH}/merged/{wildcards.ASSEMBLY}/{wildcards.SAMPLE}.merged.bam"
    else:
        assembly_adj = wildcards.ASSEMBLY.split('v2')[0]
        sample_id = SMP2ALT[wildcards.SAMPLE]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in wildcards.ASSEMBLY else assembly_adj,)


wildcard_constraints:
    SAMPLE = "[0-9]+",
    ASSEMBLY = "CHM13v2|GRCh38",
    CHROM = "chr[0-9]{1,2}|chrX|chrY",
    USE_NEW_BAM = "TOPUP|ORIGINAL"


CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])
# CHROMS = ["chr1"]

ALL_SAMPLES = ped["sample_id"].unique()
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
USE_NEW_BAMS = ["ORIGINAL"]
ASSEMBLIES = ["GRCh38", "CHM13v2"]

CHILDREN = [c for c in CHILDREN if c.startswith("200")]

# CHILDREN = ["200084", "200101"]
# ALL_SAMPLES = ["2216", "2189", "200100", "200080"]

rule all:
    input:
        dnms = expand("csv/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trgt-denovo.csv", CHROM=CHROMS, SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        phased_trgt = expand("data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        trio_vcf = expand("data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz.tbi", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        phased_snv = expand("data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS)


rule align_ccs_bam:
    input:
        ref = lambda wildcards: ASSEMBLY2REF[wildcards.ASSEMBLY],
        bam = NEW_PATH + "/K1463_{SAMPLE}/{SEQRUN}.hifi_reads.{BCID}.bam"
    output:
        bam = "data/bam/raw/{ASSEMBLY}/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam",
        bai = "data/bam/raw/{ASSEMBLY}/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam.bai",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    resources:
        mem_mb = 64_000
    shell:
        """
        ~/bin/pbmm2 align \
                    --num-threads {threads} \
                    --sort \
                    --sort-memory 4G \
                    --preset CCS \
                    --sample {params.alt_sample_id} \
                    --bam-index BAI \
                    --unmapped \
                    {input.ref} \
                    {input.bam} \
                    {output.bam}
    """


rule combine_sample_bams:
    input:
        new_bams = get_raw_ccs_bams,
        old_bam = get_orig_sample_bam
    output:
        bam_fh = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam",
    threads: 8
    shell:
        """
        module load samtools
        
        samtools merge --threads {threads} -O BAM -o {output.bam_fh} {input.new_bams} {input.old_bam}
        """


rule index_sample_bams:
    input:
        bam_fh = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam"
    output:
        bam_idx = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam.bai",
    threads: 4
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam_fh}
        """


rule create_chrom_bed:
    input:
        repeats = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY]
    output:
        "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    shell:
        """
        zgrep -w {wildcards.CHROM} {input.repeats} > {output}
        """


rule create_chrom_ref:
    input:
        reference = lambda wildcards: ASSEMBLY2REF[wildcards.ASSEMBLY],
    output: contig = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz"
    shell:
        """
        module load samtools
        
        samtools faidx {input.reference} {wildcards.CHROM} | bgzip > {output.contig}
        """


rule index_chrom_ref:
    input:
        contig = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz"
    output: idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai"
    shell:
        """
        module load samtools
        
        samtools faidx {input.contig}
        """


rule run_trgt:
    input:
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v4.0.0-x86_64-unknown-linux-gnu/trgt",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    output:
        vcf = temp("data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz"),
        bam = temp("data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.bam")
    params:
        karyotype_cmd = lambda wildcards: "--karyotype XY" if SMP2SEX[wildcards.SAMPLE] == "male" and wildcards.CHROM in ("chrX", "chrY") else ""
    threads: 16
    resources:
        mem_mb = 64_000
    shell:
        """
        {input.trgt_binary} genotype --threads {threads} \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            {params.karyotype_cmd} \
                            --output-prefix "data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}" \
                            -v
        """


rule sort_chrom_vcfs:
    input: 
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
    output:
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        bcftools view -t {wildcards.CHROM} {input.vcf} | grep 'fileformat\|FILTER\|INFO\|FORMAT\|trgt\|bcftools\|{wildcards.CHROM},\|^{wildcards.CHROM}' | \
        bcftools sort -Ob -o {output.vcf}
        bcftools index --tbi {output.vcf}
        """


rule sort_bam:
    input: "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.bam"
    output: temp("data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam")
    shell:
        """
        module load samtools
        
        samtools sort -o {output} {input}
        samtools index {output}
        """


rule call_snvs:
    input:
        ref = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        ref_idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai",
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    resources:
        mem_mb = 64_000,
    script:
        "scripts/run_deepvariant.sh"


rule sort_snv_chrom_vcfs:
    input: 
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    threads: 4
    shell:
        """
        module load bcftools

        bcftools sort -Ob -o {output} {input.vcf}
        """


rule index_snv_chrom_vcfs:
    input: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi"
    threads: 4
    shell:
        """
        module load bcftools
        
        bcftools index --tbi --threads {threads} {input}
        """


# run hiphase on sample-level VCFs
rule run_hiphase:
    input:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        snv_vcf_idx = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        str_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        str_vcf_idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
    output:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.phased.vcf.gz",
        str_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.phased.vcf.gz",
    threads: 4
    resources:
        mem_mb = 32_000
    shell:
        """
        ~/bin/hiphase-v1.4.4-x86_64-unknown-linux-gnu/hiphase --threads {threads} \
                      --bam {input.bam} \
                      --vcf {input.snv_vcf} \
                      --vcf {input.str_vcf} \
                      --reference {input.reference} \
                      --output-vcf {output.snv_vcf} \
                      --output-vcf {output.str_vcf} \

        """


rule combine_phased_snv_vcfs:
    input:
        vcfs = expand("data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -f 'PASS' | bgzip > {output}
        """


rule index_phased_snv_vcfs:
    input: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


rule combine_phased_trgt_vcfs:
    input:
        vcfs = expand("data/trgt/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_phased_trgt_vcfs:
    input: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


# joint-genotype the individual chromosomes in each trio we'll
# use to do informative site finding (no need to phase these)
rule combine_trio_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        kid_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.g.vcf.gz",
        dad_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.g.vcf.gz",
        mom_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.g.vcf.gz",
    output: "data/vcf/trios/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trio.vcf.gz"
    threads: 4
    resources:
        mem_mb = 32_000
    shell:
        """
        module load singularity

        export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

        singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
                                            {input.sif} \
                                            /usr/local/bin/glnexus_cli \
                                            --dir gl_nexus_dbs/{wildcards.SAMPLE}_{wildcards.ASSEMBLY}_{wildcards.USE_NEW_BAM}_{wildcards.CHROM} \
                                            --config DeepVariant_unfiltered \
                                            --mem-gbytes 32 \
                                            --threads {threads} \
                                            {input.kid_snp_gvcf} \
                                            {input.mom_snp_gvcf} \
                                            {input.dad_snp_gvcf} > {output}
        """


rule merge_trio_vcfs:
    input:
        vcfs = expand("data/vcf/trios/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.trio.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_merged_vcf:
    input: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz"
    output: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


rule run_trgt_denovo:
    input:
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        repeats =  "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        kid_vcf_idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        mom_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        mom_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        dad_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        dad_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        kid_bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam",
        mom_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
        dad_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
    output:
        output = "csv/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trgt-denovo.csv"
    threads: 32
    resources:
        mem_mb = 32_000,
        runtime = 720,
        slurm_account = "ucgd-rw",
        slurm_partition = "ucgd-rw"
    params:
        kid_pref = lambda wildcards: f"data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        mom_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        dad_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
    shell:
        """
        {input.trgt_denovo_binary} trio --reference {input.reference} \
                                        --bed {input.repeats} \
                                        --father {params.dad_pref} \
                                        --mother {params.mom_pref} \
                                        --child {params.kid_pref} \
                                        --out {output} \
                                        -@ {threads}
        """
