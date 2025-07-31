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

CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()

# map sample names to ALT (Coriell) IDs, suffixes, and parents
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

ORIG_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/CHM13"
NEW_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025"

# map samples to list of sub-bams
# samples, seqruns, bcids = [], [], []
SMP2SEQRUNS, SMP2BCIDS = defaultdict(list), defaultdict(list)
for fh in glob.glob(f"{NEW_PATH}/K**/*.bam"):
    sample = fh.split("/")[-2].split("_")[1]
    seqrun = fh.split("/")[-1].split(".")[0]
    bcid = fh.split("/")[-1].split(".")[-2]
    # samples.append(sample)
    # seqruns.append(seqrun)
    # bcids.append(bcid)
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

CHILDREN = ["2211"]

def get_raw_ccs_bams(wildcards):
    seqruns = SMP2SEQRUNS[wildcards.SAMPLE]
    bcids = SMP2BCIDS[wildcards.SAMPLE]
    o = []
    for s,b in zip(seqruns, bcids):
        o.append(f"data/bam/raw/{wildcards.SAMPLE}.{s}.{b}.sorted.bam")
    return o


def get_complete_bams(wildcards, sample):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """
    # if wildcards.SAMPLE in ("200100", "2189", "2216", "200080"):
    #     return NEW_PATH + f"/merged/{wildcards.SAMPLE}.merged.bam"
    if sample in ("200100", "2189", "2216", "200080"):
        return NEW_PATH + f"/merged/{sample}.merged.bam"
    else:
        assembly_adj = wildcards.ASSEMBLY.split('v2')[0]
        #sample_id = SMP2ALT[wildcards.SAMPLE]
        sample_id = SMP2ALT[sample]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in wildcards.ASSEMBLY else assembly_adj,)


wildcard_constraints:
    SAMPLE = "[0-9]+",
    ASSEMBLY = "CHM13v2",
    CHROM = "chr[0-9]{1,2}|chrX|chrY"


rule all:
    input:
        dnms = expand("csv/{SAMPLE}.{ASSEMBLY}.trgt-denovo.csv", SAMPLE=CHILDREN, ASSEMBLY=["CHM13v2"]),
        vcf = expand("data/trgt/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz", SAMPLE=ALL_SAMPLES, ASSEMBLY=["CHM13v2"]),
        trio_vcf = expand("data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi", SAMPLE=CHILDREN, ASSEMBLY=["CHM13v2"]),
        phased_snv = expand("data/vcf/snv/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz", SAMPLE=CHILDREN, ASSEMBLY=["CHM13v2"])


rule align_ccs_bam:
    input:
        ref = REF_FH,
        bam = NEW_PATH + "/K1463_{SAMPLE}/{SEQRUN}.hifi_reads.{BCID}.bam"
    output:
        bam = "data/bam/raw/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam",
        bai = "data/bam/raw/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam.bai",
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
        new_bams = lambda wildcards: get_raw_ccs_bams(wildcards),
        old_bam = lambda wildcards: smp2orig[wildcards.SAMPLE]
    output:
        bam_fh = NEW_PATH + "/merged/{SAMPLE}.merged.bam",
    threads: 8
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

rule create_chrom_bed:
    input:
        repeats = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY]
    output:
        "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed.gz"
    shell:
        """
        zgrep -w {wildcards.CHROM} {input.repeats} | bgzip > {output}
        """


rule run_trgt:
    input:
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = REF_FH,
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v3.0.0-x86_64-unknown-linux-gnu/trgt",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed.gz"
    output:
        vcf = temp("data/trgt/per-chrom/{CHROM}.{SAMPLE}.{ASSEMBLY}.vcf.gz"),
        bam = temp("data/trgt/per-chrom/{CHROM}.{SAMPLE}.{ASSEMBLY}.spanning.bam")
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        {input.trgt_binary} genotype --threads 8 \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            --output-prefix "data/trgt/per-chrom/{wildcards.CHROM}.{wildcards.SAMPLE}.{wildcards.ASSEMBLY}"
        """

rule merge_trgt_bams:
    input:
        bams = expand("data/trgt/per-chrom/{CHROM}.{{SAMPLE}}.{{ASSEMBLY}}.spanning.bam", CHROM=CHROMS)
    output:
        temp("data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.bam")
    threads: 8
    shell:
        """
        module load samtools
        
        samtools merge -O BAM -o {output} -@ {threads} {input.bams}
        """


rule merge_trgt_vcfs:
    input:
        vcfs = expand("data/trgt/per-chrom/{CHROM}.{{SAMPLE}}.{{ASSEMBLY}}.vcf.gz", CHROM=CHROMS)
    output:
        "data/trgt/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -t ^chrEBV,phix | grep -v 'ID=chrEB\|ID=phix' | bgzip > {output}
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
        idx = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools sort -Ob -o {output.vcf} {input}
        bcftools index --tbi {output.vcf}
        """

# https://github.com/PacificBiosciences/wdl-common/blob/2ae390d4ed6b80dd2a2ef10c960832ffa8c7d1d3/wdl/workflows/deepvariant/deepvariant.wdl
# rule call_snvs:
#     input:
#         ref = REF_FH,
#         bam = lambda wildcards: get_complete_bams(wildcards),
#         bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
#         sif = "deepvariant_1.9.0.sif"
#     output:
#         vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
#         gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#     threads: 8
#     params:
#         alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
#     resources:
#         mem_mb = 64_000,
#         slurm_account = "quinlan-rw",
#         slurm_partition = "quinlan-rw",
#     shell:
#         """
#         module load singularity

#         export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

#         singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
#                                             {input.sif} \
#                                             run_deepvariant \
#                                                 --model_type PACBIO \
#                                                 --num_shards {threads} \
#                                                 --output_vcf {output.vcf} \
#                                                 --output_gvcf {output.gvcf} \
#                                                 --reads {input.bam} \
#                                                 --ref {input.ref} \
#                                                 --regions {wildcards.CHROM} \
#                                                 --make_examples_extra_args "vsc_min_fraction_indels=0.5,min_mapping_quality=1" \
#                                                 --sample_name {params.alt_sample_id}
#         """

rule call_trio_snvs:
    input:
        ref = REF_FH,
        kid_bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
        kid_bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
        dad_bam = lambda wildcards: get_complete_bams(wildcards, SMP2DAD[wildcards.SAMPLE]),
        dad_bam_idx = lambda wildcards: get_complete_bams(wildcards, SMP2DAD[wildcards.SAMPLE]) + ".bai",
        mom_bam = lambda wildcards: get_complete_bams(wildcards, SMP2MOM[wildcards.SAMPLE]),
        mom_bam_idx = lambda wildcards: get_complete_bams(wildcards, SMP2MOM[wildcards.SAMPLE]) + ".bai",
        sif = "deepvariant_deeptrio-1.9.0.sif"
    output:
        kid_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.kid.{ASSEMBLY}.{CHROM}.vcf.gz",
        kid_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.kid.{ASSEMBLY}.{CHROM}.g.vcf.gz",
        dad_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.dad.{ASSEMBLY}.{CHROM}.vcf.gz",
        dad_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.dad.{ASSEMBLY}.{CHROM}.g.vcf.gz",
        mom_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.mom.{ASSEMBLY}.{CHROM}.vcf.gz",
        mom_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.mom.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        kid_name = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
        dad_name = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_name = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
    resources:
        mem_mb = 64_000,
        slurm_account = "quinlan-rw",
        slurm_partition = "quinlan-rw",
    script:
        "scripts/run_deeptrio.sh"


# create sample-level VCFs that we'll use for hiphase
rule combine_sample_vcfs:
    input:
        vcfs = expand("data/vcf/snv/per-chrom/{{SAMPLE}}.kid.{{ASSEMBLY}}.{CHROM}.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -f 'PASS' -t ^chrEBV,phix | grep -v 'ID=chrEB\|ID=phix' | bgzip > {output}
        """


rule index_sample_vcfs:
    input: "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    output: "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi"
    threads: 4
    shell:
        """
        module load bcftools
        
        bcftools index --tbi --threads {threads} {input}
        """


# run hiphase on sample-level VCFs
rule run_hiphase:
    input:
        snv_vcf = "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz",
        snv_vcf_idx = "data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi",
        str_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz",
        str_vcf_idx = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz.tbi",
        bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
        bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
        reference = REF_FH
    output:
        snv_vcf = "data/vcf/snv/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz",
        str_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.phased.vcf.gz",
    threads: 4
    resources:
        mem_mb = 32_000
    shell:
        """
        ~/bin/hiphase-v1.4.4-x86_64-unknown-linux-gnu/hiphase --threads 4 \
                      --bam {input.bam} \
                      --vcf {input.snv_vcf} \
                      --vcf {input.str_vcf} \
                      --reference {input.reference} \
                      --output-vcf {output.snv_vcf} \
                      --output-vcf {output.str_vcf} \

        """

# joint-genotype the individual chromosomes in each trio we'll
# use to do informative site finding (no need to phase these)
rule combine_trio_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        kid_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.kid.{ASSEMBLY}.{CHROM}.g.vcf.gz",
        dad_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.dad.{ASSEMBLY}.{CHROM}.g.vcf.gz",
        mom_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.mom.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    output: "data/vcf/trios/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.trio.vcf.gz"
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
                                            --dir gl_nexus_dbs/{wildcards.SAMPLE}_{wildcards.ASSEMBLY}_{wildcards.CHROM} \
                                            --config DeepVariant_unfiltered \
                                            --mem-gbytes 32 \
                                            --threads {threads} \
                                            {input.kid_snp_gvcf} \
                                            {input.mom_snp_gvcf} \
                                            {input.dad_snp_gvcf} > {output}
        """

rule merge_trio_vcfs:
    input:
        vcfs = expand("data/vcf/trios/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{CHROM}.trio.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -t ^chrEBV,phix | grep -v 'ID=chrEB\|ID=phix' | bgzip > {output}
        """

rule index_merged_vcf:
    input:
        "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz"
    output:
        "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


rule run_trgt_denovo:
    input:
        reference = REF_FH,
        repeats = "data/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = "data/trgt/{SAMPLE}.{ASSEMBLY}.sorted.vcf.gz",
        mom_vcf = lambda wildcards: "data/trgt/" + SMP2MOM[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".sorted.vcf.gz",
        dad_vcf = lambda wildcards: "data/trgt/" + SMP2DAD[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".sorted.vcf.gz",
        kid_bam = "data/trgt/{SAMPLE}.{ASSEMBLY}.spanning.sorted.bam",
        mom_bam = lambda wildcards: "data/trgt/" + SMP2MOM[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".spanning.sorted.bam",
        dad_bam = lambda wildcards: "data/trgt/" + SMP2DAD[wildcards.SAMPLE] + "." + wildcards.ASSEMBLY + ".spanning.sorted.bam",
    output:
        output = "csv/{SAMPLE}.{ASSEMBLY}.trgt-denovo.csv"
    threads: 32
    resources:
        mem_mb = 64_000,
        runtime = 1_920
    params:
        kid_pref = lambda wildcards: f"data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}",
        mom_pref = lambda wildcards: f"data/trgt/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}",
        dad_pref = lambda wildcards: f"data/trgt/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}",
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
