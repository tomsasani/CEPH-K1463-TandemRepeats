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
ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )
SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))

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
    SMP2SEQRUNS[sample].append(seqrun)
    SMP2BCIDS[sample].append(bcid)

# add orig
smp2orig = {}
for sample in ("200080", "200100", "2189", "2216"):
    smp2orig[sample] = f"{ORIG_PATH}/{SMP2ALT[sample]}.chm13.haplotagged.bam"


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
    if sample in ("200100", "2189", "2216", "200080"):
        return NEW_PATH + f"/merged/{sample}.merged.bam"
    else:
        assembly_adj = wildcards.ASSEMBLY.split('v2')[0]
        sample_id = SMP2ALT[sample]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in wildcards.ASSEMBLY else assembly_adj,)


wildcard_constraints:
    SAMPLE = "[0-9]+",
    ASSEMBLY = "CHM13v2",
    CHROM = "chr[0-9]{1,2}|chrX|chrY",
    VERSION = "v[0-9].[0-9]"

CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

ALL_SAMPLES = ped["sample_id"].unique()
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
VERSIONS = ["v4.0"]

rule all:
    input:
        dnms = expand("csv/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.trgt-denovo.csv", CHROM=CHROMS, SAMPLE=CHILDREN, ASSEMBLY=["CHM13v2"], VERSION=VERSIONS),
        phased_trgt = expand("data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz.tbi", SAMPLE=ALL_SAMPLES, ASSEMBLY=["CHM13v2"], VERSION=VERSIONS),
        trio_vcf = expand("data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz.tbi", SAMPLE=CHILDREN, ASSEMBLY=["CHM13v2"]),
        phased_snv = expand("data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz.tbi", SAMPLE=ALL_SAMPLES, ASSEMBLY=["CHM13v2"], VERSION=VERSIONS)


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
        "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    shell:
        """
        zgrep -w {wildcards.CHROM} {input.repeats} > {output}
        """


rule run_trgt:
    input:
        bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
        bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
        reference = REF_FH,
        trgt_binary = lambda wildcards: "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v4.0.0-x86_64-unknown-linux-gnu/trgt" if wildcards.VERSION == "v4.0" else "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v0.7.0-linux_x86_64",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    output:
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.vcf.gz",
        bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.spanning.bam"
    params:
        genotype_cmd = lambda wildcards: "" if wildcards.VERSION == "v0.7" else "genotype",
        karyotype_cmd = lambda wildcards: "" if wildcards.VERSION == "v4.0" else "--karyotype XX" if SMP2SEX[wildcards.SAMPLE] == "female" else "--karyotype XY"
    threads: 16
    resources:
        mem_mb = 64_000
    shell:
        """
        {input.trgt_binary} {params.genotype_cmd} --threads {threads} \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            {params.karyotype_cmd} \
                            --output-prefix "data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}" \
                            -v
        """


rule sort_chrom_vcfs:
    input: "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.vcf.gz"
    output:
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz",
        idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        bcftools view -t ^chrEBV,phix {input} | grep -v 'ID=chrEB\|ID=phix' | \
        bcftools sort -Ob -o {output.vcf}

        bcftools index --tbi {output.vcf}
        """


rule sort_bam:
    input: "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.spanning.bam"
    output: "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.spanning.sorted.bam"
    shell:
        """
        module load samtools
        
        samtools sort -o {output} {input}
        samtools index {output}
        """

# https://github.com/PacificBiosciences/wdl-common/blob/2ae390d4ed6b80dd2a2ef10c960832ffa8c7d1d3/wdl/workflows/deepvariant/deepvariant.wdl
rule call_snvs:
    input:
        ref = REF_FH,
        bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
        bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    resources:
        mem_mb = 64_000,
    script:
        "scripts/run_deepvariant.sh"


# rule call_trio_snvs:
#     input:
#         ref = REF_FH,
#         kid_bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
#         kid_bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
#         dad_bam = lambda wildcards: get_complete_bams(wildcards, SMP2DAD[wildcards.SAMPLE]),
#         dad_bam_idx = lambda wildcards: get_complete_bams(wildcards, SMP2DAD[wildcards.SAMPLE]) + ".bai",
#         mom_bam = lambda wildcards: get_complete_bams(wildcards, SMP2MOM[wildcards.SAMPLE]),
#         mom_bam_idx = lambda wildcards: get_complete_bams(wildcards, SMP2MOM[wildcards.SAMPLE]) + ".bai",
#         sif = "deepvariant_deeptrio-1.9.0.sif"
#     output:
#         kid_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.kid.{ASSEMBLY}.{CHROM}.vcf.gz",
#         kid_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.kid.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#         dad_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.dad.{ASSEMBLY}.{CHROM}.vcf.gz",
#         dad_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.dad.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#         mom_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.mom.{ASSEMBLY}.{CHROM}.vcf.gz",
#         mom_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.mom.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#     threads: 16
#     params:
#         kid_name = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
#         dad_name = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
#         mom_name = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
#     resources:
#         mem_mb = 64_000,
#         slurm_account = "ucgd-rw",
#         slurm_partition = "ucgd-shared-rw",
#         runtime = 1440,
#     script:
#         "scripts/run_deeptrio.sh"


rule sort_snv_chrom_vcfs:
    input: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz"
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.sorted.vcf.gz"
    threads: 4
    shell:
        """
        module load bcftools

        bcftools view -t ^chrEBV,phix {input} | grep -v 'ID=chrEB\|ID=phix' | \
        bcftools sort -Ob -o {output}
        """


rule index_snv_chrom_vcfs:
    input: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.sorted.vcf.gz"
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.sorted.vcf.gz.tbi"
    threads: 4
    shell:
        """
        module load bcftools
        
        bcftools index --tbi --threads {threads} {input}
        """


# run hiphase on sample-level VCFs
rule run_hiphase:
    input:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.sorted.vcf.gz",
        snv_vcf_idx = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.sorted.vcf.gz.tbi",
        str_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz",
        str_vcf_idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz.tbi",
        bam = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE),
        bam_idx = lambda wildcards: get_complete_bams(wildcards, wildcards.SAMPLE) + ".bai",
        reference = REF_FH
    output:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.phased.vcf.gz",
        str_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.phased.vcf.gz",
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
        vcfs = expand("data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{VERSION}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -f 'PASS' -t ^chrEBV,phix | grep -v 'ID=chrEB\|ID=phix' | bgzip > {output}
        """

rule index_phased_snv_vcfs:
    input: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz"
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input}
        """


rule combine_phased_trgt_vcfs:
    input:
        vcfs = expand("data/trgt/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{VERSION}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -t ^chrEBV,phix | grep -v 'ID=chrEB\|ID=phix' | bgzip > {output}
        """

rule index_phased_trgt_vcfs:
    input: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz"
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz.tbi"
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
        kid_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
        dad_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.CHROM}.g.vcf.gz",
        mom_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.CHROM}.g.vcf.gz",
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

rule create_chrom_definitions:
    input:
        "data/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed",
    output:
        "data/catalogs/per-chrom/{CHROM}.CHM13v2.trgt.bed"
    shell:
        """
        grep -w {wildcards.CHROM} {input} > {output}
        """

rule run_trgt_denovo:
    input:
        reference = REF_FH,
        repeats =  "data/catalogs/per-chrom/{CHROM}.CHM13v2.trgt.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz",
        kid_vcf_idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.sorted.vcf.gz.tbi",
        mom_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.sorted.vcf.gz",
        mom_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        dad_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.sorted.vcf.gz",
        dad_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        kid_bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.spanning.sorted.bam",
        mom_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.spanning.sorted.bam",
        dad_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}.spanning.sorted.bam",
    output:
        output = "csv/{SAMPLE}.{ASSEMBLY}.{VERSION}.{CHROM}.trgt-denovo.csv"
    threads: 32
    resources:
        mem_mb = 32_000,
        runtime = 720,
        slurm_account = "ucgd-rw",
        slurm_partition = "ucgd-rw"
    params:
        kid_pref = lambda wildcards: f"data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}",
        mom_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}",
        dad_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.{wildcards.CHROM}",
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
