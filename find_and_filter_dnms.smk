import pandas as pd
from typing import Dict
from collections import Counter
import glob

ASSEMBLY2CATALOG = {"GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
                    "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz"}


AWS_PREF = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws"
CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/"

HPRC_SAMPLES = []
for fh in glob.glob(HPRC_PREF + "GRCh38_v1.0_50bp_merge/1.1.2-69937d83/hprc/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)

def get_transmission_fh(wildcards):
    assembly_adj = wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY
    suffix = SMP2SUFF[wildcards.SAMPLE]

    if wildcards.SAMPLE not in ("2216", "2189", "2209", "2188"):
        return "UNK"
    else:
        return f"{AWS_PREF}/{wildcards.ASSEMBLY}_v1.0_50bp_merge/493ef25/trgt-denovo/transmission/{wildcards.SAMPLE}_{suffix}_{assembly_adj}_50bp_merge_trgtdn.transmission.csv.gz"


def get_annotation_fh(wildcards):
    if wildcards.ASSEMBLY == "GRCh38":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz"
    elif wildcards.ASSEMBLY == "CHM13v2":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"


def get_transmission_fh(wildcards):
    assembly_adj = wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY
    suffix = SMP2SUFF[wildcards.SAMPLE]

    if wildcards.SAMPLE not in ("2216", "2189", "2209", "2188"):
        return "UNK"
    else:
        return f"{AWS_PREF}/{wildcards.ASSEMBLY}_v1.0_50bp_merge/493ef25/trgt-denovo/transmission/{wildcards.SAMPLE}_{suffix}_{assembly_adj}_50bp_merge_trgtdn.transmission.csv.gz"


def get_bam_for_validation(sample, assembly, tech):

    assembly_adj = assembly.split('v2')[0]
    sample_id = SMP2ALT[sample] if tech in ("ont", "hifi") else sample

    tech2path = {
        "ont": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/ont-bams/{0}/{1}.minimap2.bam".format(assembly_adj, sample_id,), 
        "element": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/element/{0}_bams/{1}-E.{0}.merged.sort.bam".format(assembly_adj, sample_id,),
        "illumina": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH/cram/{0}.cram".format(sample_id),
        "hifi": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in assembly else assembly_adj,),
        }

    return tech2path[tech]

NEW_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025"

def get_complete_hifi_bams(sample_id: str, assembly: str):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """
    if sample_id in ("200100", "2189", "2216", "200080"):
        return NEW_PATH + f"/merged/{sample_id}.merged.bam"
    else:
        assembly_adj = assembly.split('v2')[0]
        alt_sample_id = SMP2ALT[sample_id]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, alt_sample_id, assembly_adj.lower() if 'CHM' in assembly else assembly_adj,)


def get_children_vcfs(wildcards):
    if wildcards.SAMPLE == "2216":
        children = ["200081", "200082", "200084", "200085", "200086", "200087"]
    elif wildcards.SAMPLE == "2189":
        children = ["200101", "200102", "200103", "200104", "200105", "200106"]
    else:
        children = []
    
    return [f"data/trgt/phased/{s}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz" for s in children]

def get_grandparent_vcfs(wildcards):
    if wildcards.SAMPLE.startswith("200"):
        grandparents = ["2209", "2188"]
    else:
        grandparents = ["2281", "2280", "2213", "2214"]
        
    return [f"data/trgt/phased/{s}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz" for s in grandparents]


wildcard_constraints:
    SAMPLE = r"[0-9]{4,6}",
    ASSEMBLY = r"[A-Z]{3}[h]{0,1}[0-9]{0,2}[v]{0,1}[0-9]{1,2}",
    TECH = r"[a-z]+",
    VERSION = "v[0-9].[0-9]"


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))
SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"], ped["suffix"]))

CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

ASSEMBLIES = ["CHM13v2", "GRCh38"]
VERSIONS = ["v4.0"]

rule all:
    input:
        expand("csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, VERSION=VERSIONS,),
        expand("csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, VERSION=VERSIONS, TECH=["element", "hifi"]),
        expand("csv/denominators/{SAMPLE}.{ASSEMBLY}.{VERSION}.denominator.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, VERSION=VERSIONS,),
        expand("csv/recurrent/{ASSEMBLY}.{VERSION}.recurrent.tsv", ASSEMBLY=ASSEMBLIES, VERSION=VERSIONS,),
        # expand("csv/hprc/combined.{ASSEMBLY}.heterozygosity.tsv", ASSEMBLY=["GRCh38"]),
        # expand("csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, VERSION=VERSIONS, TECH=["element", "hifi"])


rule combine_trgt_denovo:
    input:
        fhs = expand("csv/{{SAMPLE}}.{{ASSEMBLY}}.{{VERSION}}.{CHROM}.trgt-denovo.csv", CHROM=CHROMS)
    output: fh = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.trgt-denovo.csv"
    shell:
        """
        
        grep '^chrom' {input.fhs[0]} > {output.fh}

        cat {input.fhs} | grep -v '^chrom'  >> {output.fh}
        """

rule prefilter_all_loci:
    input:
        kid_mutation_df = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.trgt-denovo.csv",
        ped = "tr_validation/data/file_mapping.csv",
        annotations = lambda wildcards: get_annotation_fh(wildcards)
    output: 
        fh = "csv/prefiltered/all_loci/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    params:
        filtering_mode = "lenient"
    script:
        "scripts/filter_mutation_dataframe.py"


rule prefilter_denovos:
    input:
        ped = "tr_validation/data/file_mapping.csv",
        kid_mutation_df = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.trgt-denovo.csv",
        annotations = lambda wildcards: get_annotation_fh(wildcards),
        utils = "scripts/utils.py"
    output: 
        fh = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    params:
        filtering_mode = "strict"
    script:
        "scripts/filter_mutation_dataframe.py"


rule calculate_denominator:
    input:
        loci = "csv/prefiltered/all_loci/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        censat = "data/t2t.censat.bed"
    output:  out = "csv/denominators/{SAMPLE}.{ASSEMBLY}.{VERSION}.denominator.tsv"
    script:
        "scripts/calculate_denominator.py"



rule annotate_with_informative_sites:
    input:
        cohort_snp_vcf = "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.vcf.gz",
        kid_phased_snp_vcf = "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz",
        kid_phased_snp_vcf_idx = "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz.tbi",
        kid_mutation_df = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        kid_phased_str_vcf = "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.vcf.gz",
        py_script = "scripts/annotate_with_informative_sites.py",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.annotated.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
    shell:
        """
        python {input.py_script} --cohort_snp_vcf {input.cohort_snp_vcf} \
                                 --kid_phased_snp_vcf {input.kid_phased_snp_vcf} \
                                 --focal_alt {params.focal_alt_id} \
                                 --focal {wildcards.SAMPLE} \
                                 --out {output.out} \
                                 --mutations {input.kid_mutation_df} \
                                 --str_vcf {input.kid_phased_str_vcf} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id} \
                                 -slop 500000
        """


rule annotate_with_parental_haplotype:
    input:
        py_script = "scripts/annotate_with_parental_haplotype.py",
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.annotated.2gen.tsv",
        dad_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
        dad_phased_snv_vcf = lambda wildcards: f"data/vcf/snv/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
        mom_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
        mom_phased_snv_vcf = lambda wildcards: f"data/vcf/snv/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.haplotyped.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --dad_phased_str_vcf {input.dad_phased_str_vcf} \
                                 --dad_phased_snv_vcf {input.dad_phased_snv_vcf} \
                                 --mom_phased_str_vcf {input.mom_phased_str_vcf} \
                                 --mom_phased_snv_vcf {input.mom_phased_snv_vcf} \
                                 --out {output.out} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id} \
        """


rule phase:
    input:
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.haplotyped.2gen.tsv",
        py_script = "scripts/phase_by_ps.py",
    output: "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """

rule annotate_with_parental_alleles:
    input:
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.phased.2gen.tsv",
        dad_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
        mom_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
        kid_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.VERSION}.phased.vcf.gz",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.allele_sequences.tsv"
    script:
        "scripts/annotate_with_parental_allele_sequences.py"



rule validate_dnms_with_orthogonal_tech:
    input:
        py_script = "scripts/annotate_with_orthogonal_evidence.py",
        kid_mutation_df = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    output: "csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv",
    params:
        kid_bam = lambda wildcards: get_bam_for_validation(wildcards.SAMPLE, wildcards.ASSEMBLY, wildcards.TECH),
        mom_bam = lambda wildcards: get_bam_for_validation(SMP2MOM[wildcards.SAMPLE], wildcards.ASSEMBLY, wildcards.TECH),
        dad_bam = lambda wildcards: get_bam_for_validation(SMP2DAD[wildcards.SAMPLE], wildcards.ASSEMBLY, wildcards.TECH),
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --kid_bam {params.kid_bam} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam}
        """

rule combine_raw_dnms:
    input:
        kid_mutation_dfs = expand("csv/prefiltered/denovos/{SAMPLE}.{{ASSEMBLY}}.{{VERSION}}.tsv", SAMPLE=CHILDREN)
    output:
        fh = "csv/all_dnms/{ASSEMBLY}.{VERSION}.tsv"
    run:
        import pandas as pd
        res = []
        for fh in input.kid_mutation_dfs:
            df = pd.read_csv(fh, sep="\t")
            res.append(df)
        res = pd.concat(res)
        res.to_csv(output.fh, sep="\t", index=False)
        

rule validate_all_with_orthogonal_tech:
    input:
        py_script = "scripts/annotate_with_orthogonal_evidence.py",
        mutation_df = "csv/all_dnms/{ASSEMBLY}.{VERSION}.tsv"
    output: "csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv",
    params:
        kid_bam = lambda wildcards: get_bam_for_validation(wildcards.SAMPLE, wildcards.ASSEMBLY, wildcards.TECH),
    shell:
        """
        python {input.py_script} --mutations {input.mutation_df} \
                                 --kid_bam {params.kid_bam} \
                                 --out {output} \
                                 --tech {wildcards.TECH} 
        """

rule add_orthogonal_filter:
    input:
        mutations = "csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        orthogonal_evidence = "csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv"
    output:
        out = "csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.{VERSION}.tsv"
    script:
        "scripts/add_orthogonal_filter.py"
    

rule add_transmission_evidence:
    input:
        mutations = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        other_vcfs = get_children_vcfs
    output: fh = "csv/transmission/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    params:
        generation_to_query = "children"
    script:
        "scripts/assess_presence_in_cohort.py"


rule add_grandparental_evidence:
    input:
        mutations = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        other_vcfs = get_grandparent_vcfs
    output: fh = "csv/grandparents/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    params:
        generation_to_query = "grandparents"
    script:
        "scripts/assess_presence_in_cohort.py"


# take all of the de novo mutation metadata and merge into one combined file
rule merge_all_dnm_files:
    input:
        raw_denovos = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        phasing = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{VERSION}.allele_sequences.tsv",
        denominator = "csv/denominators/{SAMPLE}.{ASSEMBLY}.{VERSION}.denominator.tsv",
        transmission = "csv/transmission/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv",
        grandparental = "csv/grandparents/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    output: out = "csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{VERSION}.tsv"
    script:
        "scripts/merge_mutations_with_metadata.py"


rule find_recurrents:
    input:
        fhs = expand("csv/filtered_and_merged/{SAMPLE}.{{ASSEMBLY}}.{{VERSION}}.tsv", SAMPLE=CHILDREN)
    output:
        out = "csv/recurrent/{ASSEMBLY}.{VERSION}.recurrent.tsv"
    script:
        "scripts/find_recurrent_trids.py"


rule annotate_with_hprc_heterozygosity:
    input:
        vcf = lambda wildcards: HPRC_PREF + f"{wildcards.ASSEMBLY}_v1.0_50bp_merge/1.1.2-69937d83/hprc/{wildcards.HPRC_SAMPLE}_{wildcards.ASSEMBLY}_50bp_merge.sorted.vcf.gz",
        denovos = "{ASSEMBLY}.filtered.tsv",
        catalog = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY],
    output:
        mutations = "data/hprc/{HPRC_SAMPLE}.{ASSEMBLY}.{KIND}.heterozygosity.tsv"     
    script:
        "scripts/calculate_hprc_heterozygosity.py"


rule merge_hprc:
    input:
        fhs = expand("data/hprc/{HPRC_SAMPLE}.{{ASSEMBLY}}.{KIND}.heterozygosity.tsv", HPRC_SAMPLE=HPRC_SAMPLES, KIND=["wt", "denovo"])
    output: fh = "csv/hprc/combined.{ASSEMBLY}.heterozygosity.tsv"
    shell:
        """
        cat {input.fhs} | grep -v 'chrom' > {output.fh}
        """
