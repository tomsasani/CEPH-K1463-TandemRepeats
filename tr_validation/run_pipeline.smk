import pandas as pd
from typing import Dict
import utils

def get_bam_fh(sample: str, bam_pref: str, bam_suff: str, smp2alt: Dict[str, str], use_alt: bool = False):
    sample_id = smp2alt[sample] if use_alt else sample
    return bam_pref + sample_id + bam_suff

def get_children(sample_id: str, ped: pd.DataFrame):
    kid_fhs = []
    parent_id = "maternal_id" if sample_id == "2216" else "paternal_id"
    for kid_id in ped[ped[parent_id] == sample_id]["sample_id"].to_list():
        kid_fhs.append(f"tr_validation/data/insufficient_depth/{kid_id}.tsv")
    return kid_fhs


def get_grandparents(sample_id: str, ped: pd.DataFrame):
    grandparent_fhs = []
    dad, mom = (
        ped[ped["sample_id"] == sample_id]["paternal_id"].values[0],
        ped[ped["sample_id"] == sample_id]["maternal_id"].values[0],
    )
    dad_dad, dad_mom = (
        ped[ped["sample_id"] == dad]["paternal_id"].values[0],
        ped[ped["sample_id"] == dad]["maternal_id"].values[0],
    )
    mom_dad, mom_mom = (
        ped[ped["sample_id"] == mom]["paternal_id"].values[0],
        ped[ped["sample_id"] == mom]["maternal_id"].values[0],
    )
    for sid in (dad_dad, dad_mom, mom_dad, mom_mom):
        grandparent_fhs.append(f"tr_validation/data/insufficient_depth/{sid}.tsv")

    return grandparent_fhs

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

ASSEMBLY = "GRCh38"
# ASSEMBLY = "CHM13v2"
BAM_TECH = "element"

TECH2PREF = {"ont": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/ont-bams/{ASSEMBLY.split('v2')[0]}/", 
             "element": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/element/{ASSEMBLY.split('v2')[0]}_bams/",
             "illumina": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/CEPH/cram/",
             "hifi": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/hifi-bams/{ASSEMBLY.split('v2')[0]}/"}
TECH2SUFF = {"ont": ".minimap2.bam", 
             "element": f"-E.{ASSEMBLY.split('v2')[0]}.merged.sort.bam",
             "illumina": ".cram",
             "hifi": f".{ASSEMBLY.split('v2')[0].lower() if 'CHM' in ASSEMBLY else ASSEMBLY.split('v2')[0]}.haplotagged.bam"}

ILLUMINA_REF = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

ASSEMBLY2CATLOG = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"}

bam_pref, bam_suff = TECH2PREF[BAM_TECH], TECH2SUFF[BAM_TECH]

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

AWS_PREF = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25/"

rule all:
    input:
        expand("tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/csv/rates/{SAMPLE}.{ASSEMBLY}.denominator.tsv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.3gen.tsv", SAMPLE=["2216", "2189"], ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.tsv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY], TECH=[BAM_TECH]),
        # expand("tr_validation/data/all_loci/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY], TECH=["illumina", "element", "ont", "hifi"]),
        # expand("tr_validation/data/all_loci/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY], TECH=["element"]),

rule prefilter_all_loci:
    input:
        py_script = "tr_validation/utils.py",
    output: 
        fh = "tr_validation/data/all_loci/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    params:
        kid_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY]
    run:
        import pandas as pd
        annotations = pd.read_csv(params.annotations, sep="\t")
        # merge annotations with kid mutation dataframe
        mutations = pd.read_csv(params.kid_mutation_df, sep="\t")
        mutations = utils.filter_mutation_dataframe(mutations, 
                                                    denovo_coverage_min=0,
                                                    remove_duplicates=True,
                                                    remove_complex=True,
                                                    depth_min=10, 
                                                    parental_overlap_frac_max=1)
        # filter to just the homozygous sites
        # mutations["is_homozygous"] = mutations["child_AL"].apply(lambda a: a.split(",")[0] == a.split(",")[1])
        # mutations = mutations[mutations["is_homozygous"] == True]
        mutations["sample_id"] = wildcards.SAMPLE
        merged_mutations = mutations.merge(annotations, left_on="trid", right_on="TRid")
        # filter to homopolymers
        # merged_mutations = merged_mutations[merged_mutations["motifs"].str.len() == 1].sample(100)
        merged_mutations[["sample_id", "trid", "index", "genotype", "motifs", "n_motifs", "child_AL"]].to_csv(output.fh, sep="\t", index=False)


rule validate_all_loci:
    input:
        kid_mutation_df="tr_validation/data/all_loci/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        py_script="tr_validation/validate_with_bams.py",
    output: "tr_validation/data/all_loci/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
    params:
        vcf = lambda wildcards: AWS_PREF + "hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.phased.vcf.gz",
        kid_bam=lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam=lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2MOM[wildcards.SAMPLE] != "missing" else None,
        dad_bam=lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2DAD[wildcards.SAMPLE] != "missing" else None,
        max_al=lambda wildcards: 10_000 if wildcards.TECH in ("ont", "hifi") else 120,
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --vcf {params.vcf} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --variant_type dnm \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                                 -max_al {params.max_al} \                                 
        """


rule prefilter_denovos:
    input:
        "tr_validation/utils.py"
    output: 
        fh = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    params:
        kid_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
    run:
        import pandas as pd
        mutations = pd.read_csv(params.kid_mutation_df, sep="\t")
        print ("####", wildcards.ASSEMBLY, mutations[mutations["child_MC"].str.contains("_")].shape)

        # basic filtering, only remove stuff that's not a putative de novo
        mutations = utils.filter_mutation_dataframe(mutations,
                                                    remove_complex = False,
                                                    remove_duplicates = False,
                                                    remove_gp_ev = False,
                                                    remove_inherited = True,
                                                    parental_overlap_frac_max = 1,
                                                    denovo_coverage_min = 2,
                                                    depth_min = 10,)
        print ("####", wildcards.ASSEMBLY, mutations[mutations["child_MC"].str.contains("_")].shape)
        mutations["sample_id"] = wildcards.SAMPLE
        mutations.to_csv(output.fh, sep="\t", index=False)


# merge raw de novo calls with grandparental evidence
rule merge_with_grandparental_evidence:
    input: 
        py_script = "tr_validation/annotate_with_grandparental_evidence.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        ped = "tr_validation/data/file_mapping.csv"
    output:
        out = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.grandparental.tsv"
    params: 
        vcf_pref = AWS_PREF,
        vcf_suff = lambda wildcards: f"{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.vcf.gz"
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --ped {input.ped} \
                                 --focal {wildcards.SAMPLE} \
                                 --out {output.out} \
                                 --vcf_pref {params.vcf_pref} \
                                 --vcf_suff {params.vcf_suff}
        """

# merge raw de novo calls with transmission evidence
rule annotate_with_transmission:
    input:
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        py_script = "tr_validation/annotate_with_transmission.py",
    output: "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.transmission.tsv"
    params:
        kid_transmission_df = lambda wildcards: AWS_PREF + "trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.transmission." + f"{'v1' if 'CHM' in wildcards.ASSEMBLY else 'v1.T1'}" + ".csv.gz" if wildcards.SAMPLE in ("2216", "2189") else "UNK",
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --transmission {params.kid_transmission_df} \
                                 --out {output}
        """

# take all of the de novo mutation metadata and merge into one combined file
rule merge_all_dnm_files:
    input:
        raw_denovos = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        grandparental_evidence = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.grandparental.tsv",
        transmission_evidence = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.transmission.tsv",
        phasing = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.2gen.tsv",
        py_script = "tr_validation/merge_mutations_with_metadata.py",
    output: "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv"
    shell:
        """
        python {input.py_script} --raw_denovos {input.raw_denovos} \
                                 --grandparental_evidence {input.grandparental_evidence} \
                                 --transmission_evidence {input.transmission_evidence} \
                                 --phasing {input.phasing} \
                                 --out {output} \
        """


# query orthogonal BAMs for evidence at each of the raw denovo loci
rule validate_with_orthogonal_tech:
    input:
        py_script="tr_validation/validate_with_bams.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    output: "tr_validation/data/denovos/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
    params:
        vcf = lambda wildcards: AWS_PREF + "hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.phased.vcf.gz",
        kid_bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, bam_pref, bam_suff, SMP2ALT, use_alt=True if BAM_TECH != "element" else False),
        mom_bam = lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], bam_pref, bam_suff, SMP2ALT, use_alt=True if BAM_TECH != "element" else False) if SMP2MOM[wildcards.SAMPLE] != "missing" else None,
        dad_bam = lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], bam_pref, bam_suff, SMP2ALT, use_alt=True if BAM_TECH != "element" else False) if SMP2DAD[wildcards.SAMPLE] != "missing" else None,
        max_al = lambda wildcards: 10_000 if wildcards.TECH != "element" else 120,
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --vcf {params.vcf} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --variant_type dnm \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                                 -max_al {params.max_al}                            
        """

# merge raw de novo calls with orthogonal read evidence
rule merge_with_orthogonal_evidence:
    input:
        py_script = "tr_validation/annotate_with_orthogonal_evidence.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        kid_orthogonal = "tr_validation/data/denovos/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv"
    output: "tr_validation/csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.tsv"
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --orthogonal_evidence {input.kid_orthogonal} \
                                 --out {output}
        """


rule calculate_grouped_denominator:
    input:
        py_script = "tr_validation/calculate_grouped_denominator.py",
        insufficient_depth_sites = expand("tr_validation/data/insufficient_depth/{KID_ID}.tsv", KID_ID=ped["sample_id"].unique()),
    output: "tr_validation/csv/rates/{SAMPLE}.{ASSEMBLY}.denominator.tsv"
    params:
        children_sites = lambda wildcards: get_children(wildcards.SAMPLE, ped) if wildcards.SAMPLE in ("2189", "2216") else None,
        grandparent_sites = lambda wildcards: get_grandparents(wildcards.SAMPLE, ped),
        kid_unfiltered_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY],
    shell:
        """
        python {input.py_script} --denominator {params.kid_unfiltered_mutation_df} \
                                 --annotations {params.annotations} \
                                 --out {output} \
                                 --sample_id {wildcards.SAMPLE} \
                                -children_sites {params.children_sites} \
                                -grandparent_sites {params.grandparent_sites}

        """




rule find_trids_with_insufficient_depth:
    input:
    output: fh = "tr_validation/data/insufficient_depth/{SAMPLE}.tsv"
    params:
        kid_vcf = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge.sorted.phased.vcf.gz",
    run:
        from cyvcf2 import VCF
        import tqdm
        import numpy as np
        from collections import Counter

        res = [] 
        vcfh = VCF(params.kid_vcf, gts012=True)
        for i, v in tqdm.tqdm(enumerate(vcfh)):
            gt = v.gt_types[0]
            trid = v.INFO.get("TRID")
            if gt == 3:
                res.append(trid)
            else:
                depth = np.sum(v.format("SD"))
                if depth < 10: 
                    res.append(trid)
        res_counts = Counter(res).most_common()
        res_df = pd.DataFrame(res_counts, columns=["trid", "count"])
        res_df["sample_id"] = wildcards.SAMPLE
        res_df.to_csv(output.fh, sep="\t", index=False)


rule combine_trio_vcfs:
    input: 
    params:
        kid_phased_snp_vcf = lambda wildcards: AWS_PREF + "hiphase/" + SMP2ALT[wildcards.SAMPLE] + f".{ASSEMBLY.lower() if 'CHM' in ASSEMBLY else ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz",
        dad_phased_snp_vcf = lambda wildcards: AWS_PREF + "hiphase/" + SMP2ALT[SMP2DAD[wildcards.SAMPLE]] + f".{ASSEMBLY.lower() if 'CHM' in ASSEMBLY else ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz",
        mom_phased_snp_vcf = lambda wildcards: AWS_PREF + "hiphase/" + SMP2ALT[SMP2MOM[wildcards.SAMPLE]] + f".{ASSEMBLY.lower() if 'CHM' in ASSEMBLY else ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz",
    output:
        temp("tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz")
    threads: 4
    shell:
        """
        bcftools merge {params.kid_phased_snp_vcf} \
                       {params.dad_phased_snp_vcf} \
                       {params.mom_phased_snp_vcf} \
                       --threads 4 \
                       | \
        bcftools view -m2 -M2 \
                       --include 'INFO/AN == 6' \
                       -o {output} \
                       -Oz \
        """


rule index_trio_vcf:
    input: "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz"
    output: temp("tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz.csi")

    shell:
        """
        bcftools index {input}
        """


rule annotate_with_informative_sites:
    input:
        phased_snp_vcf = "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz",
        phased_snp_vcf_idx = "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz.csi",
        py_script = "tr_validation/annotate_with_informative_sites.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.annotated.2gen.tsv"
    params:
        kid_phased_str_vcf = lambda wildcards: AWS_PREF + "hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{ASSEMBLY.lower() if 'CHM' in ASSEMBLY else ASSEMBLY}_50bp_merge.sorted.phased.vcf.gz",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]

    shell:
        """
        python {input.py_script} --joint_snp_vcf {input.phased_snp_vcf} \
                                 --focal_alt {params.focal_alt_id} \
                                 --focal {wildcards.SAMPLE} \
                                 --out {output.out} \
                                 --mutations {input.kid_mutation_df} \
                                 --str_vcf {params.kid_phased_str_vcf} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """


rule phase:
    input:
        annotated_dnms = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.annotated.2gen.tsv",
        py_script = "tr_validation/phase.py",
    output: "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """


rule phase_three_gen:
    input:
        ped = "tr_validation/data/file_mapping.csv",
        py_script = "tr_validation/phase_three_gen.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.3gen.tsv"
    params:
        phased_snp_vcf = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/" + "CEPH-1463.joint." + f"{'chm13' if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + ".deepvariant.glnexus.phased.vcf.gz",
        transmission_df = lambda wildcards: AWS_PREF + "trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.transmission.csv.gz" if wildcards.SAMPLE in ("2216", "2189") else "UNK",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        kid_id = lambda wildcards: wildcards.SAMPLE
    shell:
        """
        python {input.py_script} --joint_snp_vcf {params.phased_snp_vcf} \
                                 --focal {params.kid_id} \
                                 --out {output.out} \
                                 --mutations {input.kid_mutation_df} \
                                 --ped {input.ped} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id} \
                                 -transmission {params.transmission_df} \
                                 -mutation_type str

        """