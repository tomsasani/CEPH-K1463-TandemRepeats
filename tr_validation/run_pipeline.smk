import pandas as pd
from typing import Dict
import utils
from collections import Counter

def get_bam_fh(sample: str, bam_pref: str, bam_suff: str, smp2alt: Dict[str, str], use_alt: bool = False):
    sample_id = smp2alt[sample] if use_alt else sample
    return bam_pref + sample_id + bam_suff

def get_children(sample_id: str, ped: pd.DataFrame):
    kid_fhs = []
    parent_id = "maternal_id" if sample_id == "2216" else "paternal_id"
    for kid_id in ped[ped[parent_id] == sample_id]["sample_id"].to_list():
        kid_fhs.append(f"tr_validation/data/insufficient_depth/{kid_id}.tsv")
    return kid_fhs

def assign_pedigree_id(row: pd.Series):
    comb_id = []
    if row["sample_id"] in ("2280", "2281") or row["paternal_id"] == "2281":
        comb_id.append("G2A")
    if row["sample_id"] in ("2214", "2213") or row["paternal_id"] == "2214":
        comb_id.append("G2B")
    if row["sample_id"] in ("2209", "2188") or row["paternal_id"] == "2209":
        comb_id.append("G3")
    if row["sample_id"] in ("2216", "200080") or row["paternal_id"] == "200080":
        comb_id.append("G4A")
    if row["sample_id"] in ("2189", "200100") or row["paternal_id"] == "2189":
        comb_id.append("G4B")

    return ",".join(comb_id)


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

ASSEMBLY = "CHM13v2"
# ASSEMBLY = "GRCh38"

TECH2PREF = {"ont": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/ont-bams/{ASSEMBLY.split('v2')[0]}/", 
             "element": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/element/{ASSEMBLY.split('v2')[0]}_bams/",
             "illumina": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/CEPH/cram/",
             "hifi": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/hifi-bams/{ASSEMBLY.split('v2')[0]}/"}
TECH2SUFF = {"ont": ".minimap2.bam", 
             "element": f"-E.{ASSEMBLY.split('v2')[0]}.merged.sort.bam",
             "illumina": ".cram",
             "hifi": f".{ASSEMBLY.split('v2')[0].lower() if 'CHM' in ASSEMBLY else ASSEMBLY.split('v2')[0]}.haplotagged.bam",}

ILLUMINA_REF = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

ASSEMBLY2CATLOG = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"}


POLYMORPHIC_TRIDS = "/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/T2T/50bp_merge/493ef25/xy_added/polymorphic_per_AL_TRIDs.txt"

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
PARENTS = ["2209", "2188"]
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

AWS_PREF = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25/"

rule all:
    input:
        expand("tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.tsv", SAMPLE=SAMPLES, ASSEMBLY=[ASSEMBLY], TECH=["hifi", "element"]),
        expand("tr_validation/data/denovos/orthogonal_support_recurrent/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY], TECH=["hifi"]),
        expand("tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.3gen.tsv", SAMPLE=["2189", "2216"], ASSEMBLY=[ASSEMBLY]),

rule prefilter_all_loci:
    input:
        py_script = "tr_validation/utils.py",
        polymorphic_trids = POLYMORPHIC_TRIDS#"tr_validation/data/polymorphic.txt"

    output: 
        fh = "tr_validation/data/all_loci/{SAMPLE}.{ASSEMBLY}.{METHOD}.prefiltered.tsv"
    params:
        kid_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY]
    run:
        import pandas as pd
        annotations = pd.read_csv(params.annotations, sep="\t")
        # merge annotations with kid mutation dataframe
        mutations = pd.read_csv(params.kid_mutation_df, sep="\t", dtype={"child_AL": str, 
                                                                         "mother_AL": str, 
                                                                         "father_AL": str, 
                                                                         "per_allele_reads_mother": str, 
                                                                         "per_allele_reads_father": str, 
                                                                         "per_allele_reads_child": str,
                                                                         "father_overlap_coverage": str,
                                                                         "mother_overlap_coverage": str})
        mutations = utils.filter_mutation_dataframe(mutations, 
                                                    denovo_coverage_min=0,
                                                    remove_duplicates=True,
                                                    remove_complex=False,
                                                    depth_min=10, 
                                                    parental_overlap_frac_max=1,)

        mutations = mutations[mutations["denovo_coverage"] == 0]

        # if desired, filter to just a sample of homozygous sites
        if wildcards.METHOD == "sample":
            mutations["is_homozygous"] = mutations["child_AL"].apply(lambda a: (len(a.split(",")) == 1) or (a.split(",")[0] == a.split(",")[1]))
            mutations = mutations[mutations["is_homozygous"] == True]
    
        mutations["sample_id"] = wildcards.SAMPLE
        merged_mutations = mutations.merge(annotations, left_on="trid", right_on="TRid")

        # filter to homopolymers if desired
        if wildcards.METHOD == "sample":
            merged_mutations = merged_mutations[merged_mutations["motifs"].str.len() == 1].sample(1_000)
        
        # if we're validating everything, we'll stick to STRs
        elif wildcards.METHOD == "all":

            merged_mutations["motif_size"] = merged_mutations.apply(lambda row: utils.determine_motif_size(row), axis=1)
            merged_mutations["simple_motif_size"] = merged_mutations.apply(lambda row: utils.determine_simplified_motif_size(row), axis=1)

            merged_mutations = merged_mutations[merged_mutations["simple_motif_size"] == "STR"]

        merged_mutations[["sample_id", "trid", "index", "genotype", "motifs", "n_motifs", "child_AL"]].to_csv(output.fh, sep="\t", index=False)


rule validate_all_loci:
    input:
        kid_mutation_df="tr_validation/data/all_loci/{SAMPLE}.{ASSEMBLY}.{METHOD}.prefiltered.tsv",
        py_script="tr_validation/validate_with_bams.py",
    output: "tr_validation/data/all_loci/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.{METHOD}.read_support.csv",
    params:
        kid_bam=lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam=lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2MOM[wildcards.SAMPLE] != "missing" else None,
        dad_bam=lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2DAD[wildcards.SAMPLE] != "missing" else None,
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam}
        """



rule prefilter_denovos:
    input:
        py_script = "tr_validation/utils.py",
        ped = "tr_validation/data/file_mapping.csv",
        polymorphic_trids = POLYMORPHIC_TRIDS,
    output: 
        fh = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    params:
        kid_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY],
    run:
        import pandas as pd
        mutations = pd.read_csv(params.kid_mutation_df, sep="\t", dtype={"child_AL": str, 
                                                                         "mother_AL": str, 
                                                                         "father_AL": str, 
                                                                         "per_allele_reads_mother": str, 
                                                                         "per_allele_reads_father": str, 
                                                                         "per_allele_reads_child": str,
                                                                         "father_overlap_coverage": str,
                                                                         "mother_overlap_coverage": str})
        mutations["sample_id"] = wildcards.SAMPLE

        # trids_to_use = []
        # with open(input.polymorphic_trids, "r") as infh:
        #     for l in infh:
        #         trids_to_use.append(l.strip())
        # mutations = mutations[mutations["trid"].isin(trids_to_use)]

        file_mapping = pd.read_csv(input.ped, dtype={"sample_id": str})
        mutations = mutations.merge(file_mapping)

        mutations = utils.filter_mutation_dataframe(mutations,
                                                    remove_complex = False,
                                                    remove_duplicates = False,
                                                    remove_gp_ev = False,
                                                    remove_inherited = True,
                                                    parental_overlap_frac_max = 0.05,
                                                    denovo_coverage_min = 2,
                                                    depth_min = 10,
                                                    child_ratio_min = 0.2)

        annotations = pd.read_csv(params.annotations, sep="\t")
        merged_mutations = mutations.merge(annotations, left_on="trid", right_on="TRid")

        merged_mutations["motif_size"] = merged_mutations.apply(lambda row: utils.determine_motif_size(row), axis=1)
        merged_mutations["simple_motif_size"] = merged_mutations.apply(lambda row: utils.determine_simplified_motif_size(row), axis=1)
        merged_mutations["sample_id"] = wildcards.SAMPLE
        merged_mutations.to_csv(output.fh, sep="\t", index=False)


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
        vcf_suff = lambda wildcards: f"{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.vcf.gz",
        lineage = lambda wildcards: "maternal" if wildcards.SAMPLE.startswith("20008") else "paternal" if wildcards.SAMPLE.startswith("2001") else "both"
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --ped {input.ped} \
                                 --focal {wildcards.SAMPLE} \
                                 --out {output.out} \
                                 --vcf_pref {params.vcf_pref} \
                                 --vcf_suff {params.vcf_suff} \
                                 --lineage {params.lineage}
        """

# merge raw de novo calls with transmission evidence
rule annotate_with_transmission:
    input:
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        py_script = "tr_validation/annotate_with_transmission.py",
    output: "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.transmission.tsv"
    params:
        kid_transmission_df = lambda wildcards: AWS_PREF + "trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.transmission." + f"{'T1' if 'CHM' in wildcards.ASSEMBLY else 'v1.T1'}" + ".csv.gz" if wildcards.SAMPLE in ("2216", "2189") else "UNK",
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
        transmission_evidence = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.transmission.tsv",
        phasing = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.phased.2gen.tsv",
        py_script = "tr_validation/merge_mutations_with_metadata.py",
        recurrents = "tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv",
        denominator = "tr_validation/csv/rates/{SAMPLE}.{ASSEMBLY}.denominator.tsv",
    output: "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv"
    shell:
        """
        python {input.py_script} --raw_denovos {input.raw_denovos} \
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
        # vcf = lambda wildcards: AWS_PREF + "hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.phased.vcf.gz",
        kid_bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam = lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2MOM[wildcards.SAMPLE] != "missing" else None,
        dad_bam = lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2DAD[wildcards.SAMPLE] != "missing" else None,
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam}
        """

# merge raw de novo calls with orthogonal read evidence
rule merge_with_orthogonal_evidence:
    input:
        py_script = "tr_validation/annotate_with_orthogonal_evidence.py",
        kid_mutation_df = "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv",
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
        utils = "tr_validation/utils.py",
        polymorphic_trids = POLYMORPHIC_TRIDS,
    output: "tr_validation/csv/rates/{SAMPLE}.{ASSEMBLY}.denominator.tsv"
    params:
        kid_unfiltered_mutation_df = lambda wildcards: AWS_PREF + "trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.csv.gz",
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY],

    shell:
        """
        python {input.py_script} --denominator {params.kid_unfiltered_mutation_df} \
                                 --annotations {params.annotations} \
                                 --out {output} \
                                 --sample_id {wildcards.SAMPLE} \
                                 -children_sites None \
                                 -trids None
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
        "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz"
    threads: 4
    shell:
        """
        bcftools merge {params.kid_phased_snp_vcf} \
                       {params.dad_phased_snp_vcf} \
                       {params.mom_phased_snp_vcf} \
                       --threads 4 \
                       | \
        bcftools view -m2 -M2 \
                       --include 'INFO/AN >= 4' \
                       -o {output} \
                       -Oz \
        """


rule index_trio_vcf:
    input: "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz"
    output: "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz.csi"

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
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
        phased_snp_vcf = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/" + "CEPH-1463.joint." + f"{'chm13' if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + ".deepvariant.glnexus.phased.vcf.gz",


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
        transmission_df = lambda wildcards: AWS_PREF + "trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.transmission.csv.gz" if wildcards.SAMPLE in ("2216", "2189", "2209", "2188") else "UNK",
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
                                 --transmission {params.transmission_df} \
                                 -mutation_type str

        """

rule get_recurrent_trids:
    input:
        fhs = expand("tr_validation/data/denovos/{SAMPLE}.{{ASSEMBLY}}.prefiltered.tsv", SAMPLE=CHILDREN),
    output:
        fh = "tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv"
    run:
        import pandas as pd

        dfs = []
        for fh in input.fhs:
            df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
            dfs.append(df)
        dfs = pd.concat(dfs)

        samples_with_denovo = dfs.groupby("trid").agg(samples_with_denovo=("sample_id", lambda s: ",".join(s))).reset_index()
        dfs = dfs.merge(samples_with_denovo)

        dfs["generation"] = dfs["sample_id"].apply(lambda s: "G4" if s.startswith("200") else "G2" if s in ("2209", "2188") else "G3")

        # figure out which TRIDs are observed as DNs
        # multiple times in a single generation, rather than multiple
        # times across generations
        generational = dfs.groupby("trid").agg(inter=("generation", lambda g: len(set(g))), intra=("generation", lambda g: sum([gi == "G3" for gi in g])),).reset_index()

        intra_trids = generational[generational["intra"] > 1]["trid"].unique()
        inter_trids = generational[generational["inter"] > 1]["trid"].unique()
        combined_trids = set(intra_trids).union(set(inter_trids))

        print ("INTRA-GENERATIONAL", len(intra_trids))
        print ("INTER-GENERATIONAL", len(inter_trids))
        # output a subsetted dataframe with the information necessary to do orthogonal validation
        dfs = dfs[dfs["trid"].isin(combined_trids)].drop_duplicates("trid")
        dfs = dfs.merge(generational)

        dfs[["trid", "#chrom", "start", "end", "index", "motifs", "child_AL", "samples_with_denovo"]].drop_duplicates("trid").to_csv(output.fh, index=False, sep="\t")


# query orthogonal BAMs for evidence at each of the raw denovo loci
rule validate_recurrent_with_orthogonal_tech:
    input:
        py_script="tr_validation/validate_with_bams.py",
        trids = "tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv"
    output: "tr_validation/data/denovos/orthogonal_support_recurrent/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
    params:
        # vcf = lambda wildcards: AWS_PREF + "hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.phased.vcf.gz",
        kid_bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam = None,
        dad_bam = None,
    shell:
        """
        python {input.py_script} --mutations {input.trids} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                        
        """