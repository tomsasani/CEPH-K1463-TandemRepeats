import pandas as pd 
from typing import Dict
import utils

ASSEMBLY = "GRCh38"

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
PARENTS = ["2209", "2188"]
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

ASSEMBLY2CATLOG = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"}



TECH2PREF = {"ont": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/ont-bams/{ASSEMBLY.split('v2')[0]}/", 
             "element": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/element/{ASSEMBLY.split('v2')[0]}_bams/",
             "illumina": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/CEPH/cram/",
             "hifi": f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/hifi-bams/{ASSEMBLY.split('v2')[0]}/"}

TECH2SUFF = {"ont": ".minimap2.bam", 
             "element": f"-E.{ASSEMBLY.split('v2')[0]}.merged.sort.bam",
             "illumina": ".cram",
             "hifi": f".{ASSEMBLY.split('v2')[0].lower() if 'CHM' in ASSEMBLY else ASSEMBLY.split('v2')[0]}.haplotagged.bam",}

def snv_fh(wildcards):

    sample_id = None
    if wildcards.MEMBER == "kid":
        sample_id = SMP2ALT[wildcards.SAMPLE] 
    elif wildcards.MEMBER == 'mom':
        sample_id = SMP2ALT[SMP2MOM[wildcards.SAMPLE]]
    elif wildcards.MEMBER == 'dad':
        sample_id = SMP2ALT[SMP2DAD[wildcards.SAMPLE]]
    
    pref = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{wildcards.ASSEMBLY}_v1.0_50bp_merge/493ef25/hiphase/"
    
    return pref + f"{sample_id}.{wildcards.ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz"


def get_bam_fh(sample: str, bam_pref: str, bam_suff: str, smp2alt: Dict[str, str], use_alt: bool = False):
    sample_id = smp2alt[sample] if use_alt else sample
    return bam_pref + sample_id + bam_suff

CHILDREN = ["200106"]

rule all:
    input:
        # expand("tr_validation/definitions/{SAMPLE}.{ASSEMBLY}.tsv", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY]),
        # expand("tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{TECH}.{MEMBER}.vcf.gz", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY], TECH=["hifi"], MEMBER=["mom", "dad", "kid"])
        expand("tr_validation/trgt_denovo_out/{SAMPLE}.{ASSEMBLY}.denovo.csv", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.subsampled.phased.vcf.gz", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY]),
        expand("tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.phased.vcf.gz", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY], MEMBER=["kid", "mom", "dad"])

rule create_definition_file:
    input:
        mutations = "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv",
    output:
        definitions = "tr_validation/definitions/{SAMPLE}.{ASSEMBLY}.tsv",
        mean_depths = "tr_validation/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
        subsample_regions = "tr_validation/csv/{SAMPLE}.{ASSEMBLY}.bed"
    run:
        import pandas as pd

        mutations = pd.read_csv(input.mutations, sep="\t")

        mean_dad_depth = mutations["dad_dp"].mean()
        mean_mom_depth = mutations["mom_dp"].mean()

        with open(output.mean_depths, "w") as outfh:
            print ("\t".join(["dad_dp_mean", "mom_dp_mean", "maternal_ratio"]), file=outfh)
            print ("\t".join(list(map(str, [mean_dad_depth, mean_mom_depth, mean_mom_depth / mean_dad_depth]))), file=outfh)
        outfh.close()

        with open(output.subsample_regions, "w") as outfh:
            for i, row in mutations.iterrows():
                chrom, start, end = row["#chrom"], int(row["start"]), int(row["end"])
                print ("\t".join(["#chrom", "start", "end"]), file=outfh)
                print ("\t".join(list(map(str, [chrom, start - 50_000, end + 50_000]))), file=outfh)
            outfh.close()

        with open(output.definitions, "w") as outfh:
            for i, row in mutations.iterrows():
                chrom, start, end = row["#chrom"], str(row["start"]), str(row["end"])
                ident = row["trid"]
                motifs = row["motifs"]
                struc = row["struc"]

                info = f"ID={ident};MOTIFS={motifs};STRUC={struc}"
                out = "\t".join([chrom, start, end, info])
                print(out, file=outfh)

        outfh.close()


rule subsample_snp_vcf:
    input:
        subsample_regions = "tr_validation/csv/{SAMPLE}.{ASSEMBLY}.bed"
    output:
        "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.vcf.gz"
    params:
        snv_vcf = lambda wildcards: snv_fh(wildcards)
    shell:
        """
        module load bcftools
        
        bcftools view -R {input.subsample_regions} -Oz -o {output} {params.snv_vcf}
        bcftools index {output}
        """


rule downsample_bam:
    input:
        py_script = "tr_validation/create_downsampling_cmd.py",
        mean_depths = "tr_validation/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
        subsample_regions = "tr_validation/csv/{SAMPLE}.{ASSEMBLY}.bed"
    output: 
        dad_bam = temp("tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.dad.subsampled.bam"),
        mom_bam = temp("tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.mom.subsampled.bam"),
        kid_bam = temp("tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.kid.subsampled.bam"),
    params:
        kid_bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF["hifi"], TECH2SUFF["hifi"], SMP2ALT, use_alt=True),
        mom_bam = lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], TECH2PREF["hifi"], TECH2SUFF["hifi"], SMP2ALT, use_alt=True),
        dad_bam = lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], TECH2PREF["hifi"], TECH2SUFF["hifi"], SMP2ALT, use_alt=True),

    threads: 8
    shell:
        """
        module load samtools

        CMD=$(python {input.py_script} --depth {input.mean_depths} \
                                 --bed {input.subsample_regions} \
                                 --dad_input_bam {params.dad_bam} \
                                 --mom_input_bam {params.mom_bam} \
                                 --kid_input_bam {params.kid_bam} \
                                 --dad_output_bam {output.dad_bam} \
                                 --mom_output_bam {output.mom_bam} \
                                 --kid_output_bam {output.kid_bam})

        echo $CMD
        echo $CMD | sh
        """

rule index_bam:
    input:
        bam = "tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.bam"
    output:
        bam_idx = "tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.bam.bai"
    shell:
        """
        module load samtools
        
        samtools index {input.bam}
        """


rule run_trgt:
    input:
        bam = "tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.bam",
        bam_idx = "tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.bam.bai",
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta",
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt",
        repeats = "tr_validation/definitions/{SAMPLE}.{ASSEMBLY}.tsv"
    output:
        vcf = temp("tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz"),
        bam = temp("tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.spanning.bam")
    shell:
        """
        {input.trgt_binary} genotype --threads 1 \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            --output-prefix tr_validation/data/trgt_out/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.MEMBER}
        """

rule index_trgt_vcf:
    input:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz"
    output:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz.csi"
    shell:
        """
        module load bcftools

        bcftools index {input}
        """



rule run_hiphase:
    input:
        snv_vcf = "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.subsampled.vcf.gz",
        str_vcf = "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz",
        str_vcf_idx = "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz.csi",
        bam = "tr_validation/data/subsampled/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.bam",
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"
    output:
        snv_vcf = "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.{MEMBER}.subsampled.phased.vcf.gz",
        str_vcf = "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.phased.vcf.gz",
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


rule sort_bam:
    input:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.spanning.bam"
    output:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.spanning.sorted.bam"
    shell:
        """
        module load samtools
        
        samtools sort -o {output} {input}
        samtools index {output}
        """

rule sort_vcf:
    input:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz"
    output:
        "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.{MEMBER}.sorted.vcf.gz"
    shell:
        """
        module load bcftools
        
        bcftools sort -Ob -o {output} {input}
        bcftools index {output}
        """


rule run_trgt_denovo:
    input:
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta",
        repeats = "tr_validation/definitions/{SAMPLE}.{ASSEMBLY}.tsv",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo-v0.1.2-linux_x86_64",
        vcfs = expand("tr_validation/data/trgt_out/{{SAMPLE}}.{{ASSEMBLY}}.{MEMBER}.sorted.vcf.gz", MEMBER=["kid", "mom", "dad"]),
        bams = expand("tr_validation/data/trgt_out/{{SAMPLE}}.{{ASSEMBLY}}.{MEMBER}.spanning.sorted.bam", MEMBER=["kid", "mom", "dad"])
    output:
        output = "tr_validation/trgt_denovo_out/{SAMPLE}.{ASSEMBLY}.denovo.csv"
    params:
        kid_pref = lambda wildcards: f"tr_validation/data/trgt_out/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.kid",
        mom_pref = lambda wildcards: f"tr_validation/data/trgt_out/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom",
        dad_pref = lambda wildcards: f"tr_validation/data/trgt_out/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad",
    shell:
        """
        {input.trgt_denovo_binary} trio --reference {input.reference} \
                                        --bed {input.repeats} \
                                        --father {params.dad_pref} \
                                        --mother {params.mom_pref} \
                                        --child {params.kid_pref} \
                                        --out {output}
        """


rule prefilter_denovos:
    input:
        py_script = "tr_validation/utils.py",
        ped = "tr_validation/data/file_mapping.csv",
        kid_mutation_df = "tr_validation/trgt_denovo_out/{SAMPLE}.{ASSEMBLY}.hifi.denovo.annotated.csv"
    output: 
        fh = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    params:
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY],
    run:
        import pandas as pd
        mutations = pd.read_csv(input.kid_mutation_df, sep="\t", dtype={"child_AL": str, 
                                                                         "mother_AL": str, 
                                                                         "father_AL": str, 
                                                                         "per_allele_reads_mother": str, 
                                                                         "per_allele_reads_father": str, 
                                                                         "per_allele_reads_child": str,
                                                                         "father_overlap_coverage": str,
                                                                         "mother_overlap_coverage": str})
        mutations["sample_id"] = wildcards.SAMPLE

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




# merge raw de novo calls with transmission evidence
rule annotate_with_transmission:
    input:
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        py_script = "tr_validation/annotate_with_transmission.py",
    output: "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.transmission.tsv"
    params:
        kid_transmission_df = lambda wildcards: AWS_PREF + "trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.ASSEMBLY}" + "_50bp_merge_trgtdn.transmission." + "v1.T1" + ".csv.gz" if wildcards.SAMPLE in ("2216", "2189") else "UNK",
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
    output: "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv"
    shell:
        """
        python {input.py_script} --raw_denovos {input.raw_denovos} \
                                 --transmission_evidence {input.transmission_evidence} \
                                 --phasing {input.phasing} \
                                 --out {output} \
        """

rule combine_trio_vcfs:
    input: 
        kid_phased_snv_vcf = "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.kid.subsampled.phased.vcf.gz",
        mom_phased_snv_vcf = "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.mom.subsampled.phased.vcf.gz",
        dad_phased_snv_vcf = "tr_validation/data/snv_vcf/{SAMPLE}.{ASSEMBLY}.dad.subsampled.phased.vcf.gz",
    output:
        "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz"
    threads: 4
    shell:
        """
        bcftools merge {params.kid_phased_snv_vcf} \
                       {params.dad_phased_snv_vcf} \
                       {params.mom_phased_snv_vcf} \
                       --threads 4 \
                       | \
        bcftools view -m2 -M2 \
                       --include 'INFO/AN >= 4' \
                       -o {output} \
                       -Oz \
        """


rule index_trio_vcf:
    input: "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz"
    output: "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz.csi"

    shell:
        """
        bcftools index {input}
        """


rule annotate_with_informative_sites:
    input:
        phased_snp_vcf = "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz",
        phased_snp_vcf_idx = "tr_validation/data/vcf/{SAMPLE}.{ASSEMBLY}.trio.vcf.gz.csi",
        py_script = "tr_validation/annotate_with_informative_sites.py",
        kid_mutation_df = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        phase_str_vcf = "tr_validation/data/trgt_out/{SAMPLE}.{ASSEMBLY}.kid.phased.vcf.gz",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.{ASSEMBLY}.annotated.2gen.tsv"
    params:
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
                                 --str_vcf {params.phased_str_vcf} \
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
