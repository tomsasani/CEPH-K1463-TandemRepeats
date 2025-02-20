import pandas as pd 
from typing import Dict
import utils


CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX"]

ASSEMBLY = "GRCh38"

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str,})

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
PARENTS = ["2209", "2188"]
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

ASSEMBLY2CATLOG = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz"}
ASSEMBLY2DEFINITIONS = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz",
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


def get_bam_fh(sample: str, bam_pref: str, bam_suff: str, smp2alt: Dict[str, str], use_alt: bool = False, member: str = "kid"):
    if member == "mom":
        sample = SMP2MOM[sample]
    elif member == "dad":
        sample = SMP2DAD[sample]
    sample_id = smp2alt[sample] if use_alt else sample
    return bam_pref + sample_id + bam_suff

# CHILDREN = ["200106"]

CHILDREN = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()[:2]
CHILDREN = ped[ped["paternal_id"].isin(["200080", "2189"])]["sample_id"].to_list()

CHILDREN = ["200087", "200106", "2187"]
print (CHILDREN)

PREF = "tr_validation/downsampled"
AWS_PREF = f"/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25/"


KID_TARGETS, MOM_TARGETS, DAD_TARGETS = [], [], []
# kid is always "downsampled" to about the same as they originally were
# mom and dad are downsampled such that one parent is constant at 50X and the other varies
for mom_target in (10, 20, 30, 40):
    KID_TARGETS.append(50)
    DAD_TARGETS.append(50)
    MOM_TARGETS.append(mom_target)
for dad_target in (10, 20, 30, 40):
    KID_TARGETS.append(50)
    MOM_TARGETS.append(50)
    DAD_TARGETS.append(dad_target)

KID_TARGETS, MOM_TARGETS, DAD_TARGETS = (list(range(10, 60, 10)), list(range(10, 70, 10)), list(range(10, 70, 10)),)

KID_TARGETS = [50, 50, 50, 50, 50]
DAD_TARGETS = [70, 80, 90, 70, 70]
MOM_TARGETS = [70, 70, 70, 80, 90]

print (list(zip(KID_TARGETS, MOM_TARGETS, DAD_TARGETS)))

rule all:
    input:
        # expand("tr_validation/downsampled/csv/2187.GRCh38.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.prefiltered.tsv", zip, KID_TARGET=KID_TARGETS, MOM_TARGET=MOM_TARGETS, DAD_TARGET=DAD_TARGETS),
        expand("tr_validation/downsampled/csv/2187.GRCh38.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.phased.2gen.tsv", zip, KID_TARGET=KID_TARGETS, MOM_TARGET=MOM_TARGETS, DAD_TARGET=DAD_TARGETS)
    # expand(PREF + "/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY]),
    # expand(PREF + "/csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.tsv", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY], TECH=["element"])
    # "tr_validation/mosdepth.html",
    # expand(PREF + "/data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.{MEMBER}.vcf.gz", SAMPLE=CHILDREN, ASSEMBLY=[ASSEMBLY], MEMBER=["dad"]),
    # expand("tr_validation/data/mosdepth/{SAMPLE}.mosdepth.global.dist.txt", SAMPLE=ALL_SAMPLES),

# rule intersect_originals:
#     input:
#         original = "tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv",
#         new = PREF + "/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv",
#     output:
#         out = PREF + "/csv/filtered_and_merged/intersected/{SAMPLE}.{ASSEMBLY}.tsv"
#     run:
#         import pandas

#         columns = ["trid", "sample_id", "phase_summary", "denovo_coverage", "child_coverage", "dad_dp", "mom_dp", "denovo_al"]

#         original = pd.read_csv(input.original, sep="\t")[columns]
#         new = pd.read_csv(input.new, sep="\t")[columns]

#         intersected = original.merge(new, how="cross")
#         intersected = intersected[intersected["trid_x"] == intersected["trid_y"]]

#         intersected.to_csv(output.out, index=False, sep="\t")



rule run_mosdepth:
    input:
        mosdepth = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/mosdepth"
    output:
        dist = "tr_validation/data/mosdepth/{SAMPLE}.mosdepth.global.dist.txt",
        by_region = "tr_validation/data/mosdepth/{SAMPLE}.regions.bed.gz"
    params:
        bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF["hifi"], TECH2SUFF["hifi"], SMP2ALT, use_alt=True)
    shell:
        """
        {input.mosdepth} --by 1000 --fast-mode -n tr_validation/data/mosdepth/{wildcards.SAMPLE} {params.bam}
        """
        
rule plot_mosdepth:
    input:
        fhs = expand("tr_validation/data/mosdepth/{SAMPLE}.mosdepth.global.dist.txt", SAMPLE=ALL_SAMPLES),
        py_script = "tr_validation/data/mosdepth/plot.py"
    output:
        "tr_validation/mosdepth.html"
    shell:
        """
        python {input.py_script} -o {output} {input.fhs}
        """


rule create_definition_file:
    input:
        mutations = "tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    output:
        mean_depths = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
        definitions = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.definitions.txt",
        bed = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.bed"
    params:
        kid_depth = lambda wildcards: f"tr_validation/data/mosdepth/{wildcards.SAMPLE}.mosdepth.summary.txt",
        mom_depth = lambda wildcards: f"tr_validation/data/mosdepth/{SMP2MOM[wildcards.SAMPLE]}.mosdepth.summary.txt",
        dad_depth = lambda wildcards: f"tr_validation/data/mosdepth/{SMP2DAD[wildcards.SAMPLE]}.mosdepth.summary.txt"
    run:
        import pandas as pd

        mutations = pd.read_csv(input.mutations, sep="\t")
        mutations["start_adj"] = mutations["start"].values - 50_000
        mutations["start_adj"] = mutations["start_adj"].apply(lambda s: 1 if s < 0 else s)
        mutations["end_adj"] = mutations["end"].values + 50_000
        # create a BED file and defintiions file
        mutations[["#chrom", "start_adj", "end_adj"]].to_csv(output.bed, header=False, index=False, sep="\t")

        with open(output.definitions, "w") as outfh:
            for i, row in mutations.iterrows():
                chrom, start, end = row["#chrom"], row["start"], row["end"]
                trid = row["trid"]
                struc, motifs = row["struc"], row["motifs"]
                out_formatted = f"{chrom}\t{start}\t{end}\tID={trid};MOTIFS={motifs};STRUC={struc}"
                print (out_formatted, file=outfh)
        outfh.close()

        kid_depth = pd.read_csv(params.kid_depth, sep="\t")
        mean_kid_depth = kid_depth[kid_depth["chrom"] == "total"]["mean"].values[0]

        dad_depth = pd.read_csv(params.dad_depth, sep="\t")
        mean_dad_depth = dad_depth[dad_depth["chrom"] == "total"]["mean"].values[0]

        mom_depth = pd.read_csv(params.mom_depth, sep="\t")
        mean_mom_depth = mom_depth[mom_depth["chrom"] == "total"]["mean"].values[0]
        
        with open(output.mean_depths, "w") as outfh:
            print ("\t".join(["dad_dp_mean", "mom_dp_mean", "kid_dp_mean"]), file=outfh)
            print ("\t".join(list(map(str, [mean_dad_depth, mean_mom_depth, mean_kid_depth]))), file=outfh)
        outfh.close()


rule create_annotations_bed:
    input:
        definitions = ASSEMBLY2CATLOG[ASSEMBLY]
    output:
        definitions = PREF + "/data/definitions.bed"
    shell:
        """
        gzip -d -c {input} > {output}
        """


rule downsample_bam:
    input:
        py_script = "tr_validation/create_downsampling_cmd.py",
        mean_depths = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
        subsample_regions = PREF + "/data/definitions.bed",
    output: 
        bam = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
    params:
        bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF["hifi"], TECH2SUFF["hifi"], SMP2ALT, use_alt=True, member=wildcards.MEMBER),
    threads: 16
    shell:
        """
        module load sambamba

        echo {wildcards.MEMBER}

        CMD=$(python {input.py_script} --depth {input.mean_depths} \
                                 --bed {input.subsample_regions} \
                                 --input_bam {params.bam} \
                                 --output_bam {output.bam} \
                                 -target_depth {wildcards.TARGET} \
                                 -nthreads 16 \
                                 -member {wildcards.MEMBER})

        echo $CMD
        echo $CMD | sh
        """


rule index_bam:
    input:
        bam = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam"
    output:
        bam_idx = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai"
    shell:
        """
        module load samtools
        
        samtools index {input.bam}
        """


rule run_trgt:
    input:
        bam = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        bam_idx = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai",
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta",
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt",
        repeats = ASSEMBLY2CATLOG[ASSEMBLY]
    output:
        vcf = temp(PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.vcf.gz"),
        bam = temp(PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.spanning.bam")
    threads: 8
    shell:
        """
        {input.trgt_binary} genotype --threads 8 \
                            --genome {input.reference} \
                            --repeats {input.repeats} \
                            --reads {input.bam} \
                            --output-prefix tr_validation/downsampled/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.MEMBER}.{wildcards.TARGET}
        """



rule sort_bam:
    input:
        PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.spanning.bam"
    output:
        PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.spanning.sorted.bam"    
    shell:
        """
        module load samtools
        
        samtools sort -o {output} {input}
        samtools index {output}
        """

rule sort_vcf:
    input:
        PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.vcf.gz"
    output:
        vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.sorted.vcf.gz",
        idx = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.sorted.vcf.gz.csi"
    shell:
        """
        module load bcftools
        
        bcftools sort -Ob -o {output.vcf} {input}
        bcftools index {output.vcf}
        """


rule get_orig_snv_bed:
    input:
    output:
        PREF + "/data/snv_regions/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{CHROM}.bed"
    params:
        snp_vcf = lambda wildcards: f"{AWS_PREF}hiphase/{SMP2ALT[wildcards.SAMPLE] if wildcards.MEMBER == 'kid' else SMP2ALT[SMP2MOM[wildcards.SAMPLE]] if wildcards.SAMPLE == 'mom' else SMP2ALT[SMP2MOM[wildcards.SAMPLE]]}.{ASSEMBLY.lower() if 'CHM' in ASSEMBLY else ASSEMBLY}.deepvariant.glnexus.phased.vcf.gz",
    shell:
        """
        echo {params.snp_vcf}
        bcftools view -H {params.snp_vcf} {wildcards.CHROM} | cut -f 1-2 | awk -v OFS='\t' '{{print $1, $2-1, $2}}' > {output}
        """


rule call_snvs:
    input:
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta",
        bam = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        bam_idx = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai",
        regions = PREF + "/data/snv_regions/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{CHROM}.bed"
    output:
        vcf = PREF + "/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.bcf",
    shell:
        """
        bcftools mpileup -r {wildcards.CHROM} -Ou -f {input.reference} --annotate FORMAT/DP {input.bam} | bcftools call -V indels -mv -f GQ -Ob -o {output.vcf}
        """


rule combine_snvs:
    input:
        vcfs = expand(PREF + "/data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{MEMBER}}.{{TARGET}}.{CHROM}.bcf", CHROM=CHROMS)
    output:
        PREF + "/data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.vcf.gz"
    shell:
        """
        bcftools concat -Oz -o {output} {input.vcfs}
        bcftools index {output}
        """


rule run_hiphase:
    input:
        snv_vcf = PREF + "/data/vcf/snv/raw/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.vcf.gz",
        str_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.sorted.vcf.gz",
        str_vcf_idx = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.sorted.vcf.gz.csi",
        bam = PREF + "/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"
    output:
        snv_vcf = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz",
        str_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz",
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
        reference = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta",
        repeats = PREF + "/data/definitions.bed",#ASSEMBLY2CATLOG[ASSEMBLY], # PREF + "/definitions/{SAMPLE}.{ASSEMBLY}.tsv",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.sorted.vcf.gz",
        mom_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.sorted.vcf.gz",
        dad_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.sorted.vcf.gz",
        kid_bam = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.spanning.sorted.bam",
        mom_bam = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.spanning.sorted.bam",
        dad_bam = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.spanning.sorted.bam",
    output:
        output = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv"
    threads: 8
    params:
        kid_pref = lambda wildcards: PREF + f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.kid.{wildcards.KID_TARGET}",
        mom_pref = lambda wildcards: PREF + f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}",
        dad_pref = lambda wildcards: PREF + f"/data/trgt/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}",
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


rule prefilter_denovos:
    input:
        py_script = "tr_validation/utils.py",
        ped = "tr_validation/data/file_mapping.csv",
        kid_mutation_df = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv"
    output: 
        fh = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.prefiltered.tsv"
    params:
        annotations = lambda wildcards: ASSEMBLY2DEFINITIONS[wildcards.ASSEMBLY],
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
        kid_mutation_df = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        py_script = "tr_validation/annotate_with_transmission.py",
    output: PREF + "/csv/{SAMPLE}.{ASSEMBLY}.transmission.tsv"
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
        raw_denovos = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv",
        transmission_evidence = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.transmission.tsv",
        phasing = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.phased.2gen.tsv",
        py_script = "tr_validation/merge_mutations_with_metadata.py",
    output: PREF + "/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv"
    shell:
        """
        python {input.py_script} --raw_denovos {input.raw_denovos} \
                                 --transmission_evidence {input.transmission_evidence} \
                                 --phasing {input.phasing} \
                                 --out {output} \
        """

rule combine_trio_vcfs:
    input: 
        kid_phased_snv_vcf = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz",
        mom_phased_snv_vcf = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.phased.vcf.gz",
        dad_phased_snv_vcf = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.phased.vcf.gz",
    output:
        PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.trio.vcf.gz"
    threads: 4
    shell:
        """
        bcftools merge {input.kid_phased_snv_vcf} \
                       {input.dad_phased_snv_vcf} \
                       {input.mom_phased_snv_vcf} \
                       --threads 4 \
                       | \
        bcftools view -m2 -M2 \
                       --include 'INFO/AN >= 4' \
                       -o {output} \
                       -Oz \
        """


rule index_trio_vcf:
    input: PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.trio.vcf.gz"
    output: PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.trio.vcf.gz.csi"
    shell:
        """
        bcftools index {input}
        """


rule annotate_with_informative_sites:
    input:
        phased_snp_vcf = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.trio.vcf.gz",
        phased_snp_vcf_idx = PREF + "/data/vcf/snv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.trio.vcf.gz.csi",
        py_script = "tr_validation/annotate_with_informative_sites.py",
        kid_mutation_df = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.prefiltered.tsv",
        # for the kid's phased STR VCF, it only matters what they were donsampled to
        phased_str_vcf = PREF + "/data/trgt/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz",
    output:
        out = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.annotated.2gen.tsv"
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
                                 --str_vcf {input.phased_str_vcf} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """


rule phase:
    input:
        annotated_dnms = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.annotated.2gen.tsv",
        py_script = "tr_validation/phase.py",
    output: PREF + "/csv/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """

rule validate_with_orthogonal_tech:
    input:
        py_script="tr_validation/validate_with_bams.py",
        kid_mutation_df = PREF + "/csv/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv"
    output: PREF + "/csv/read_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
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
        kid_mutation_df = PREF + "/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv",
        kid_orthogonal = PREF + "/csv/read_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv"
    output: PREF + "/csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.tsv"
    shell:
        """
        python {input.py_script} --mutations {input.kid_mutation_df} \
                                 --orthogonal_evidence {input.kid_orthogonal} \
                                 --out {output}
        """