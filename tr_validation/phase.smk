import pandas as pd

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
#SAMPLES = [s for s in SAMPLES if s not in ("2216", "2209")]
#SAMPLES = ["2216", "2209"]
#SAMPLES = ["2188", "2209", "2189", "2216"]
#SAMPLES = ["2209"]
#SAMPLES = ["2189"]
#SAMPLES = ["200081", ""]

HIPHASE_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/"

rule all:
    input: 
        "tr_validation/csv/phased/combined/phased.2gen.tsv",
        #"tr_validation/csv/michelle.combined.phased.tsv"


rule catalog_inf_sites:
    input:
        py_script = "tr_validation/collect_informative_sites.py",
        phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
    output:
        out = "tr_validation/csv/informative_sites/{SAMPLE}.csv"
    params:
        kid_phase_blocks = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/" + SMP2ALT[wildcards.SAMPLE] + ".GRCh38.hiphase.blocks.tsv",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        kid_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]

    shell:
        """
        python {input.py_script} --phase_blocks {params.kid_phase_blocks} \
                                 --phased_snp_vcf {input.phased_snp_vcf} \
                                 --focal {params.kid_id} \
                                 --out {output.out}\
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """

rule combine_trio_vcfs:
    input: 
    params:
        kid_phased_snp_vcf = lambda wildcards: HIPHASE_PREF + SMP2ALT[wildcards.SAMPLE] + ".GRCh38.deepvariant.glnexus.phased.vcf.gz",
        dad_phased_snp_vcf = lambda wildcards: HIPHASE_PREF + SMP2ALT[SMP2DAD[wildcards.SAMPLE]] + ".GRCh38.deepvariant.glnexus.phased.vcf.gz",
        mom_phased_snp_vcf = lambda wildcards: HIPHASE_PREF + SMP2ALT[SMP2MOM[wildcards.SAMPLE]] + ".GRCh38.deepvariant.glnexus.phased.vcf.gz",
    output:
        "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz"
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
    output: "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz.csi"

    shell:
        """
        bcftools index {input}
        """

rule annotate_with_informative_sites:
    input:
        phased_snp_vcf = "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz",
        #phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz",
        phased_snp_vcf_idx = "tr_validation/data/vcf/{SAMPLE}.trio.vcf.gz.csi",
        py_script = "tr_validation/annotate_with_informative_sites.py",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.annotated.2gen.tsv"
    params:
        kid_phased_str_vcf = lambda wildcards: HIPHASE_PREF + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge.sorted.phased.vcf.gz",
        kid_mutation_df = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.csv.gz",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]

    shell:
        """
        python {input.py_script} --joint_snp_vcf {input.phased_snp_vcf} \
                                 --focal_alt {params.focal_alt_id} \
                                 --focal {wildcards.SAMPLE} \
                                 --out {output.out} \
                                 --mutations {params.kid_mutation_df} \
                                 --str_vcf {params.kid_phased_str_vcf} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """

rule phase:
    input:
        annotated_dnms = "tr_validation/csv/phased/{SAMPLE}.annotated.2gen.tsv",
        py_script = "tr_validation/phase.py",
    output: "tr_validation/csv/phased/{SAMPLE}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """


rule phase_three_gen:
    input:
        phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz",
        ped = "tr_validation/data/file_mapping.csv",
        py_script = "tr_validation/phase_three_gen.py",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.phased.3gen.tsv"
    params:
        kid_mutation_df = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.csv.gz",
        transmission_df = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.transmission.csv.gz",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        kid_id = lambda wildcards: wildcards.SAMPLE
    shell:
        """
        python {input.py_script} --joint_snp_vcf {input.phased_snp_vcf} \
                                 --focal {params.kid_id} \
                                 --out {output.out} \
                                 --mutations {params.kid_mutation_df} \
                                 --ped {input.ped} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id} \
                                 -transmission {params.transmission_df} \
                                 -mutation_type str

        """

rule phase_three_gen_michelle:
    input:
        phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz",
        ped = "tr_validation/data/file_mapping.csv",
        py_script = "tr_validation/phase_three_gen.py",
        kid_mutation_df = "K1463_denovo_snvs.tsv",

    output:
        out = "tr_validation/csv/michelle.{SAMPLE}.phased.3gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        kid_id = lambda wildcards: wildcards.SAMPLE
    shell:
        """
        python {input.py_script} --joint_snp_vcf {input.phased_snp_vcf} \
                                 --focal {params.kid_id} \
                                 --out {output.out} \
                                 --mutations {input.kid_mutation_df} \
                                 --ped {input.ped} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id} \
                                 -mutation_type snv
        """



rule combine_phased:
    input:
        fhs = expand("tr_validation/csv/phased/{SAMPLE}.phased.{{GEN}}.tsv", SAMPLE = SAMPLES)
        
    output:
        fhs = "tr_validation/csv/phased/combined/phased.{GEN}.tsv"
    run:
        import pandas as pd

        res = []
        for fh in input.fhs:
            print (fh)
            df = pd.read_csv(fh, sep="\t")
            res.append(df)
        res_df = pd.concat(res, ignore_index=True)
        res_df.to_csv(output.fhs, index=False, sep="\t")


rule combine_phased_michelle:
    input:
        fhs = expand("tr_validation/csv/michelle.{SAMPLE}.phased.3gen.tsv", zip, SAMPLE = ["2188", "2209", "2189"])
    output:
        fhs = "tr_validation/csv/michelle.combined.phased.tsv"
    run:
        import pandas as pd

        res = []
        for fh in input.fhs:
            df = pd.read_csv(fh, sep="\t")
            res.append(df)
        res_df = pd.concat(res)
        res_df.to_csv(output.fhs, index=False, sep="\t")