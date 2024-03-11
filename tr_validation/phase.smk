import pandas as pd

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))


SMP2SUFF = dict(zip(ped["sample_id"].to_list(), [f.split("/")[-1].split("_")[1] for f in ped["vcf_fh"].to_list()]))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
#SAMPLES = ["2188", "2209", "2189", "2216"]
#SAMPLES = ["2209"]
#SAMPLES = ["2189"]

rule all:
    input: 
        "tr_validation/csv/combined.phased.2gen.csv",


rule catalog_inf_sites:
    input:
        py_script = "tr_validation/collect_informative_sites.py",
        phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
    output:
        out = "tr_validation/csv/informative_sites/{SAMPLE}.csv"
    params:
        kid_phase_blocks = lambda wildcards: "tr_validation/data/hiphase/" + SMP2ALT[wildcards.SAMPLE] + ".GRCh38.hiphase.blocks.tsv",
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


rule phase_two_gen:
    input:
        informative_sites = "tr_validation/csv/informative_sites/{SAMPLE}.csv",
        py_script = "tr_validation/phase.py",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.phased.2gen.csv"
    params:
        kid_phased_str_vcf = lambda wildcards: "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge.sorted.phased.vcf.gz",
        kid_phase_blocks = lambda wildcards: "tr_validation/data/hiphase/" + SMP2ALT[wildcards.SAMPLE] + ".GRCh38.hiphase.blocks.tsv",
        kid_mutation_df = lambda wildcards: "tr_validation/data/denovos/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.csv.gz",
        dad_id = lambda wildcards: SMP2DAD[wildcards.SAMPLE],
        mom_id = lambda wildcards: SMP2MOM[wildcards.SAMPLE]
    shell:
        """
        python {input.py_script} --informative_sites {input.informative_sites} \
                                 --focal {wildcards.SAMPLE} \
                                 --phase_blocks {params.kid_phase_blocks} \
                                 --out {output.out} \
                                 --mutations {params.kid_mutation_df} \
                                 --str_vcf {params.kid_phased_str_vcf} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """

rule phase_three_gen:
    input:
        phased_snp_vcf = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz",
        ped = "tr_validation/data/file_mapping.csv",
        py_script = "tr_validation/phase_three_gen.py",
    output:
        out = "tr_validation/csv/phased/{SAMPLE}.phased.3gen.csv"
    params:
        kid_mutation_df = lambda wildcards: "tr_validation/data/denovos/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.csv.gz",
        transmission_df = lambda wildcards: "tr_validation/data/transmission/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + "_GRCh38_50bp_merge_trgtdn.transmission.csv.gz",
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        kid_id = lambda wildcards: wildcards.SAMPLE
    shell:
        """
        python {input.py_script} --joint_snp_vcf {input.phased_snp_vcf} \
                                 --focal {params.kid_id} \
                                 --out {output.out} \
                                 --mutations {params.kid_mutation_df} \
                                 --transmission {params.transmission_df} \
                                 --ped {input.ped} \
                                 --dad_id {params.dad_id} \
                                 --mom_id {params.mom_id}
        """

rule combine_phased:
    input:
        fhs = expand("tr_validation/csv/phased/{SAMPLE}.phased.{{GEN}}.csv", zip, SAMPLE = SAMPLES)
    output:
        fhs = "tr_validation/csv/combined.phased.{GEN}.csv"
    run:
        import pandas as pd

        res = []
        for fh in input.fhs:
            df = pd.read_csv(fh)
            res.append(df)
        res_df = pd.concat(res)
        res_df.to_csv(output.fhs, index=False)