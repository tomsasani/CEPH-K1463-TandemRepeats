import pandas as pd
from typing import Dict
from collections import Counter
import glob


PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

ASSEMBLY = "GRCh38"

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
CHILDREN = [s for s in CHILDREN if s not in ("2188", "2209")]
ALL_SAMPLES = ped["sample_id"].to_list()

HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/1.1.2-69937d83/hprc"

HPRC_SAMPLES = []
for fh in glob.glob(HPRC_PREF + "/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)

cohort2samples = {"hprc": HPRC_SAMPLES, "ceph": ALL_SAMPLES}

COHORT = "hprc"
SAMPLES = cohort2samples[COHORT]

# only look at the trids with multiple motifs
recurrent_dnms = pd.read_csv(f"csv/recurrent/{ASSEMBLY}.v4.0.recurrent.tsv", sep="\t")
recurrent_dnms = recurrent_dnms[recurrent_dnms["motifs"].str.contains(",")]
trids = [t for t in recurrent_dnms["trid"].unique() if not t.startswith("chrX")][:2]

trids = ["chr8_2623352_2623487_trsolve"]

rule all:
    input:
       expand("trviz/{COHORT}/{TRID}.{GENOME}.png", COHORT=[COHORT], TRID=trids, GENOME = [ASSEMBLY]),


rule extract_ceph_fasta:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.v4.0.recurrent.tsv"
    output:
        fasta = temp("data/recurrent_fasta/ceph/{TRID}.{SAMPLE}.{ASSEMBLY}.fa")
    params:
        sample_vcf = lambda wildcards: f"data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.v4.0.phased.vcf.gz",
        sample_dnms = lambda wildcards: f"csv/prefiltered/denovos/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.v4.0.tsv" if wildcards.SAMPLE in CHILDREN else None
    script:
        "scripts/extract_fasta_for_trviz.py"
        

rule extract_hprc_fasta:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.v4.0.recurrent.tsv"
    output:
        fasta = temp("data/recurrent_fasta/hprc/{TRID}.{SAMPLE}.{ASSEMBLY}.fa")
    params:
        sample_vcf = lambda wildcards: f"{HPRC_PREF}/" + wildcards.SAMPLE + f"_{wildcards.ASSEMBLY.lower() if 'CHM' in wildcards.ASSEMBLY else wildcards.ASSEMBLY}" + "_50bp_merge.sorted.vcf.gz",
        sample_dnms = None,
    script:
        "scripts/extract_fasta_for_trviz.py"

rule combine_fasta:
    input:
        fastas = expand("data/recurrent_fasta/{{COHORT}}/{{TRID}}.{SAMPLE}.{{ASSEMBLY}}.fa", SAMPLE=SAMPLES)
    output:
        fasta = "data/recurrent_fasta/{COHORT}/combined/{TRID}.{ASSEMBLY}.combined.fa"
    shell:
        """
        cat {input.fastas} > {output.fasta}
        """

rule decompose:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.v4.0.recurrent.tsv",
        fasta = "data/recurrent_fasta/{COHORT}/combined/{TRID}.{ASSEMBLY}.combined.fa"
    output:
        png = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.png",
        seq_tsv = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.encoded.tsv",
        key_tsv = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.key.tsv",
        aln_fa = "trviz/{COHORT}/{TRID}.{ASSEMBLY}_alignment_output.fa"
    script:
        "scripts/decompose_trviz.py"