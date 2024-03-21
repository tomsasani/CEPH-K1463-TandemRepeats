import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cyvcf2
from typing import Dict, List
from utils import filter_mutation_dataframe
import glob
from build_roc_curve_dnms import get_dp, annotate_with_allr, check_for_overlap

def get_motif(row: pd.Series, smp2vcf: Dict[str, str]):
    trid = row["trid"]
    chrom = trid.split("_")[0]
    start, end = trid.split("_")[1:3]
    region = f"{chrom}:{start}-{end}"
    sample_id = row["sample_id"]
    vcf = cyvcf2.VCF(smp2vcf[sample_id])
    motif = None

    for v in vcf(region):
        v_trid = v.INFO.get("TRID")
        if v_trid != trid: continue
        motif = v.INFO.get("MOTIFS")
    return motif

GEN = "2gen"

# loop over phased mutation dataframes
df = pd.read_csv(
    "tr_validation/csv/phased/combined/phased.2gen.tsv",
    sep="\t",
    dtype={"sample_id": str},
)

trid_counts = df.groupby("trid").size().to_dict()
duplicate_trids = [k for k,v in trid_counts.items() if v > 1]
duplicates = df[df["trid"].isin(duplicate_trids)]

print (duplicates[duplicates["sample_id"] == "200082"])

df = filter_mutation_dataframe(
        df,
        remove_complex=False,
        remove_duplicates=True,
        mc_column="child_MC"
    )


MIN_INF_SITES = 1

if GEN == "2gen":
    df = df[(df["n_upstream"] >= MIN_INF_SITES) | (df["n_downstream"] >= MIN_INF_SITES)]
else:
    df = df[(df["most_common_freq"] >= 0.8) & (df["candidate_postzygotic"] == "N")]

df["mom_dp"] = df["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
df["dad_dp"] = df["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )
df["maternal_dp_ratio"] = df["mom_dp"] / df["dad_dp"]

df["father_overlap_coverage"] = df["father_overlap_coverage"].apply(lambda f: get_dp(f))
df["mother_overlap_coverage"] = df["mother_overlap_coverage"].apply(lambda f: get_dp(f))

df["dad_min_support"] = df["per_allele_reads_father"].apply(lambda a: min(list(map(int, a.split(",")))))
df["mom_min_support"] = df["per_allele_reads_mother"].apply(lambda a: min(list(map(int, a.split(",")))))

ped = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str},
)

MIN_ALLELE_COVERAGE = 1
MIN_DENOVO_COVERAGE = 2
MIN_ALLELE_RATIO = 0.
MIN_CHILD_RATIO = 0.
MAX_CHILD_RATIO = 1
MIN_PARENT_SUPPORT = 10

df = df.query(f"allele_coverage >= {MIN_ALLELE_COVERAGE} and \
              denovo_coverage >= {MIN_DENOVO_COVERAGE} and \
              allele_ratio >= {MIN_ALLELE_RATIO} and \
              child_ratio >= {MIN_CHILD_RATIO} and \
              child_ratio <= {MAX_CHILD_RATIO} and \
              father_dropout == 'N' and \
              mother_dropout == 'N' and \
              dad_min_support >= {MIN_PARENT_SUPPORT} and \
              mom_min_support >= {MIN_PARENT_SUPPORT}")

print (df.groupby("sample_id").size())

df["true_phase"] = df["phase_summary"].apply(lambda p: p.split(":")[0])
df["support"] = df["phase_summary"].apply(lambda p: int(p.split(":")[1]))

print (df[df["sample_id"] == "200101"].head()[["trid", "genotype", "phase_summary", "non_denovo_al_diff", "denovo_al_diff"]])

df["generation"] = df["sample_id"].apply(lambda s: "G4" if str(s).startswith("200") else "G2" if s in ("2209", "2188") else "G3")

g = sns.FacetGrid(data=df, row="generation", sharex=False, sharey=False, aspect=3)
g.map(sns.boxplot, "sample_id", "support", "true_phase", dodge=True)
g.add_legend()
g.savefig("support.png")

COLS = ["sample_id", "true_phase"]

phase_counts = df.groupby(COLS).size().reset_index().rename(columns={0: "count"})
totals = phase_counts.groupby("sample_id").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})

phase_counts = phase_counts.merge(totals)
phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]

ages = pd.read_csv("K20 parental age at birth.csv", dtype={"sample_id": str, "UGRP Lab ID (archive)": str})
phase_counts["generation"] = phase_counts["sample_id"].apply(lambda s: "G4" if str(s).startswith("200") else "G2" if s in ("2209", "2188") else "G3")

f, axarr = plt.subplots(len(phase_counts["generation"].unique()), figsize=(10, 8), sharey=False)

g = sns.FacetGrid(data=phase_counts, row="generation", sharex=False, sharey=False, aspect=3)
g.map(sns.barplot,
        "sample_id",
        "count",
        "true_phase",
    )

g.figure.suptitle("Counts of phased DNMs (using HiPhase approach)")
g.add_legend()
g.tight_layout()
g.savefig(f"phase.{GEN}.png", dpi=200)
