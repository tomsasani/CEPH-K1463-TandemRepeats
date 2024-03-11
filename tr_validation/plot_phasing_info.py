import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cyvcf2
from typing import Dict

from build_roc_curve_dnms import get_dp

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

GEN = "3gen"

df = pd.read_csv(
    f"tr_validation/csv/combined.phased.{GEN}.csv",
    dtype={"sample_id": str},
)

if GEN == "2gen":
    df = df[(df["n_consistent_upstream"] >= 1) & (df["n_consistent_downstream"] >= 1)]
else:
    df = df[(df["most_common_freq"] >= 0.8)]

df["mom_dp"] = df["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
df["dad_dp"] = df["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )

ped = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str},
)
smp2vcf = dict(zip(ped["sample_id"], ped["vcf_fh"]))

df = df.query("denovo_coverage >= 2 and \
              child_ratio >= 0.1 and \
              child_coverage >= 10 and \
              mom_dp >= 10 and \
              dad_dp >= 10")

print (df.groupby("sample_id").size())

# df["motif"] = df.apply(lambda row: get_motif(row, smp2vcf), axis=1)
# df["is_homopolymer"] = df["motif"].apply(lambda m: len(m) == 1)

phase_counts = df.groupby(["sample_id", "true_phase"]).size().reset_index().rename(columns={0: "count"})
totals = phase_counts.groupby("sample_id").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})

phase_counts = phase_counts.merge(totals, on="sample_id")

phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]

ages = pd.read_csv("K20 parental age at birth.csv", dtype={"sample_id": str, "UGRP Lab ID (archive)": str})
phase_counts = phase_counts.merge(ages, left_on="sample_id", right_on="UGRP Lab ID (archive)")

phase_counts["annotation"] = phase_counts["fraction"].apply(lambda f: round(f, 2))
phase_counts = phase_counts.astype({"sample_id": "str"})

phase_counts["generation"] = phase_counts["sample_id"].apply(lambda s: "G4" if str(s).startswith("200") else "G2" if s in ("2209", "2188") else "G3")

phase_counts = phase_counts[phase_counts["generation"].isin(["G2", "G3"])]

grp_order = phase_counts.sort_values("total", ascending=True)["sample_id"]
annot = phase_counts[phase_counts["true_phase"] == "dad"].sort_values("total")["annotation"].to_list()
val = phase_counts[phase_counts["true_phase"] == "dad"].sort_values("total")["count"].to_list()

f, axarr = plt.subplots(2, figsize=(10, 8), sharey=True)
for i, (gen, gen_df) in enumerate(phase_counts.groupby("generation")):
    grp_order = gen_df.sort_values("PaAge", ascending=True)["sample_id"]
    sns.barplot(
        data=gen_df,
        x="sample_id",
        y="count",
        hue="true_phase",
        ec="k",
        lw=1,
        order=grp_order,
        ax=axarr[i],
    )
    sns.despine(ax=axarr[i], top=True, right=True)
    if i > 0: axarr[i].legend_.remove()
    axarr[i].set_title(gen)

# for idx, (text, y) in enumerate(zip(annot, val)):
#     ax.annotate(f"{int(text * 100)}%", (idx, y * 1.1))
# ax.set_xlabel("Sample ID")
# ax.set_ylabel("DNM count")
# ax.set_title("DNM phasing in G3 (using HiPhase approach)")
# ax.legend(frameon=False, title="Inferred parent-of-origin")
# sns.despine(ax=ax, top=True, right=True)
f.tight_layout()
f.savefig(f"phase.{GEN}.png", dpi=200)


# new = phase_counts[["sample_id", "true_phase", "MaAge", "PaAge", "count"]].melt(id_vars=["sample_id", "true_phase", "count"])
# phase_counts = phase_counts[phase_counts["sample_id"].isin([2216, 2211, 2212, 2298, 2215, 2217, 2189, 2187])]
# print (phase_counts)

f, ax = plt.subplots()
ax.scatter(phase_counts[phase_counts["true_phase"] == "mom"]["MaAge"], phase_counts[phase_counts["true_phase"] == "mom"]["count"], label="paternal")
ax.scatter(phase_counts[phase_counts["true_phase"] == "dad"]["PaAge"], phase_counts[phase_counts["true_phase"] == "dad"]["count"], label="maternal")
ax.legend()
f.tight_layout()
f.savefig("o.png", dpi=200)
