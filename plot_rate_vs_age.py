import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss

def calculate_poisson_ci(ct: int, obs: int, alpha: float = 0.95):

    l, u = (1 - alpha) / 2, 1 - ((1 - alpha) / 2)

    lower_bound = ss.chi2.ppf(l, ct * 2) / 2
    upper_bound = ss.chi2.ppf(u, (ct + 1) * 2) / 2

    # lower_bound = ss.poisson.ppf(l, ct)
    # upper_bound = ss.poisson.ppf(u, ct)
    
    # Convert bounds to rates per observation
    rate_lower = lower_bound / obs
    rate_upper = upper_bound / obs
    
    return (rate_lower, rate_upper)

pd.set_option("display.precision", 8)
plt.rc("font", size=16)

# define some global variables we'll access for filtering
FILTER_TRANSMITTED = False
PER_HAPLOTYPE = True
FILTER_RECURRENT = True
FILTER_ELEMENT = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "GRCh38"
# ASSEMBLY = "CHM13v2"

ORTHOGONAL_TECH = "element"

# define the minimum size of events we want to consider
# (i.e., if MIN_SIZE = 1, don't include interruptions)
MIN_SIZE = 1

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == "2209"]
# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

# if we want to calculate mutation rates per-haplotype, rather
# than per-genome, we need to multiply the denominator by 2 to account
# for there being 2 copies of every locus in a diploid genome.
if PER_HAPLOTYPE:
    denoms["denominator"] *= 2

mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else p)

metadata = pd.read_csv("tr_validation/data/k20_parental_age_at_birth.csv", dtype={"UGRP Lab ID (archive)": str})

mutations = mutations.merge(metadata, left_on="sample_id", right_on="UGRP Lab ID (archive)")
mutations = mutations[mutations["phase"] != "unknown"]
mutations["TR_type"] = mutations["motif_size"].apply(lambda m: "homopolymer" if m == 1 else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR" if m > 6 else "?")

sample_counts = mutations.groupby(["sample_id", "PaAge", "MaAge", "phase", "TR_type"]).size().reset_index().rename(columns={0: "count"})
sample_counts["age"] = sample_counts.apply(lambda row: row["PaAge"] if row["phase"] == "dad" else row["MaAge"], axis=1)

# compute total denominators in each sample
denoms["TR_type"] = denoms["motif_size"].apply(
    lambda m: (
        "homopolymer"
        if m == 1
        else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR" if m > 6 else "?"
    )
)
per_sample_denoms = (
    denoms.groupby(["sample_id", "TR_type"]).agg({"denominator": "sum"}).reset_index()
)

sample_counts = sample_counts.merge(per_sample_denoms)

sample_counts[
    [
        "sample_id",
        "age",
        "phase",
        "TR_type",
        "count",
        "denominator",
    ]
].to_csv("age_effect.tsv", sep="\t", index=False)

sample_totals = sample_counts.groupby(["sample_id", "PaAge", "MaAge", "TR_type"]).agg(total=("count", "sum")).reset_index()
sample_counts = sample_counts.merge(sample_totals)

sample_counts = sample_counts[sample_counts["TR_type"] != "VNTR"]

sample_alphas = (
    sample_counts[sample_counts["phase"] == "dad"]
    .groupby(["sample_id", "PaAge"])
    .agg(count=("count", "sum"), total=("total", "sum"))
    .reset_index()
)
sample_alphas["alpha"] = sample_alphas["count"] / sample_alphas["total"]

sample_alphas.to_csv("alphas.tsv", sep="\t", index=False)

f, ax = plt.subplots(figsize=(8, 6))
sns.regplot(data=sample_alphas, x="PaAge", y="alpha", ax=ax, scatter_kws={"edgecolor": "w", "s": 150, "linewidths": 2})
ax.set_xlabel("Paternal age")
ax.set_ylabel(r"$\alpha$")# + "\n(fraction of DNMs phased to dad)")
ax.set_title(r"$\alpha$" + " for TR DNMs\nincreases with paternal age")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("alpha.png", dpi=200)

sample_counts.rename(columns={"phase": "Parent-of-origin"}, inplace=True)
sample_counts["rate"] = sample_counts["count"] / sample_counts["denominator"]

g = sns.lmplot(
    data=sample_counts,
    col="Parent-of-origin",
    x="age",
    y="count",
    hue="TR_type",
    palette="colorblind",
    scatter_kws={"edgecolor": "w", "s": 150, "linewidths": 2},
)
g.set_axis_labels(x_var="Parental age", y_var="Number of TR DNMs")
g.tight_layout()
g.savefig("age_effects.png", dpi=200)
