import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss
import statsmodels.api as sm
import statsmodels.formula.api as smf


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

def get_motif_types(row):
    if row["n_motifs"] > 1:
        if row["max_motiflen"] > 6 and row["min_motiflen"] <= 6:
            return "complex"
        else:
            if row["min_motiflen"] == 1:
                return "homopolymer"
            elif row["min_motiflen"] > 1 and row["max_motiflen"] <= 6:
                return "non-homopolymer STR"
            else:
                return "VNTR"
    else:
        if row["min_motiflen"] == 1:
            return "homopolymer"
        elif 1 < row["min_motiflen"] <= 6:
            return "non-homopolymer STR"
        else:
            return "VNTR"

pd.set_option("display.precision", 8)
plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

ASSEMBLY = "CHM13v2"
PER_HAPLOTYPE = False

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", dtype={"sample_id": str, "paternal_id": str}, sep="\t")

mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4A" if s.startswith("2000") else "G4B" if s.startswith("2001") else "G3")

# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"csv/denominators/*.{ASSEMBLY}.v4.0.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

# if we want to calculate mutation rates per-haplotype, rather
# than per-genome, we need to multiply the denominator by 2 to account
# for there being 2 copies of every locus in a diploid genome.
if PER_HAPLOTYPE:
    denoms["denominator"] *= 2


metadata = pd.read_csv("tr_validation/data/k20_parental_age_at_birth.csv", dtype={"UGRP Lab ID (archive)": str})

mutations = mutations.merge(metadata, left_on="sample_id", right_on="UGRP Lab ID (archive)")

mutations = mutations[mutations["phase"] != "unknown"]

mutations["TR type"] = mutations.apply(lambda row: get_motif_types(row), axis=1)

sample_counts = (
    mutations.groupby(
        [
            "sample_id",
            "PaAge",
            "phase",
            "TR type",
            "generation",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

# compute total denominators in each sample
denoms["TR type"] = denoms["motif_size"].apply(
    lambda m: (
        "homopolymer"
        if m == 1
        else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR" if m > 6 else "?"
    )
)
per_sample_denoms = (
    denoms.groupby(["sample_id"]).agg({"denominator": "sum"}).reset_index()
)

sample_counts = sample_counts.merge(per_sample_denoms)
sample_counts["rate"] = sample_counts["count"] / sample_counts["denominator"]

sample_totals = sample_counts.groupby(["sample_id"]).agg(total=("count", "sum")).reset_index()
sample_counts = sample_counts.merge(sample_totals)

g = sns.lmplot(data=sample_counts,  x="PaAge", y="rate", row="generation", col="TR type", hue="phase")
g.savefig("test.png")


alpha_counts = (
    mutations.groupby(
        [
            "sample_id",
            "PaAge",
            "phase",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
alpha_totals = alpha_counts.groupby(["sample_id", "PaAge"]).agg(total=("count", "sum")).reset_index()
alpha_counts = alpha_counts.merge(alpha_totals)
alpha_counts = alpha_counts[alpha_counts["phase"] == "dad"]
alpha_counts["alpha"] = alpha_counts["count"] / alpha_counts["total"]


model = sm.OLS(alpha_counts["alpha"], alpha_counts["PaAge"]).fit()
print (model.summary())

g = sns.lmplot(
    data=alpha_counts,
    x="PaAge",
    y="alpha",
    scatter_kws={
        "edgecolor": "w",
        "s": 150,
        "linewidths": 2,
    },
)

g.savefig("alpha.png", dpi=200)

sample_counts = (
    mutations.groupby(
        [
            "sample_id",
            "PaAge",
            "phase",
            "TR type",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
sample_totals = sample_counts.groupby(["sample_id"]).agg(total=("count", "sum")).reset_index()
sample_counts = sample_counts.merge(sample_totals)
print (sample_counts)

model = smf.ols(
    formula="count ~ PaAge",
    data=sample_counts[(sample_counts["phase"] == "dad") & (sample_counts["TR type"] == "non-homopolymer STR")],
).fit()

print (model.summary())

model = smf.ols(
    formula="count ~ PaAge",
    data=sample_counts[(sample_counts["phase"] == "dad") & (sample_counts["TR type"] == "complex")],
).fit()

print (model.summary())


sample_counts.rename(columns={"phase": "Parent-of-origin"}, inplace=True)

g = sns.lmplot(
    data=sample_counts[sample_counts["Parent-of-origin"] != "unknown"],
    x="PaAge",
    y="count",
    hue="Parent-of-origin",
    col="TR type",
    aspect=1.1,
    col_wrap=2,
    # palette="deep",
    #scatter_kws={"edgecolor": "w", "s": 150, "linewidths": 2},
)
g.set_axis_labels(x_var="Parental age", y_var="Number of TR DNMs")
g.tight_layout()
g.savefig("age_effects.png", dpi=200)
