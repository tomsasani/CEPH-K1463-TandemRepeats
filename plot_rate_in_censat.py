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
plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

# define some global variables we'll access for filtering
PER_HAPLOTYPE = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"sample_id": str})

# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
print (sample_ids)
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))
orig2alt = {v:k for k,v in alt2orig.items()}

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"csv/denominators/*.{ASSEMBLY}.v4.0.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)

print (denoms["sample_id"].unique())

if PER_HAPLOTYPE:
    denoms["denominator"] *= 2
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

group_cols = ["simple_motif_size", "overlaps_censat"]


# compute total denominators in each sample
per_sample_denoms = (
    denoms.groupby(group_cols).agg({"denominator": "sum"}).reset_index()
)
# compute total mutation counts in each sample
per_sample_mutation_counts = (
    mutations.groupby(group_cols).size().reset_index().rename(columns={0: "numerator"})
)
print (per_sample_mutation_counts)
per_sample_rates = per_sample_mutation_counts.merge(
    per_sample_denoms,
    # left_on=["sample_id", "overlaps_censat"],
    # right_on=["sample_id", "overlaps_censat"],
)
per_sample_rates["rate"] = per_sample_rates["numerator"] / per_sample_rates["denominator"]
per_sample_rates["ci_lo"] = per_sample_rates.apply(
    lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[0],
    axis=1,
)
per_sample_rates["ci_hi"] = per_sample_rates.apply(
    lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[1],
    axis=1,
)


per_sample_rates["overlaps_censat"] = per_sample_rates["overlaps_censat"].replace(
    {
        "no": "Non-repetitive",
        "gsat": r"$\gamma$",
        "censat": "Other",
        "hsat1B": "Classical human",
        "mon": "Monomeric",
        "ct": "Centomeric transition",
        "bsat": r"$\beta$",
    }
)


f, ax = plt.subplots(figsize=(8, 6))
per_sample_rates = per_sample_rates.sort_values("rate")
order2idx = dict(
    zip(
        per_sample_rates["overlaps_censat"].unique(),
        range(per_sample_rates["overlaps_censat"].nunique()),
    )
)

per_sample_rates["idx"] = per_sample_rates["overlaps_censat"].apply(lambda s: order2idx[s])

adj = -0.2

for tr, tr_df in per_sample_rates.groupby("simple_motif_size"):

    ax.scatter(
        x=tr_df["idx"] + adj,
        y=tr_df["rate"],
        label=tr

    )
    ax.errorbar(
        x=tr_df["idx"] + adj,
        y=tr_df["rate"],
        yerr=[tr_df["rate"] - tr_df["ci_lo"], tr_df["ci_hi"] - tr_df["rate"]],
        linestyle="none",
    )
    adj += 0.2

for censat, censat_df in per_sample_rates.groupby("overlaps_censat"):

    ax.text(
        censat_df["idx"].values[0] - 0.15,
        per_sample_rates["ci_hi"].max() + 0.1,
        censat_df["numerator"].sum(),
        family="monospace",
    )

ax.legend(frameon=True, shadow=True, title="Motifs in TR locus", loc="upper left")
ax.set_xticks(range(len(order2idx)))
ax.set_xticklabels(order2idx.keys(), rotation=45)
ax.set_yscale("log")
ax.set_ylabel("Mutation rate (per locus, per haplotype\nper generation) +/- 95% CI")
ax.set_xlabel("CenSat overlap")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("censat.png", dpi=200)

# chi-square test
totals = per_sample_rates.groupby("overlaps_censat").agg(total_denom=("denominator", "sum"), total_num=("numerator", "sum")).reset_index()
a_back, a_fore = totals[totals["overlaps_censat"] == 0].values[0][1:]
b_back, b_fore = totals[totals["overlaps_censat"] == 1].values[0][1:]

a_back -= a_fore
b_back -= b_fore

print (a_fore, b_fore, a_back, b_back)

print (ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]]))
