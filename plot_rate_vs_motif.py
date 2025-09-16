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
REMOVE_COMPLEX = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"sample_id": str})

mutations["is_complex"] = mutations["n_motifs"] > 1


# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))
orig2alt = {v:k for k,v in alt2orig.items()}

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"csv/denominators/*.{ASSEMBLY}.v4.0.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
if PER_HAPLOTYPE:
    denoms["denominator"] *= 2
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

if REMOVE_COMPLEX:
    denoms = denoms[denoms["is_complex"] == False]
    mutations = mutations[mutations["is_complex"] == False]


# figure out the total number of BP affected by TRs
print ("### Total BP affected by TRs")
total_bp_aff = mutations.groupby("sample_id").agg(total_bp=("likely_denovo_size", lambda s: sum(abs(s)))).reset_index()
print (total_bp_aff["total_bp"].sum() / mutations.shape[0])
print (total_bp_aff["total_bp"].sum() / len(sample_ids))

# compute total denominators in each sample
per_sample_denoms = (
    denoms.groupby("sample_id").agg({"denominator": "sum"}).reset_index()
)
# compute total mutation counts in each sample
per_sample_mutation_counts = (
    mutations.groupby("sample_id").size().reset_index().rename(columns={0: "numerator"})
)
per_sample_rates = per_sample_mutation_counts.merge(per_sample_denoms)
per_sample_rates["rate"] = per_sample_rates["numerator"] / per_sample_rates["denominator"]

# group by sample, calculate overall mutation rate
print ("## OVERALL MUTATION RATES ##")
per_sample_rates["rate"] = per_sample_rates["numerator"] / per_sample_rates["denominator"]
per_sample_rates["ci_lo"] = per_sample_rates.apply(lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[0], axis=1)
per_sample_rates["ci_hi"] = per_sample_rates.apply(lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[1], axis=1)
rate_mean = per_sample_rates["rate"].mean()
rate_sem = 1.96 * ss.sem(per_sample_rates["rate"].values)
print (per_sample_rates)
print (f"Mean rate = {rate_mean}, 95% CI = {rate_mean - rate_sem, rate_mean + rate_sem}")

# group denominators by simple motif size
per_motif_denoms = (
    denoms.groupby("simple_motif_size")
    .agg({"denominator": "sum"})
    .reset_index()
)
# group mutation counts by simple motif size
per_motif_mutation_counts = (
    mutations.groupby("simple_motif_size")
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)
per_motif_rates = per_motif_mutation_counts.merge(per_motif_denoms)

per_motif_rates["rate"] = per_motif_rates["numerator"] / per_motif_rates["denominator"]
per_motif_rates["ci_lo"] = per_motif_rates.apply(
    lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[0],
    axis=1,
)
per_motif_rates["ci_hi"] = per_motif_rates.apply(
    lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[1],
    axis=1,
)

print ("## PER-CLASS MUTATION RATES ##")
per_motif_rates["mean"] = per_motif_rates["numerator"] / len(sample_ids)
print (per_motif_rates["numerator"].sum() / len(sample_ids))
print (per_motif_rates)

# decide whether to group by motif size or reference length
mutations["reflen"] = mutations["end"] - mutations["start"]
max_reflen = mutations["reflen"].max()

reflen_bins = np.arange(10, max_reflen, 10)
bins, labels = pd.cut(mutations["reflen"].values, bins=reflen_bins, retbins=True)

mutations["reflen_bin"] = bins.astype(str)

def plot_mutation_rate_vs(
    mutations: pd.DataFrame,
    denoms: pd.DataFrame,
    colname: str = "motif_size",
    plot_strs: bool = False,
    outname: str = "out.png",
    include_labels: bool = True,
):

    # group by more straightforward motif size definition
    per_col_denoms = (
        denoms.groupby(colname)
        .agg({"denominator": "sum"})
        .reset_index()
    )
    
    per_col_mutation_counts = (
        mutations.groupby(colname)
        .size()
        .reset_index()
        .rename(columns={0: "numerator"})
    )
    per_col_rates = per_col_mutation_counts.merge(per_col_denoms)


    if colname == "reflen_bin":
        per_col_rates["reflen_bin"] = per_col_rates["reflen_bin"].apply(lambda b: int(b.split(",")[0].lstrip("()")))
        per_col_rates = per_col_rates.sort_values("reflen_bin")
        # get the reflen bin that contains 95% of data
        total_denom = per_col_rates["denominator"].sum()
        per_col_rates["denominator_frac"] = per_col_rates["denominator"] / total_denom
        per_col_rates["denominator_cumsum"] = np.cumsum(per_col_rates["denominator_frac"])
        
        per_col_rates = per_col_rates[per_col_rates["denominator_cumsum"] < 0.97]
    
    per_col_rates["rate"] = per_col_rates["numerator"] / per_col_rates["denominator"]
    per_col_rates["ci_lo"] = per_col_rates.apply(
        lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[0],
        axis=1,
    )
    per_col_rates["ci_hi"] = per_col_rates.apply(
        lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[1],
        axis=1,
    )

    # plot rates across samples
    if colname == "motif_size":
        # f, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
        # f = plt.figure(figsize=(18, 8))
        f, ax1 = plt.subplots(figsize=(14, 7))

        # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4], wspace=0.1)
        # ax2 = plt.subplot(gs[0, 0])
        # ax1 = plt.subplot(gs[0, 1])#, sharey=ax2)
    else:
        f, ax1 = plt.subplots(figsize=(11, 5))
    ycol = "rate"

    #per_col_rates = per_col_rates[per_col_rates["motif_size"] <= 1]
    print (per_col_rates)

    x, y = per_col_rates[colname], per_col_rates[ycol].values
    err = per_col_rates[["ci_lo", "ci_hi"]].values.T
    ax1.fill_between(x, y1=err[0], y2=err[1], color="gainsboro", alpha=0.5, zorder=0)
    ax1.plot(x, y, c="dodgerblue", zorder=1, lw=0.75)
    ax1.scatter(x, y, c="dodgerblue", ec="w", lw=0.25, s=50)

    # second axis for the denominator
    ax_denom = ax1.twinx()
    per_col_denoms = per_col_denoms[per_col_denoms[colname] <= per_col_rates[colname].max()]
    denom_y = per_col_denoms["denominator"] // mutations["sample_id"].nunique()
    if PER_HAPLOTYPE:
        denom_y //= 2

    ax_denom.plot(per_col_denoms[colname], denom_y, c="darkgrey", zorder=-1, lw=0.75)
    ax_denom.scatter(per_col_denoms[colname], denom_y, c="darkgrey", ec="w", lw=0.25, s=50)
    ax_denom.set_yscale("log")

    # subplot for per-motif type rates
    # group by more straightforward motif size definition
    if colname == "motif_size":
        per_col_denoms_simple = (
            denoms.groupby("simple_motif_size")
            .agg({"denominator": "sum"})
            .reset_index()
        )
        per_col_mutation_counts_simple = (
            mutations.groupby("simple_motif_size")
            .size()
            .reset_index()
            .rename(columns={0: "numerator"})
        )
        per_col_rates_simple = per_col_mutation_counts_simple.merge(per_col_denoms_simple)

        per_col_rates_simple["rate"] = per_col_rates_simple["numerator"] / per_col_rates_simple["denominator"]
        per_col_rates_simple["ci_lo"] = per_col_rates_simple.apply(
            lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[0],
            axis=1,
        )
        per_col_rates_simple["ci_hi"] = per_col_rates_simple.apply(
            lambda row: calculate_poisson_ci(row["numerator"], row["denominator"])[1],
            axis=1,
        )
        # x, y = [0.25, 1, 1.75], per_col_rates_simple[ycol].values
        # err = per_col_rates_simple[["ci_lo", "ci_hi"]].values.T
        # offset = np.zeros_like(err)
        # offset[0, :] = y - err[0]
        # offset[1, :] = err[1] - y
        # ax2.set_xticks([0.25, 1, 1.75])
        # ax2.set_xlim(0, 2)
        # ax2.errorbar(x, y, yerr=offset, c="dodgerblue", fmt="o", ms=10)
        # sns.despine(ax=ax2, left=False, right=True, top=True)

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    sns.despine(ax=ax1, top=True)
    sns.despine(ax=ax_denom, right=False, top=True)
    # sns.despine(ax=ax_denom, left=False, right=False, top=True)
    ax1.set_xlim(0, per_col_rates[colname].max() + 10)
    # ax_denom.set_xlim(0, per_col_rates[colname].max() + 1)

    if include_labels:
        ax1.set_ylabel("Mutation rate\n(per locus, per haplotype\n per generation) +/- 95% CI", c="dodgerblue")
        
        # ax2.set_xticklabels(per_col_rates_simple["simple_motif_size"].to_list())
        ax_denom.set_ylabel("Total number of loci in reference genome", c="darkgrey")

        if colname == "motif_size":
            ax1.set_xlabel("Minimum motif size in locus (bp)")
        else:
            ax1.set_xlabel("Length of reference allele (bp)")
    else:
        plt.setp(ax1.get_xticklabels(), visible=False)
        # plt.setp(ax2.get_xticklabels(), visible=False)

    # plt.setp(ax1.get_yticklabels(), visible=include_labels, c="dodgerblue")
    # plt.setp(ax_denom.get_yticklabels(), visible=include_labels, c="darkgrey")

    # if colname == "motif_size":
    #     plt.setp(ax2.get_yticklabels(), visible=include_labels, c="dodgerblue")

    f.tight_layout()
    f.savefig(outname, dpi=300 if outname.endswith('png') else None)

plot_mutation_rate_vs(mutations, denoms, colname="motif_size", plot_strs=False, include_labels=True, outname="motif_size.rates.png")
