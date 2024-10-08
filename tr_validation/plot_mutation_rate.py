import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss
from FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS

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
plt.rc("font", size=12)

# define some global variables we'll access for filtering
FILTER_TRANSMITTED = False
PER_HAPLOTYPE = True
FILTER_RECURRENT = True
PER_HAPLOTYPE = True
FILTER_ELEMENT = True
FILTER_SV = True
FILTER_VNTR = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

ORTHOGONAL_TECH = "element"

# define the minimum size of events we want to consider
# (i.e., if MIN_SIZE = 1, don't include interruptions)
MIN_SIZE = 1

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == 2209]
mutations = mutations[np.abs(mutations["likely_denovo_size"]) >= MIN_SIZE]

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

# if desired, remove recurrent sites that are recurrent in G3
if FILTER_RECURRENT:
    recurrent = pd.read_csv(f"tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv", sep="\t")
    recurrent["contains_G3"] = recurrent["samples_with_denovo"].apply(lambda samples: any([s in sample_ids for s in samples.split(",")]))
    recurrent = recurrent[recurrent["contains_G3"] == True]
    recurrent_trids = recurrent["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

# if desired, filter DNMs that didn't pass our Element validation checks.
# element validation data should already be in this dataframe as an extra annotation.
if FILTER_ELEMENT:
    # filter to STRs that had orthogonal data
    orthogonal = mutations[
        (mutations["simple_motif_size"] == "STR")
        & (mutations["validation_status"] != "no_data")
    ]
    # make sure we're only looking at sites with max AL <= 120
    orthogonal["max_al"] = orthogonal["child_AL"].apply(lambda a: max(map(int, a.split(","))))
    orthogonal = orthogonal[orthogonal["max_al"] <= 120]

    fail_trids = orthogonal[orthogonal["validation_status"] != "pass"]["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(fail_trids)]

# if desired, filter SVs that failed manual inspection
if FILTER_SV:
    mutations = mutations[
        (
            (np.abs(mutations["likely_denovo_size"]) >= 50)
            & (~mutations["trid"].isin(FAIL_SVS))
        )
        | (np.abs(mutations["likely_denovo_size"]) < 50)
    ]

# if desired, filter VNTRs that failed manual inspection
if FILTER_VNTR:
    mutations = mutations[
        (
            (mutations["motif_size"] >= 6)
            & (~mutations["trid"].isin(FAIL_VNTRS))
        )
        | (mutations["motif_size"] < 6)
    ]

# if desired, remove untransmitted DNMs in the samples for whom we can assess that
if FILTER_TRANSMITTED:

    has_transmission = mutations[mutations["sample_id"].isin(["2189", "2216"])]
    is_transmitted = has_transmission[has_transmission["children_with_denovo_allele"] != "unknown"]

    mutations = mutations[
        (~mutations["sample_id"].isin(["2189", "2216"]))
        | (mutations["children_with_denovo_allele"] != "unknown")
    ]

mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else p)

# output filtered mutations to CSV
mutations.to_csv(f"tr_validation/csv/mutations.{ASSEMBLY}.filtered.csv", index=False)


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
    plot_strs: bool = True,
    outname: str = "out.png"
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
    f, ax = plt.subplots(figsize=(11, 5))

    ycol = "rate"

    x, y = per_col_rates[colname], per_col_rates[ycol].values
    err = per_col_rates[["ci_lo", "ci_hi"]].values.T

    ax.fill_between(x, y1=err[0], y2=err[1], color="gainsboro", alpha=0.5)
    ax.scatter(x, y, c="dodgerblue", ec="w", lw=0.5)
    ax.plot(x, y, c="dodgerblue")

    if plot_strs:
        ax_new = inset_axes(
            ax,
            1.75,
            1.2,
            loc="upper left",
            bbox_to_anchor=(0.1, 0.8, .3, .3),
            bbox_transform=ax.transAxes,
            axes_kwargs={"yscale": "log"},
        )

        str_df = per_col_rates.query("motif_size >= 1 and motif_size <= 6")

        x, y = str_df["motif_size"], str_df[ycol].values
        err = str_df[["ci_lo", "ci_hi"]].values.T

        ax_new.set_yscale("log")
        ax_new.set_ylim(7e-8, 7e-5)

        ax_new.fill_between(x, y1=err[0], y2=err[1], color="gainsboro", alpha=0.5)
        ax_new.scatter(x, y, c="dodgerblue", ec="w", lw=0.5)
        ax_new.plot(x, y, c="dodgerblue")
        ax_new.set_title("STR DNMs", size=10)

    ax.set_yscale("log")
    ax.set_ylabel("Mutation rate\n(per locus, per haplotype\n per generation) +/- 95% CI")
    if colname == "motif_size":
        ax.set_xlabel("Minimum motif size in locus (bp)")
    else:
        ax.set_xlabel("Length of reference allele (bp)")
    sns.despine(ax=ax)
    # f.tight_layout()
    f.savefig(outname)#, dpi=200)

plot_mutation_rate_vs(mutations, denoms, colname="motif_size", plot_strs=True, outname="motif_size.rates.pdf")
plot_mutation_rate_vs(mutations, denoms, colname="reflen_bin", plot_strs=False, outname="reflen.rates.svg")
