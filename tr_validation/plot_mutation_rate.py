import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss
pd.set_option("display.precision", 8)


plt.rc("font", size=13)


# define some global variables we'll access for filtering
FILTER_RECURRENT = True
FILTER_TRANSMITTED = True
PER_HAPLOTYPE = True
FILTER_ELEMENT = True
FILTER_SV = True


# define the assembly we're sourcing our DNMs from
ASSEMBLY = "CHM13v2"

# define the minimum size of events we want to consider
# (i.e., if MIN_SIZE = 1, don't include interruptions)
MIN_SIZE = 1

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.element.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))

sample_ids = mutations["sample_id"].unique()

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == 2209]
mutations = mutations[np.abs(mutations["likely_denovo_size"]) >= MIN_SIZE]


# if desired, remove recurrent sites
if FILTER_RECURRENT:
    recurrent_trids = [
        k
        for k, v in mutations.groupby("trid").size().sort_values().to_dict().items()
        if v > 1
    ]
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

# if desired, remove untransmitted DNMs in the samples for whom we can assess that
if FILTER_TRANSMITTED:
    mutations = mutations[
        (~mutations["sample_id"].isin(["2189", "2216"]))
        | (mutations["children_with_denovo_allele"] != "unknown")
    ]

# if desired, filter DNMs that didn't pass our Element validation checks
if FILTER_ELEMENT:
    # get mutations
    ortho = []
    for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.element.tsv"):
        df = pd.read_csv(
            fh,
            sep="\t",
            dtype={"sample_id": str},
        )
        df = df[df["simple_motif_size"] == "STR"]
        ortho.append(df)

    ortho = pd.concat(ortho)

    ortho = ortho[np.abs(ortho["likely_denovo_size"]) >= MIN_SIZE]
    ortho = ortho[ortho["paternal_id"] == 2209]
    # make sure we're only looking at sites with max AL <= 120
    ortho["max_al"] = ortho["child_AL"].apply(lambda a: max(map(int, a.split(","))))
    ortho = ortho[ortho["max_al"] <= 120]
    ortho = ortho[ortho["validation_status"] != "no_data"]

    fail_trids = ortho[ortho["validation_status"] != "pass"]["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(fail_trids)]


# total bp in SVs in each sample
if FILTER_SV:
    # read in manual validation CSVs
    manual_sv = pd.read_csv(f"tr_validation/data/TRs.50bp.{ASSEMBLY}.tsv", sep="\t")
    pass_sv = manual_sv[manual_sv["hifi_validation"] == "YES"]

    per_sample_sv_rates = pass_sv.groupby("alt_sample_id").size().reset_index().rename(columns={0: "numerator"})
    per_sample_sv_rates["sample_id"] = per_sample_sv_rates["alt_sample_id"].apply(lambda s: alt2orig[s])
    print (per_sample_sv_rates["numerator"].sum(), per_sample_sv_rates["numerator"].mean())

    sv_denoms = denoms.groupby(["sample_id"]).agg({"denominator": "sum"}).reset_index()
    per_sample_sv_rates = per_sample_sv_rates.merge(sv_denoms)
    print (per_sample_sv_rates)
    per_sample_sv_rates["rate"] = per_sample_sv_rates["numerator"] / per_sample_sv_rates["denominator"]
    per_sample_sv_rates["ci_lo"] = per_sample_sv_rates.apply(lambda row: (ss.chi2.ppf(0.025, 2 * row["numerator"]) / 2.) / row["denominator"], axis=1)
    per_sample_sv_rates["ci_hi"] = per_sample_sv_rates.apply(lambda row: (ss.chi2.ppf(0.975, 2 * (row["numerator"] + 1)) / 2.) / row["denominator"], axis=1)
    rate_mean = per_sample_sv_rates["rate"].mean()
    rate_sem = 1.96 * ss.sem(per_sample_sv_rates["rate"].values)

    print (f"Mean SV rate = {rate_mean}, 95% CI = {rate_mean - rate_sem, rate_mean + rate_sem}")

    pass_trids = pass_sv["trid"].unique()
    mutations = mutations[
        (
            (np.abs(mutations["likely_denovo_size"] >= 50))
            & (mutations["trid"].isin(pass_trids))
        )
        | (np.abs(mutations["likely_denovo_size"]) < 50)
    ]
print ("after sv", mutations.shape)
print ("### Total BP affected by TRs")
total_bp_aff = mutations.groupby("sample_id").agg(total_bp=("likely_denovo_size", lambda s: sum(abs(s)))).reset_index()
print (total_bp_aff)


# sum denominators across chromosomes if necessary
if "chrom" in denoms:
    denoms = (
        denoms.groupby(["sample_id", "motif_size", "simple_motif_size"])
        .agg({"denominator": "sum"})
        .reset_index()
    )
    if not PER_HAPLOTYPE:
        denoms["denominator"] /= 2

per_sample_denoms = (
    denoms.groupby(
        [
            "sample_id",
            "simple_motif_size",
        ]
    )
    .agg({"denominator": "sum"})
    .reset_index()
)
per_sample_denoms = per_sample_denoms[per_sample_denoms["sample_id"].isin(mutations["sample_id"].to_list())]

per_sample_mutation_counts = (
    mutations.groupby(
        [
            "sample_id",
            "simple_motif_size",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)

per_sample = per_sample_mutation_counts.merge(per_sample_denoms)
per_sample["rate"] = per_sample["numerator"] / per_sample["denominator"]

# group by sample, calculate overall mutation rate
print ("MUTATION RATES PER SAMPLE")
per_sample_mutation_counts = (
    mutations.groupby(
        [
            "sample_id",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)

per_sample_denoms = per_sample_denoms.groupby("sample_id").agg({"denominator": "sum"}).reset_index()
per_sample_rates = per_sample_mutation_counts.merge(per_sample_denoms)

per_sample_rates["rate"] = per_sample_rates["numerator"] / per_sample_rates["denominator"]
per_sample_rates["ci_lo"] = per_sample_rates.apply(lambda row: (ss.chi2.ppf(0.025, 2 * row["numerator"]) / 2.) / row["denominator"], axis=1)
per_sample_rates["ci_hi"] = per_sample_rates.apply(lambda row: (ss.chi2.ppf(0.975, 2 * (row["numerator"] + 1)) / 2.) / row["denominator"], axis=1)
rate_mean = per_sample_rates["rate"].mean()
rate_sem = 1.96 * ss.sem(per_sample_rates["rate"].values)

print (f"Mean rate = {rate_mean}, 95% CI = {rate_mean - rate_sem, rate_mean + rate_sem}")


# group by simple motif size
per_motif_denoms = (
    denoms.groupby(["sample_id", "simple_motif_size"])
    .agg({"denominator": "sum"})
    .reset_index()
)
per_motif_mutation_counts = (
    mutations.groupby(["sample_id", "simple_motif_size"])
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)
per_motif_rates = per_motif_mutation_counts.merge(per_motif_denoms)
per_motif_rates = (
    per_motif_rates.groupby(["simple_motif_size"])
    .agg(
        num_total=("numerator", "sum"),
        denom_total=("denominator", "sum"),
    )
    .reset_index()
)

per_motif_rates["rate"] = per_motif_rates["num_total"] / per_motif_rates["denom_total"]
per_motif_rates["ci_lo"] = per_motif_rates.apply(
    lambda row: (ss.chi2.ppf(0.025, 2 * row["num_total"]) / 2.0) / row["denom_total"],
    axis=1,
)
per_motif_rates["ci_hi"] = per_motif_rates.apply(
    lambda row: (ss.chi2.ppf(0.975, 2 * (row["num_total"] + 1)) / 2.0)
    / row["denom_total"],
    axis=1,
)
print (per_motif_rates)

# GROUP BY MOTIF AND TRANSMISSION
per_motif_denoms = (
    denoms.groupby(["sample_id", "motif_size"])
    .agg({"denominator": "sum"})
    .reset_index()
)
per_motif_mutation_counts = (
    mutations.groupby(["sample_id", "motif_size"])
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)
per_motif_rates = per_motif_mutation_counts.merge(per_motif_denoms)
per_motif_rates["has_transmission"] = per_motif_rates["sample_id"].apply(lambda s: s in ("2189", "2216"))
per_motif_rates = (
    per_motif_rates.groupby(["motif_size", "has_transmission"])
    .agg(
        num_total=("numerator", "sum"),
        denom_total=("denominator", "sum"),
    )
    .reset_index()
)
print (per_motif_rates)

per_motif_rates["rate"] = per_motif_rates["num_total"] / per_motif_rates["denom_total"]
per_motif_rates["ci_lo"] = per_motif_rates.apply(lambda row: (ss.chi2.ppf(0.025, 2 * row["num_total"]) / 2.) / row["denom_total"], axis=1)
per_motif_rates["ci_hi"] = per_motif_rates.apply(lambda row: (ss.chi2.ppf(0.975, 2 * (row["num_total"] + 1)) / 2.) / row["denom_total"], axis=1)
print (per_motif_rates.groupby("has_transmission").agg({"num_total": "sum", "denom_total": "sum"}))


# plot rates across samples
f, ax = plt.subplots(figsize=(12, 5))

xcol, ycol = "motif_size", "rate"
palette = {True: "goldenrod", False: "cornflowerblue"}

params = {"linestyle": "-", "capthick":1, "ms":8, "capsize":2, "mec":"w", "fmt":"o"}

for tm, tm_df in per_motif_rates.query("motif_size >= 1").groupby("has_transmission"):

    x, y = tm_df[xcol], tm_df[ycol].values
    err = tm_df[["ci_lo", "ci_hi"]].values.T
    offsets = np.abs(err - y[None, :])

    ax.fill_between(
        x,
        y1=err[0],
        y2=err[1],
        color=palette[tm],
        alpha=0.1,
    )

    ax.scatter(
        x,
        y,
        c=palette[tm],
        ec="w",
        lw=0.5,
        label=f"Transmitted" if tm else f"Not transmitted",
    )

    ax.plot(
        x,
        y,
        c=palette[tm],
    )

# create inset ax
ax_strs = inset_axes(
    ax,
    2.5,
    1.5,
    loc="upper left",
    bbox_to_anchor=(0.65,0.3,.3,.3),
    bbox_transform=ax.transAxes,
    axes_kwargs={"yscale": "log"},
)

for tm, tm_df in per_motif_rates.query("motif_size >= 1 and motif_size <= 6").groupby("has_transmission"):

    x, y = tm_df[xcol], tm_df[ycol].values
    err = tm_df[["ci_lo", "ci_hi"]].values.T
    offsets = np.abs(err - y[None, :])

    ax_strs.fill_between(
        x,
        y1=err[0],
        y2=err[1],
        color=palette[tm],
        alpha=0.1,
    )

    ax_strs.scatter(
        x,
        y,
        c=palette[tm],
        ec="w",
        lw=0.5,
    )

    ax_strs.plot(
        x,
        y,
        c=palette[tm],
    )


ax_strs.set_xlabel("Motif size (bp)")
ax_strs.set_ylabel(None)
ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per haplotype\n per generation) +/- 95% CI")
ax.set_xlabel("Motif size (bp)")
ax.legend(fancybox=True, shadow=True, title="Validated via transmission to G4?")
sns.despine(ax=ax)
f.savefig("rates.poisson.png", dpi=200)
f.savefig("rates.poisson.pdf")
