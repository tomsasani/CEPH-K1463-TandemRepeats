import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import functools
import math
import scipy.stats as ss

def find_lcm(denoms: pd.Series):
    return functools.reduce(math.lcm, denoms)

plt.rc("font", size=13)

FILTER_RECURRENT = True
FILTER_TRANSMITTED = True

ASSEMBLY = "CHM13v2"
# ASSEMBLY = "GRCh38"

mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# remove recurrent
recurrent_trids = [k for k,v in mutations.groupby("trid").size().sort_values().to_dict().items() if v > 1]

if FILTER_RECURRENT:
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

# annotate samples that theoretically have transmission information
mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))

mutations["complex"] = mutations["child_MC"].str.contains("_")
mutations["motif_size"] = mutations.apply(lambda row: len(row["motif"]) if row["complex"] is False else -1, axis=1)

mutations = mutations.query("likely_denovo_size != 0")

mutations["is_transmitted"] = mutations["children_with_denovo_allele"].apply(lambda c: c != "unknown")
print (mutations[mutations["has_transmission"]].groupby("is_transmitted").size())


# filter to sites that have transmission information if desired
if FILTER_TRANSMITTED:
    mutations = mutations[(mutations["has_transmission"] == False) | (mutations["children_with_denovo_allele"] != "unknown")]

# read in denominators
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    denoms.append(df)
denoms = pd.concat(denoms)

# sum across chromosomes if we've already grouped by that
if "chrom" in denoms:
    denoms = denoms.groupby(["sample_id", "motif_size"]).agg({"denominator": "sum"}).reset_index()

mutations = mutations.merge(denoms)
mutations["count"] = 1

mutations_grouped = mutations.groupby(["motif_size", "has_transmission"]).agg(numerator = ("count", "sum"), denominator = ("denominator", "sum")).reset_index()
mutations_grouped["rate"] = mutations_grouped["numerator"] / mutations_grouped["denominator"]

mutation_counts = (
    mutations.groupby(
        [
            "sample_id",
            "motif_size",
            "has_transmission",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "numerator"})
)

res_df = mutation_counts.merge(denoms)

res_df["rate"] = res_df["numerator"] / res_df["denominator"]
res_df["ci"] = res_df.apply(lambda row: 1.96 * np.sqrt(row["numerator"] / row["denominator"]), axis=1)
res_df["ci_lo"] = res_df.apply(lambda row: (row["numerator"] - row["ci"]) / row["denominator"], axis=1)
res_df["ci_hi"] = res_df.apply(lambda row: (row["numerator"] + row["ci"]) / row["denominator"], axis=1)


# use poisson CIs
res_df_grouped = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index().rename(columns={"numerator": "num_total", "denominator": "denom_total"})
res_df_grouped["rate"] = res_df_grouped["num_total"] / res_df_grouped["denom_total"]
res_df_grouped["ci_lo"] = res_df_grouped.apply(lambda row: (ss.chi2.ppf(0.025, 2 * row["num_total"]) / 2.) / row["denom_total"], axis=1)
res_df_grouped["ci_hi"] = res_df_grouped.apply(lambda row: (ss.chi2.ppf(0.975, 2 * (row["num_total"] + 1)) / 2.) / row["denom_total"], axis=1)

print (res_df_grouped)

boop = res_df.groupby(["sample_id", "has_transmission"]).agg({"denominator": "sum", "numerator": "sum"}).reset_index()
boop["rate"] = boop["numerator"] / boop["denominator"]
print (boop.groupby("has_transmission").agg({"rate": ["mean", "std"]}))

# plot rates across samples
f, ax = plt.subplots(figsize=(12, 5))

xcol, ycol = "motif_size", "rate"
palette = {True: "goldenrod", False: "cornflowerblue"}

params = {"linestyle": "-", "capthick":1, "ms":8, "capsize":2, "mec":"w", "fmt":"o"}

for tm, tm_df in res_df_grouped.query("motif_size >= 1 and motif_size <= 100").groupby("has_transmission"):
    print (tm_df[["motif_size", "has_transmission", "rate", "ci_lo", "ci_hi"]])

    x, y = tm_df[xcol], tm_df[ycol].values

    err = tm_df[["ci_lo", "ci_hi"]].values.T

    offsets = np.abs(err - y[None, :])

    ax.errorbar(
        x,
        y,
        yerr=offsets,
        c=palette[tm],
        mfc=palette[tm],
        **params,
        label=f"Transmitted" if tm else f"Not transmitted",
    )

# create inset ax
ax_strs = inset_axes(
    ax,
    3,
    2,
    loc="center right",
    axes_kwargs={"yscale": "log"},
)

for tm, tm_df in res_df_grouped.query("motif_size >= 1 and motif_size <= 6").groupby("has_transmission"):

    x, y = tm_df[xcol], tm_df[ycol].values

    err = tm_df[["ci_lo", "ci_hi"]].values.T

    offsets = np.abs(err - y[None, :])

    ax_strs.errorbar(
        x,
        y,
        yerr=offsets,
        c=palette[tm],
        mfc=palette[tm],
        **params,
        zorder=1,
    )


ax_strs.set_xlabel("Motif size (bp)")
ax_strs.set_ylabel(None)
# ax_strs.get_legend().remove()

# add watermark
# now = datetime.datetime.now()
# ax.text(0.5, 0.5, f"draft ({str(now).split(' ')[0]})", transform=ax.transAxes,
#         fontsize=20, color='gray', alpha=0.25,
#         ha='center', va='center', rotation=30)

ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)\n+/- 95% CI")
ax.set_xlabel("Motif size (bp)")
ax.legend(fancybox=True, shadow=True, title="Validated via transmission to G4?")
sns.despine(ax=ax)
f.savefig("rates.poisson.png", dpi=200)
f.savefig("rates.poisson.pdf")
