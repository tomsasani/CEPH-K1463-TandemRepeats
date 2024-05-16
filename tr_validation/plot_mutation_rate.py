import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes, inset_axes
import seaborn as sns
import utils
import datetime
from collections import Counter
import scipy.stats as ss

plt.rc("font", size=13)

ASSEMBLY = "CHM13v2"
# ASSEMBLY = "GRCh38"

mutations = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# annotate samples that theoretically have transmission information
mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))

mutations["complex"] = mutations["child_MC"].str.contains("_")
mutations["motif_size"] = mutations.apply(lambda row: len(row["motif"]) if row["complex"] is False else -1, axis=1)

# print (mutations.groupby("chrom").size())

REMOVE_COMPLEX = False
FILTER_TRANSMITTED = True

## filter mutations
mutations = utils.filter_mutation_dataframe(
    mutations,
    remove_complex=REMOVE_COMPLEX,
    remove_duplicates=False,
    remove_gp_ev=True,
    remove_inherited=True,
    parental_overlap_frac_max=0.05,
    denovo_coverage_min=2,
    depth_min=10,
    child_ratio_min=0.2,
)


# filter to sites that have transmission information if desired
if FILTER_TRANSMITTED:
    mutations = mutations[(mutations["has_transmission"] == False) | (mutations["children_with_denovo_allele"] != "unknown")]

# read in denominators
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    denoms.append(df)
denoms = pd.concat(denoms)


# remove complex loci if we're doing that
if REMOVE_COMPLEX:
    denoms = denoms[denoms["is_complex"] == False]
# sum across chromosomes if we've already grouped by that
if "chrom" in denoms:
    denoms = denoms.groupby(["sample_id", "motif_size"]).agg({"denominator": "sum"}).reset_index()

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
res_df_grouped["ci"] = res_df_grouped.apply(lambda row: 1.96 * np.sqrt(row["num_total"] / row["denom_total"]), axis=1)
res_df_grouped["ci_lo"] = res_df_grouped.apply(lambda row: (row["num_total"] - row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["ci_hi"] = res_df_grouped.apply(lambda row: (row["num_total"] + row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["rate"] = res_df_grouped["num_total"] / res_df_grouped["denom_total"]


f, ax = plt.subplots(figsize=(12, 5))

res_df_grouped = res_df_grouped.query("motif_size >= 1 and motif_size <= 100")

for tm, tm_df in res_df_grouped.groupby("has_transmission"):
    ax.errorbar(
        tm_df["motif_size"],
        tm_df["rate"],
        yerr=[
            tm_df["rate"] - tm_df["ci_lo"],
            tm_df["ci_hi"] - tm_df["rate"],
        ],
        fmt="o",
        mec="k",
        color="cornflowerblue" if tm else "grey",
        capsize=2,
        label="transmitted" if tm else "not transmitted",
    )
    ax.plot(tm_df["motif_size"], tm_df["rate"], c="cornflowerblue" if tm else "grey")

# create inset ax
ax_strs = inset_axes(
    ax,
    3,
    2,
    loc="center right",
    axes_kwargs={"yscale": "log"},
)
res_df_grouped_strs = res_df_grouped.query("motif_size <= 6")

for tm, tm_df in res_df_grouped_strs.groupby("has_transmission"):
    ax_strs.errorbar(
        tm_df["motif_size"],
        tm_df["rate"],
        yerr=[
            tm_df["rate"] - tm_df["ci_lo"],
            tm_df["ci_hi"] - tm_df["rate"],
        ],
        fmt="o",
        capsize=2,
        mec="k",
        color="cornflowerblue" if tm else "grey",

    )

    ax_strs.plot(tm_df["motif_size"], tm_df["rate"], c="cornflowerblue" if tm else "grey")


ax.set_yscale("log")
ax.legend(frameon=True, fancybox=True, shadow=True, title="Transmission to G4?")
sns.despine(ax=ax)
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Motif size (bp)")
f.savefig("rates.poisson.png", dpi=200)

boop = res_df.groupby(["sample_id", "has_transmission"]).agg({"denominator": "sum", "numerator": "sum"}).reset_index()
boop["rate"] = boop["numerator"] / boop["denominator"]
print (boop.groupby("has_transmission").agg({"rate": ["mean", "std"]}))

# plot rates across samples
f, ax = plt.subplots(figsize=(12, 5))

x, y = "motif_size", "rate"
palette = {True: "goldenrod", False: "cornflowerblue"}

res_df_grouped = res_df.groupby(["motif_size", "has_transmission"]).agg({"rate": ["mean", lambda r: np.std(r)]}).reset_index()
res_df_grouped.columns = ["motif_size", "has_transmission", "rate", "stdev"]

params = {"linestyle": "-", "capthick":1, "ms":8, "capsize":2, "mec":"w", "fmt":"o"}

for tm, tm_df in res_df_grouped.query("motif_size >= 1 and motif_size <= 100").groupby("has_transmission"):
    ax.errorbar(
        tm_df["motif_size"],
        tm_df["rate"],
        yerr=tm_df["stdev"],
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
    ax_strs.errorbar(tm_df["motif_size"], tm_df["rate"], yerr=tm_df["stdev"], c=palette[tm], mfc=palette[tm], **params)


ax_strs.set_xlabel("Motif size (bp)")
ax_strs.set_ylabel(None)
# ax_strs.get_legend().remove()

# add watermark
# now = datetime.datetime.now()
# ax.text(0.5, 0.5, f"draft ({str(now).split(' ')[0]})", transform=ax.transAxes,
#         fontsize=20, color='gray', alpha=0.25,
#         ha='center', va='center', rotation=30)

ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)\n+/- STDEV")
ax.set_xlabel("Motif size (bp)")
ax.legend(fancybox=True, shadow=True, title="Validated via transmission to G4?")
sns.despine(ax=ax)
f.savefig("rates.png", dpi=200)
f.savefig("rates.pdf")
