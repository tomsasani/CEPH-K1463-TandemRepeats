import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes, inset_axes
import seaborn as sns
import utils
import datetime
from collections import Counter
import tqdm

plt.rc("font", size=13)

mutations = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# annotate samples that theoretically have transmission information
mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))

mutations["complex"] = mutations["child_MC"].str.contains("_")
mutations["motif_size"] = mutations.apply(lambda row: len(row["motif"]) if row["complex"] is False else -1, axis=1)


REMOVE_COMPLEX = False
FILTER_TRANSMITTED = True

## filter mutations
mutations = utils.filter_mutation_dataframe(
    mutations,
    remove_complex=REMOVE_COMPLEX,
    remove_duplicates=False,
    remove_gp_ev=True,
    remove_inherited=True,
    parental_overlap_frac_max=1,
    denovo_coverage_min=2,
    depth_min=10,
)

# filter to sites that have transmission information if desired
if FILTER_TRANSMITTED:
    mutations = mutations[(mutations["has_transmission"] == False) | (mutations["children_with_denovo_allele"] != "unknown")]

# read in denominators
denoms = []
for fh in glob.glob("tr_validation/csv/rates/*.denominator.tsv"):
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


# use poisson CIs
res_df_grouped = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index().rename(columns={"numerator": "num_total", "denominator": "denom_total"})
res_df_grouped["ci"] = res_df_grouped.apply(lambda row: 1.96 * np.sqrt(row["num_total"] / row["denom_total"]), axis=1)
res_df_grouped["ci_lo"] = res_df_grouped.apply(lambda row: (row["num_total"] - row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["ci_hi"] = res_df_grouped.apply(lambda row: (row["num_total"] + row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["rate"] = res_df_grouped["num_total"] / res_df_grouped["denom_total"]

res_df = res_df.merge(res_df_grouped)

# plot rates across samples
f, ax = plt.subplots(figsize=(12, 5))

x, y = "motif_size", "rate"
palette = {True: "goldenrod", False: "cornflowerblue"}

sns.lineplot(
    data=res_df.query("motif_size <= 100 and motif_size >= 1"),
    x=x,
    y=y,
    hue="has_transmission",
    errorbar=("ci", 95),
    ax=ax,
    palette=palette,
)

# create inset ax
ax_strs = inset_axes(
    ax,
    3,
    2,
    loc="center right",
    axes_kwargs={"yscale": "log"},
)

sns.lineplot(
    data=res_df.query("motif_size <= 6 and motif_size >= 1"),
    x=x,
    y=y,
    hue="has_transmission",
    errorbar=("ci", 95),
    ax=ax_strs,
    palette=palette,
)

ax_strs.set_xlabel("Motif size (bp)")
ax_strs.set_ylabel(None)
ax_strs.get_legend().remove()

# add watermark
now = datetime.datetime.now()
ax.text(0.5, 0.5, f"draft ({str(now).split(' ')[0]})", transform=ax.transAxes,
        fontsize=20, color='gray', alpha=0.25,
        ha='center', va='center', rotation=30)

ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Motif size (bp)")
ax.legend(fancybox=True, shadow=True, title="Validated via transmission to G4?")
sns.despine(ax=ax)
f.savefig("rates.png", dpi=200)
