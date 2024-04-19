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

def gc_frac(motif: str):
    gc_bp, bp = 0, 0
    for m in motif.split(","):
        gc_bp += sum([b in ("G", "C") for b in list(m)])
        bp += len(m)
    return gc_bp / bp

ped = pd.read_csv("tr_validation/data/file_mapping.csv", dtype={"sample_id": str, "alt_sample_id": str})
smp2alt = dict(zip(ped["sample_id"], ped["alt_sample_id"]))

mutations = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# annotate samples that theoretically have transmission information
mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))
# annotate mutations with motif size

mutations["complex"] = mutations["child_MC"].str.contains("_")

mutations["motif_size"] = mutations.apply(lambda row: len(row["motif"]) if row["complex"] is False else -1, axis=1)

print (mutations.query("motif_size == 1").groupby("motif").size())


REMOVE_COMPLEX = False

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

# filter to sites that have transmission information if relevant
mutations = mutations[(mutations["has_transmission"] == False) | (mutations["children_with_denovo_allele"] != "unknown")]

print (mutations.query("motif_size == 1").groupby("motif").size())

# mutations["chrom"] = mutations["trid"].apply(lambda t: t.split("_")[0])
# mutations["start"] = mutations["trid"].apply(lambda t: t.split("_")[1])
# mutations["end"] = mutations["trid"].apply(lambda t: t.split("_")[2])
# mutations["inheritance"] = mutations["phase_summary"].apply(lambda p: "paternal" if p.split(":")[0] == "dad" else "maternal" if p.split(":")[0] == "mom" else "cannot_determine")
# mutations["children_with_denovo_allele"] = mutations["children_with_denovo_allele"].apply(lambda c: "|".join(c.split(",")))
# mutations["variant_type"] = "de novo"
# mutations["sample_id"] = mutations["sample_id"].apply(lambda s: smp2alt[s])
# mutations[
#     [
#         "sample_id",
#         "chrom",
#         "start",
#         "end",
#         "motif",
#         "variant_type",
#         "inheritance",
#         "children_with_denovo_allele",
#     ]
# ].to_csv(
#     "combined.tr.for_michelle.tsv",
#     sep="\t",
#     index=False,
# )

# read in denominators
denoms = []
for fh in glob.glob("tr_validation/csv/rates/*.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    denoms.append(df)
denoms = pd.concat(denoms)

# denoms["motif_size"] = denoms.apply(lambda row: -1 if row["is_complex"] else row["motif_size"], axis=1)

# remove complex loci if we're doing that
if REMOVE_COMPLEX:
    denoms = denoms[denoms["is_complex"] == False]

# sum across chromosomes if we've already grouped by that
if "chrom" in denoms:
    denoms = denoms.groupby(["sample_id", "motif_size"]).agg({"denominator": "sum"}).reset_index()

# bootstrap within groups to get confidence intervals
bootstrapped_counts = []
n = 100
for _ in tqdm.tqdm(range(n)):
    # bootstrap resample the mutations
    mutations_resampled = mutations.sample(frac=1, replace=True)
    mutation_counts = (
        mutations_resampled.groupby(
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
    mutation_counts["bs"] = _
    bootstrapped_counts.append(mutation_counts)

bootstrapped_counts = pd.concat(bootstrapped_counts)

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

# merge mutations with denominators
res_df = bootstrapped_counts.merge(denoms)
res_df["rate"] = res_df["numerator"] / res_df["denominator"]

# use poisson CIs
# res_df_grouped = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index().rename(columns={"numerator": "num_total", "denominator": "denom_total"})
# res_df_grouped["ci"] = res_df_grouped.apply(lambda row: 1.96 * np.sqrt(row["num_total"] / row["denom_total"]), axis=1)
# res_df_grouped["ci_lo"] = res_df_grouped.apply(lambda row: (row["num_total"] - row["ci"]) / row["denom_total"], axis=1)
# res_df_grouped["ci_hi"] = res_df_grouped.apply(lambda row: (row["num_total"] + row["ci"]) / row["denom_total"], axis=1)
# res_df_grouped["rate"] = res_df_grouped["num_total"] / res_df_grouped["denom_total"]

# res_df = res_df.merge(res_df_grouped)

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

# ax_complex = inset_axes(
#     ax,
#     1.5,
#     1.,
#     loc="upper left",
#     axes_kwargs={"yscale": "log"},
# )

# sns.barplot(data=res_df.query("motif_size == -1"), x=x, y=y, hue="has_transmission", errorbar=None, palette=palette)

# ax_complex.get_legend().remove()
# ax_complex.set_ylabel(None)
# sns.despine(ax=ax_complex)
# ax_complex.set_xticklabels(["Complex (2+ motifs)"])
# ax_complex.set_xlabel(None)

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
# f.tight_layout()
f.savefig("rates.png", dpi=200)
