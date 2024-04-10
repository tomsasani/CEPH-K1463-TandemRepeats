import pandas as pd
from utils import filter_mutation_dataframe
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import seaborn.objects as so
import tqdm 

MIN_INF_SITES = 1

# read in mutations
mutations = pd.read_csv(
    "tr_validation/csv/combined.annotated_gp.transmission.element.tsv",
    sep="\t",
    dtype={"sample_id": str},
).fillna(value={"children_with_denovo_allele": "unknown"})

# merge filtered mutations with phase
phases = pd.read_csv("tr_validation/csv/phased/combined/phased.2gen.tsv", sep="\t", dtype={"sample_id": str})
phases = phases[(phases["n_upstream"] >= MIN_INF_SITES) & (phases["n_downstream"] >= MIN_INF_SITES)]

mutations = mutations.merge(phases, how="left").fillna(
    value={
        "phase_summary": "unknown",
        "n_upstream": 0,
        "n_downstream": 0,
    }
)
mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else "unknown")

print (mutations.groupby(["sample_id", "phase"]).size())

# annotate samples that theoretically have transmission information
mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))
# annotate with motif size
mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

# filter to only the DNMs that are transmitted in the relevant samples
mutations = mutations[
    (
        (mutations["has_transmission"] == True)
        & (mutations["children_with_denovo_allele"] != "unknown")
    )
    | (mutations["has_transmission"] == False)
]


REMOVE_COMPLEX = False
REMOVE_DUPLICATES = True
REMOVE_GRANDPARENTAL = True
PARENTAL_OVERLAP_MAX = 1

MIN_INF_SITES = 1

# apply baseline filters to DNMs
mutations_filtered = filter_mutation_dataframe(
    mutations,
    remove_complex=REMOVE_COMPLEX,
    remove_duplicates=REMOVE_DUPLICATES,
    remove_gp_ev=REMOVE_GRANDPARENTAL,
    parental_overlap_frac_max=PARENTAL_OVERLAP_MAX,
    denovo_coverage_min=1,
    depth_min=10,
)

print (mutations_filtered.groupby(["sample_id", "phase"]).size())

### validation with Element ###
mutations_filtered["is_validated"] = mutations_filtered["parental_overlap"].apply(lambda v: "fail" if v != "pass" else "pass")


print (mutations_filtered.groupby("sample_id").size())


mutations_validated = (
    mutations_filtered[mutations_filtered["Total read support"] >= 10]
    .groupby(["sample_id", "is_validated"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
mutations_validated_totals = mutations_validated.groupby(["sample_id"]).agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})
mutations_validated = mutations_validated.merge(mutations_validated_totals)
mutations_validated["frac"] = mutations_validated["count"] / mutations_validated["total"]

f, ax = plt.subplots()
sns.barplot(
    data=mutations_validated,
    x="sample_id",
    y="count",
    hue="is_validated",
    errorbar=None,
    ax=ax,
    ec="w",
    lw=1,
)
f.savefig("validation.png", dpi=200)
print ('done with validation')


# read in annotations file
annotations = pd.read_csv("/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz", sep="\t")
# read in site stats
site_stats = pd.read_csv("tr_validation/csv/combined.site_stats.tsv", sep="\t")

annotations["motif_size"] = annotations["motifs"].apply(lambda m: len(m))

if REMOVE_COMPLEX:
    annotations = annotations.query("n_motifs == 1")

res = []
# loop over samples, grouped by relevant categories
for (sample_id, motif_size, phase, has_transmission), motif_df in tqdm.tqdm(
    mutations_filtered.groupby(
        [
            "sample_id",
            "motif_size",
            "phase",
            "has_transmission",
        ]
    )
):
    # denominator is the number of annotated loci of the specified motif size
    denom = annotations[annotations["motif_size"] == motif_size]
    # denominator also has to take into account the sites we couldn't access
    # in this trio
    bad_sites_in_trio = site_stats[site_stats["sample_id"] == sample_id]["trid"].unique()
    denom = denom[~denom["TRid"].isin(bad_sites_in_trio)]
    # numerator is the number of unique TRID de novos in the child
    numer = motif_df.drop_duplicates(["trid", "genotype"]).shape[0]
    rate = numer / denom.shape[0]
    res.append(
        {
            "sample_id": sample_id,
            "motif_size": motif_size,
            "phase": phase,
            "has_transmission": has_transmission,
            "numerator": numer,
            "denominator": denom.shape[0],
            "rate": rate,
        }
    )

res_df = pd.DataFrame(res)

# plot rates across samples
f, ax = plt.subplots(figsize=(6, 4))
sns.pointplot(
    data=res_df.query("motif_size <= 6"),
    x="motif_size",
    y="rate",
    hue="has_transmission",
    capsize=0.05,
    lw=1,
    markersize=5,
)

ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Motif length")
ax.set_title("Autosomal, non-pathogenic STR mutation rates in K1463 (G3)")
ax.legend(frameon=False, title="Transmission to G4?")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("rates.png", dpi=200)


# aggregate across samples
res_df_totals = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index()
res_df_totals["rate"] = res_df_totals["numerator"] / res_df_totals["denominator"]
# calculate poisson confidence intervals
res_df_totals["ci"] = res_df_totals.apply(lambda row: 1.96 * np.sqrt(row["numerator"] / row["denominator"]), axis=1)
res_df_totals["ci_lo"] = res_df_totals.apply(lambda row: (row["numerator"] - row["ci"]) / row["denominator"], axis=1)
res_df_totals["ci_hi"] = res_df_totals.apply(lambda row: (row["numerator"] + row["ci"]) / row["denominator"], axis=1)

# plot rates after aggregating
f, ax = plt.subplots()
for has_transmission, sub_df in res_df_totals.query("motif_size <= 6").groupby("has_transmission"):
    ax.plot(sub_df["motif_size"], sub_df["rate"], label=has_transmission)
    ax.fill_between(x=sub_df["motif_size"], y1=sub_df["ci_lo"], y2=sub_df["ci_hi"], alpha=0.25, color="gainsboro")
ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Motif length")
ax.set_title("STR mutation rates in K1463 (G3)")
sns.despine(ax=ax)
ax.legend(frameon=False, title="Transmission to G4?")
f.tight_layout()
f.savefig("rates.mean.png", dpi=200)


f, ax = plt.subplots(figsize=(8, 4))
res_df_totals = res_df.groupby(["sample_id", "phase"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index()
res_df_totals["rate"] = res_df_totals["numerator"] / res_df_totals["denominator"]

ind = np.arange(res_df_totals["sample_id"].nunique())
smp2idx = dict(zip(res_df_totals["sample_id"].unique(), ind))
res_df_totals["ind"] = res_df_totals["sample_id"].apply(lambda s: smp2idx[s])

bottom = np.zeros(ind.shape[0])
for phase, phase_df in res_df_totals.groupby("phase"):
    #if phase == "unknown": continue
    phase_df_sorted = phase_df.sort_values("ind", ascending=True)
    vals = phase_df_sorted["numerator"].values
    ax.bar(ind, vals, 1, bottom=bottom, label=phase, ec="w", lw=1)
    bottom += vals
ax.set_xticks(ind)
ax.set_xticklabels(smp2idx.keys())
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Sample ID")
ax.legend(frameon=False)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("rates_bar.png", dpi=200)

# merge with ped
ped = pd.read_csv(
    "tr_validation/data/k20_parental_age_at_birth.csv", dtype={"UGRP Lab ID (archive)": str}
)

res_df_merged = res_df.merge(ped, left_on="sample_id", right_on="UGRP Lab ID (archive)").query("motif_size <= 6")

res_df_merged_total = res_df_merged.groupby(["sample_id", "PaAge"]).agg({"rate": "sum"}).reset_index()

f, ax = plt.subplots()
sns.regplot(data=res_df_merged_total, x="PaAge", y="rate")
sns.despine(ax=ax)
f.savefig("age_rates.png", dpi=200)
