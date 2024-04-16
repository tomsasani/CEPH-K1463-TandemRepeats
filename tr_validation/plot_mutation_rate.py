import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes, inset_axes
import seaborn as sns
import utils

mutations = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.tsv"):
    df = pd.read_csv(fh, sep="\t")
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown"})

denoms = []
for fh in glob.glob("tr_validation/csv/rates/*.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t")
    denoms.append(df)
denoms = pd.concat(denoms)

res_df = mutations.merge(denoms)
print (res_df)

res_df["numerator"] = 1

res_df_grouped = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index().rename(columns={"numerator": "num_total", "denominator": "denom_total"})
res_df_grouped["ci"] = res_df_grouped.apply(lambda row: 1.96 * np.sqrt(row["num_total"] / row["denom_total"]), axis=1)
res_df_grouped["ci_lo"] = res_df_grouped.apply(lambda row: (row["num_total"] - row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["ci_hi"] = res_df_grouped.apply(lambda row: (row["num_total"] + row["ci"]) / row["denom_total"], axis=1)
res_df_grouped["rate"] = res_df_grouped["num_total"] / res_df_grouped["denom_total"]
res_df = res_df.merge(res_df_grouped)

# plot rates across samples
f, ax = plt.subplots(figsize=(12, 5))

x, y = "motif_size", "rate"

for group, group_df in res_df_grouped.query("motif_size <= 100").groupby("has_transmission"):
    group_df = group_df.sort_values("motif_size", ascending=True)
    ax.plot(
        group_df[x],
        group_df[y],
        label="Transmitted to G4" if group else "Untransmitted",
        lw=2,
    )
    # ax.scatter(
    #     group_df[x],
    #     group_df[y],
    # )
    ax.fill_between(
        group_df[x],
        group_df["ci_lo"],
        group_df["ci_hi"],
        color="gainsboro",
        alpha=0.5,
    )

axins = inset_axes(
    ax,
    3,
    2,
    loc="center right",
    axes_kwargs={"yscale": "log"},
)

for group, group_df in res_df_grouped.query("motif_size <= 6").groupby("has_transmission"):
    group_df = group_df.sort_values("motif_size", ascending=True)
    axins.plot(
        group_df[x],
        group_df[y],
        lw=2,
    )
    axins.fill_between(
        group_df[x],
        group_df["ci_lo"],
        group_df["ci_hi"],
        color="gainsboro",
        alpha=0.5,
    )
axins.set_xlim(1, 6)
axins.set_xlabel("Motif size (bp)")
# axins.set_title("Short tandem repeats (STRs)")

# mark_inset(ax, axins, loc1=1, loc2=3, fc="none", lw=1, ec="k")#, ec="0.5")

ax.text(0.5, 0.5, 'draft (2024/4/11)', transform=ax.transAxes,
        fontsize=20, color='gray', alpha=0.25,
        ha='center', va='center', rotation=30)

ax.set_yscale("log")
ax.set_ylabel("Mutation rate\n(per locus, per generation)")
ax.set_xlabel("Motif size (bp)")
ax.set_title("Autosomal, non-pathogenic STR mutation rates in K1463 (G3)")
ax.legend(frameon=False)#, title="Transmission to G4?")
sns.despine(ax=ax)
f.savefig("rates.png", dpi=200)


# # aggregate across samples
# res_df_totals = res_df.groupby(["motif_size", "has_transmission"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index()
# res_df_totals["rate"] = res_df_totals["numerator"] / res_df_totals["denominator"]
# # calculate poisson confidence intervals
# res_df_totals["ci"] = res_df_totals.apply(lambda row: 1.96 * np.sqrt(row["numerator"] / row["denominator"]), axis=1)
# res_df_totals["ci_lo"] = res_df_totals.apply(lambda row: (row["numerator"] - row["ci"]) / row["denominator"], axis=1)
# res_df_totals["ci_hi"] = res_df_totals.apply(lambda row: (row["numerator"] + row["ci"]) / row["denominator"], axis=1)

# # plot rates after aggregating
# f, ax = plt.subplots()
# for has_transmission, sub_df in res_df_totals.query("motif_size <= 6").groupby("has_transmission"):
#     ax.plot(sub_df["motif_size"], sub_df["rate"], label=has_transmission)
#     ax.fill_between(x=sub_df["motif_size"], y1=sub_df["ci_lo"], y2=sub_df["ci_hi"], alpha=0.25, color="gainsboro")
# ax.set_yscale("log")
# ax.set_ylabel("Mutation rate\n(per locus, per generation)")
# ax.set_xlabel("Motif length")
# ax.set_title("STR mutation rates in K1463 (G3)")
# sns.despine(ax=ax)
# ax.legend(frameon=False, title="Transmission to G4?")
# f.tight_layout()
# f.savefig("rates.mean.png", dpi=200)


# # f, ax = plt.subplots(figsize=(8, 4))
# # res_df_totals = res_df.groupby(["sample_id", "phase"]).agg({"numerator": "sum", "denominator": "sum"}).reset_index()
# # res_df_totals["rate"] = res_df_totals["numerator"] / res_df_totals["denominator"]

# # ind = np.arange(res_df_totals["sample_id"].nunique())
# # smp2idx = dict(zip(res_df_totals["sample_id"].unique(), ind))
# # res_df_totals["ind"] = res_df_totals["sample_id"].apply(lambda s: smp2idx[s])

# # bottom = np.zeros(ind.shape[0])
# # for phase, phase_df in res_df_totals.groupby("phase"):
# #     #if phase == "unknown": continue
# #     phase_df_sorted = phase_df.sort_values("ind", ascending=True)
# #     vals = phase_df_sorted["numerator"].values
# #     ax.bar(ind, vals, 1, bottom=bottom, label=phase, ec="w", lw=1)
# #     bottom += vals
# # ax.set_xticks(ind)
# # ax.set_xticklabels(smp2idx.keys())
# # ax.set_ylabel("Mutation rate\n(per locus, per generation)")
# # ax.set_xlabel("Sample ID")
# # ax.legend(frameon=False)
# # sns.despine(ax=ax)
# # f.tight_layout()
# # f.savefig("rates_bar.png", dpi=200)

# # # merge with ped
# # ped = pd.read_csv(
# #     "tr_validation/data/k20_parental_age_at_birth.csv", dtype={"UGRP Lab ID (archive)": str}
# # )

# # res_df_merged = res_df.merge(ped, left_on="sample_id", right_on="UGRP Lab ID (archive)").query("motif_size <= 6")

# # res_df_merged_total = res_df_merged.groupby(["sample_id", "PaAge"]).agg({"rate": "sum"}).reset_index()

# # f, ax = plt.subplots()
# # sns.regplot(data=res_df_merged_total, x="PaAge", y="rate")
# # sns.despine(ax=ax)
# # f.savefig("age_rates.png", dpi=200)