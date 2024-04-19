import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cyvcf2
from typing import Dict, List
from utils import filter_mutation_dataframe
import glob
import datetime

plt.rc("font", size=14)


def calculate_phase_counts(df: pd.DataFrame, group_cols: List[str] = ["sample_id", "phase"]):
    phase_counts = df.groupby(group_cols).size().reset_index().rename(columns={0: "count"})
    totals = phase_counts.groupby("sample_id").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})

    phase_counts = phase_counts.merge(totals)
    phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]
    return phase_counts

GEN = "2gen"

mutations = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

mutations = filter_mutation_dataframe(
    mutations,
    remove_complex=False,
    remove_duplicates=False,
    remove_gp_ev=True,
    remove_inherited=True,
    parental_overlap_frac_max=0.05,
    denovo_coverage_min=2,
    depth_min=10,
)
mutations["pass_inf_sites"] = mutations.apply(
    lambda row: (
        ((float(row["n_upstream"]) > 1) | (float(row["n_downstream"]) > 1))
        if row["phase_summary"] != "unknown"
        else True
    ),
    axis=1,
)
mutations = mutations[mutations["pass_inf_sites"] == True]

mutations["is_transmitted"] = mutations["children_with_denovo_allele"].apply(lambda c: c != "unknown")

# # remove untransmitted
# mutations = mutations[
#     (
#         (mutations["sample_id"].isin(["2216", "2189"]))
#         & (mutations["children_with_denovo_allele"] != "unknown")
#     )
#     | (
#         (~mutations["sample_id"].isin(["2216", "2189"]))
#         & (mutations["children_with_denovo_allele"] == "unknown")
#     )
# ]

print (mutations.query("sample_id == '2189'")[["children_with_denovo_allele"]])

ped = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str},
)

sample2dad = dict(zip(ped["sample_id"], ped["paternal_id"]))
mutations["dad_sample_id"] = mutations["sample_id"].apply(lambda s: sample2dad[s])


mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0])

phase_counts = calculate_phase_counts(mutations, group_cols=["sample_id", "phase", "is_transmitted"])
ages = pd.read_csv("tr_validation/data/k20_parental_age_at_birth.csv", dtype={"sample_id": str, "UGRP Lab ID (archive)": str})

order = phase_counts.groupby("sample_id").agg({"count": "sum"}).sort_values("count").reset_index()["sample_id"].to_list()
order2idx = dict(zip(order, range(len(order))))

phase_counts["order_idx"] = phase_counts["sample_id"].apply(lambda s: order2idx[s])

ind = np.arange(phase_counts["sample_id"].nunique())

bottom = np.zeros(ind.shape[0])

f, ax = plt.subplots(figsize=(8, 5))

# for phase, phase_df in only_phased.groupby("phase"):
#     sorted_phase_df = phase_df.sort_values("order_idx")
#     vals = sorted_phase_df["fraction"].values
#     counts = sorted_phase_df["count"].values
#     ax2.barh(
#         ind,
#         vals,
#         1,
#         left=bottom,
#         ec="w",
#         lw=2,
#         color="goldenrod" if phase == "dad" else "cornflowerblue",
#         label="Paternal" if phase == "dad" else "Maternal"
#     )
#     bottom += vals
# ax2.legend(ncol=2, fancybox=True, shadow=True, loc="upper center", bbox_to_anchor=(0.5, 1.075))
# ax2.set_xlabel("Fraction of phased TR DNMs")
# sns.despine(ax=ax2)
# bottom = np.zeros(ind.shape[0])

for (phase, is_transmitted), phase_df in phase_counts.groupby(["phase", "is_transmitted"]):
    sorted_phase_df = phase_df.sort_values("order_idx")
    print (sorted_phase_df)
    vals = np.zeros(ind.shape[0])
    idxs = sorted_phase_df["order_idx"].values
    vals[idxs] = sorted_phase_df["count"].values
    ax.barh(
        ind,
        vals,
        1,
        left=bottom,
        ec="w",
        #ls=":" if is_transmitted else "-",
        hatch="x" if is_transmitted else None,
        lw=1,
        color="gainsboro" if phase == "unknown" else "goldenrod" if phase == "dad" else "cornflowerblue",
        label=f'{"Paternal" if phase == "dad" else "Maternal" if phase == "mom" else "Unknown"}' if not is_transmitted else None,
    )
    # for y_i in ind:
    #     x = vals[y_i] / 3 + bottom[y_i]
    #     y = y_i - 0.1
    #     text = str(counts[y_i])
    #     if y_i == ind.shape[0] - 1: 
    #         adj = 4 if phase in ("dad", "unknown") else 5
    #         x = vals[y_i] / adj + bottom[y_i]
    #         pref = "Mat." if phase == "mom" else "Pat." if phase == "dad" else "Unk."
    #         #text = f"{pref} (n = {counts[y_i]})"
    #         #text = f"(n = {counts[y_i]})"
    #     ax1.text(x, y, text, c="w" if phase != "unknown" else "k", size=14)
    

    bottom += vals
ax.set_yticks(ind)
ax.set_yticklabels(phase_counts.drop_duplicates("sample_id").sort_values("order_idx")["sample_id"].unique())
ax.set_ylabel("Sample ID")
ax.set_xlabel("Number of TR DNMs")
ax.legend(title="Parent-of-origin", fancybox=True, shadow=True, bbox_to_anchor=(0.75, 0.35), fontsize=12)
sns.despine(ax=ax, left=True)

now = datetime.datetime.now()
ax.text(0.5, 0.5, f"draft ({str(now).split(' ')[0]})", transform=ax.transAxes,
        fontsize=20, color='gray', alpha=0.25,
        ha='center', va='center', rotation=30)

f.tight_layout()
f.savefig("phase.png", dpi=200)


# df['min_parent_support'] = df[['mom_min_support', 'dad_min_support']].min(axis=1)
# df['min_parent_support'] = df[['mom_dp', 'dad_dp']].min(axis=1)

# minimums = np.arange(0, 1, 0.01)
# minimums = np.arange(1, 25, 1)

# smp2idx = dict(zip(df["dad_sample_id"].unique(), range(df["dad_sample_id"].nunique())))


# smp2gen = dict(zip(df["dad_sample_id"], df["generation"]))
# arr = np.zeros((df["sample_id"].nunique(), minimums.shape[0]))
# counts = np.zeros((df["sample_id"].nunique(), minimums.shape[0]))

# N_BOOT = 1000

# bootstrap = np.zeros((N_BOOT, df["sample_id"].nunique(), minimums.shape[0]))

# res = []
# for mi, m in tqdm.tqdm(enumerate(minimums)):
#     df_sub = df[(df["dad_min_support"] >= m) & (df["mom_min_support"] >= m) & (df["kid_min_support"] >= m)]
#     #df_sub = df[(df["child_coverage"] >= m) & (df["mom_dp"] >= m) & (df["dad_dp"] >= m)]

#     # loop over samples
#     for sample, sample_df in df_sub.groupby("dad_sample_id"):

#         if sample == "2214" and m == 15:
#             print (sample_df[["trid", "denovo_coverage", "allele_coverage", "allele_ratio", "per_allele_reads_child", "per_allele_reads_father", "per_allele_reads_mother"]])
#         idx = smp2idx[sample]
#         # calculate baseline paternal fraction
#         phases = sample_df["phase"].values
#         paternal_frac = np.sum(phases) / phases.shape[0]

#         arr[idx, mi] = paternal_frac
#         counts[idx, mi] = phases.shape[0]
#         # now, bootstrap N times and resample mutations each time
#         for bi in range(N_BOOT):
#             resampled_phases = np.random.choice(phases, replace=True, size=phases.shape[0])
#             paternal_frac = np.sum(resampled_phases) / resampled_phases.shape[0]
#             bootstrap[bi, idx, mi] = paternal_frac

# f, axarr = plt.subplots(3, sharey=True, sharex=True, figsize=(6, 12))
# for s, i in smp2idx.items():
#     gen = smp2gen[s]
#     ax_i = 0 if gen == "G2" else 1 if gen == "G3" else 2

#     fracs = arr[i, ]
#     lbs = np.percentile(bootstrap[:, i, :], q=2.5, axis=0)
#     ubs = np.percentile(bootstrap[:, i, :], q=97.5, axis=0)
#     axarr[ax_i].plot(minimums, fracs, label=s)
#     axarr[ax_i].fill_between(minimums, y1=lbs, y2=ubs, alpha=0.1)
#     axarr[ax_i].set_ylim(0, 1)
#     axarr[ax_i].set_xlabel("Minimum allele depth required in all members of trio")
#     axarr[ax_i].set_ylabel("Fraction of DNMs phased to father")
#     axarr[ax_i].legend(title="Paternal ID", frameon=False)
#     axarr[ax_i].set_title(gen)
#     sns.despine(ax=axarr[ax_i])

# f.tight_layout()
# f.savefig('test.png', dpi=200)
