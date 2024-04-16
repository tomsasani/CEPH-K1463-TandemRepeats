import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cyvcf2
from typing import Dict, List
from utils import filter_mutation_dataframe
import glob
from build_roc_curve_dnms import get_dp, annotate_with_allr, check_for_overlap
from plot_depth_ratios import get_binomial_p
import tqdm

plt.rc("font", size=16)


def calculate_phase_counts(df: pd.DataFrame):
    phase_counts = df.groupby(["sample_id", "phase"]).size().reset_index().rename(columns={0: "count"})
    totals = phase_counts.groupby("sample_id").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})

    phase_counts = phase_counts.merge(totals)
    phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]
    return phase_counts

GEN = "2gen"

# loop over phased mutation dataframes
res = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    res.append(df)
df = pd.concat(res).fillna("unknown")

df = filter_mutation_dataframe(
    df,
    remove_complex=False,
    remove_duplicates=True,
    remove_gp_ev=False,
    remove_inherited=True,
    parental_overlap_frac_max=0.05,
    denovo_coverage_min=2,
    depth_min=10,
)
df["pass_inf_sites"] = df.apply(lambda row: ((float(row["n_upstream"]) > 1) | (float(row["n_downstream"]) > 1)) if row["phase_summary"] != "unknown" else True, axis=1)
#df = df[(df["phase_summary"] == "unknown") | ((df["n_upstream"].astype(float) > 1) | (df["n_downstream"].astype(float) > 1))]
df = df[df["pass_inf_sites"] == True]

ped = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str},
)

sample2dad = dict(zip(ped["sample_id"], ped["paternal_id"]))
df["dad_sample_id"] = df["sample_id"].apply(lambda s: sample2dad[s])


df["phase"] = df["phase_summary"].apply(lambda p: p.split(":")[0])

phase_counts = calculate_phase_counts(df)
ages = pd.read_csv("tr_validation/data/k20_parental_age_at_birth.csv", dtype={"sample_id": str, "UGRP Lab ID (archive)": str})
phase_counts["generation"] = phase_counts["sample_id"].apply(lambda s: "G4" if str(s).startswith("200") else "G2" if s in ("2209", "2188") else "G3")

sns.set_style("ticks")

f, ax = plt.subplots(figsize=(12, 8))

order = phase_counts[phase_counts["phase"] == "dad"].sort_values("fraction", ascending=False)["sample_id"].to_list()
order2idx = dict(zip(order, range(len(order))))

phase_counts["order_idx"] = phase_counts["sample_id"].apply(lambda s: order2idx[s])
print (phase_counts)

ind = np.arange(phase_counts["sample_id"].nunique())
bottom = np.zeros(ind.shape[0])

for phase, phase_df in phase_counts.groupby("phase"):
    sorted_phase_df = phase_df.sort_values("order_idx", ascending=False)
    vals = sorted_phase_df["fraction"].values
    counts = sorted_phase_df["count"].values
    ax.barh(
        ind,
        vals,
        1,
        left=bottom,
        ec="w",
        lw=2,
        color="gainsboro" if phase == "unknown" else "goldenrod" if phase == "dad" else "cornflowerblue",
    )
    for y_i in ind:
        x = vals[y_i] / 3 + bottom[y_i]
        y = y_i - 0.1
        text = str(counts[y_i])
        if y_i == ind.shape[0] - 1: 
            adj = 4 if phase in ("dad", "unknown") else 5
            x = vals[y_i] / adj + bottom[y_i]
            pref = "Mat." if phase == "mom" else "Pat." if phase == "dad" else "Unk."
            text = f"{pref} (n = {counts[y_i]})"
        ax.text(x, y, text, c="w" if phase != "unknown" else "k", size=18)
    

    bottom += vals
ax.set_yticks(ind)
ax.set_yticklabels(phase_counts.drop_duplicates("sample_id").sort_values("order_idx", ascending=False)["sample_id"].unique())
ax.set_ylabel("Sample ID")
ax.set_xlabel("Fraction of TR DNMs")
sns.despine(ax=ax, left=True)
ax.set_title("Inferred parent-of-origin for TR DNMs")
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
