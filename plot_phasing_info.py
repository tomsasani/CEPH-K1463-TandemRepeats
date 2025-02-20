import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List
import glob

plt.rc("font", size=14)

def check_phase(p):
    if p == "unknown": return p
    else:
        support = int(p.split(":")[1])
        if support < 1:
            return "unknown"
        else:
            return p.split(":")[0]

def calculate_phase_counts(df: pd.DataFrame, group_cols: List[str] = ["sample_id", "phase"]):
    phase_counts = df.groupby(group_cols).size().reset_index().rename(columns={0: "count"})
    totals = (
        phase_counts.groupby("alt_sample_id")
        .agg({"count": "sum"})
        .reset_index()
        .rename(columns={"count": "total"})
    )

    phase_counts = phase_counts.merge(totals)
    phase_counts["fraction"] = phase_counts["count"] / phase_counts["total"]
    return phase_counts

ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

ORTHOGONAL_TECH = "element"

orig = []
# for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.tsv"):
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    orig.append(df)
mutations = pd.concat(orig)

print (mutations)


mutations = pd.read_csv(
    f"tr_validation/csv/mutations.{ASSEMBLY}.filtered.tsv",
    sep="\t",
    dtype={"sample_id": str, "paternal_id": str},
)

mutations = mutations[~mutations["alt_sample_id"].isin(["NA12877", "NA12878"])]
mutations = mutations[mutations["paternal_id"] == "2209"]

mutations["combined_phase"] = mutations["phase_summary"].apply(lambda p: check_phase(p))
mutations["ypos"] = 1

# f, ax = plt.subplots()
# res = []
# for min_size in range(10):
#     filtered = mutations[np.abs(mutations["likely_denovo_size"]) >= min_size]
#     print (min_size, filtered.shape)
#     filtered["Minimum de novo size"] = min_size
#     filtered = filtered[filtered["combined_phase"] != "unknown"]
#     filtered["phase_int"] = filtered["combined_phase"].apply(lambda p: 1 if p == "dad" else 0)
#     res.append(filtered[["alt_sample_id", "phase_int", "Minimum de novo size"]])

# res = pd.concat(res)

# sns.pointplot(
#     data=res,
#     x="Minimum de novo size",
#     y="phase_int",
#     # hue="alt_sample_id",
#     linestyles="none",
#     dodge=True,
#     errorbar=("ci", 95),
#     ax=ax,
# )
# ax.set_ylabel("Fraction of paternal DNMs\n(mean +/- 95% CI)")
# f.tight_layout()
# f.savefig("test.png", dpi=200)

mutations["is_transmitted"] = mutations.apply(
    lambda row: (
        "Y"
        if (
            row["sample_id"] in ("2216", "2189")
            and row["children_with_denovo_allele"] != "unknown"
        )
        else "?" if (row["sample_id"] not in ("2189", "2216"))
        else "N"
    ),
    axis=1,
)

ped = pd.read_csv(
    "tr_validation/data/file_mapping.csv",
    dtype={"sample_id": str},
)

sample2dad = dict(zip(ped["sample_id"], ped["paternal_id"]))

# remove unphased
REMOVE_UNPHASED = True
if REMOVE_UNPHASED:
    mutations = mutations[~mutations["combined_phase"].isin(["unknown", "conflicting"])]


phase_counts = calculate_phase_counts(
    mutations,
    group_cols=["alt_sample_id", "paternal_id", "is_transmitted", "combined_phase"],
)

phase_counts["generation"] = phase_counts["paternal_id"].apply(
    lambda s: (
        "G4A"
        if s == "200080"
        else (
            "G4B"
            if s == "2189"
            else "G3" if s == "2209" else "G2A" if s == "2281" else "G2B"
        )
    )
)


phase_counts_totals = calculate_phase_counts(mutations, group_cols=["alt_sample_id", "combined_phase"])
ages = pd.read_csv("tr_validation/data/k20_parental_age_at_birth.csv", dtype={"sample_id": str, "UGRP Lab ID (archive)": str})

order = phase_counts.groupby("alt_sample_id").agg({"count": "sum"}).sort_values("count").reset_index()["alt_sample_id"].to_list()
order = phase_counts.sort_values("generation")["alt_sample_id"].unique()
# order = [o for o in order]
order2idx = dict(zip(list(map(str, order)), range(len(order))))


# order according to michelle's diagram
# order2idx = dict(
#     zip(
#         (
#             "NA12879",
#             "NA12886",
#             "NA12881",
#             "NA12882",
#             "NA12883",
#             "NA12884",
#             "NA12885",
#             "NA12887",
#         ),
#         list(range(8))[::-1],
#     )
# )

phase_counts["order_idx"] = phase_counts["alt_sample_id"].apply(lambda s: order2idx[str(s)])

ind = np.arange(phase_counts["alt_sample_id"].nunique())
bottom = np.zeros(ind.shape[0])

f, ax = plt.subplots(figsize=(7, 6))


phase_totals = phase_counts.groupby(["alt_sample_id"]).agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})
phase_fracs = phase_counts.groupby(["alt_sample_id", "combined_phase"]).agg({"count": "sum"}).reset_index()
phase_fracs = phase_fracs.merge(phase_totals)
phase_fracs["frac"] = phase_fracs["count"] / phase_fracs["total"]

print (phase_fracs.query('combined_phase == "dad"')["frac"].mean())

for (phase, transmission), phase_df in phase_counts.groupby(["combined_phase", "is_transmitted"]):
    
    sorted_phase_df = phase_df.sort_values("order_idx")
    # print (sorted_phase_df)
    vals = np.zeros(ind.shape[0])
    idxs = sorted_phase_df["order_idx"].values
    vals[idxs] = sorted_phase_df["count"].values

    label = f'{"Paternal" if phase == "dad" else "Maternal" if phase == "mom" else "Unknown"}'
    if transmission != "?": label = None

    ax.barh(
        ind,
        vals,
        0.75,
        left=bottom,
        ec="w",
        lw=1,
        alpha=1,# if transmission in ("Y") else 0.75,
        color="gainsboro" if phase == "unknown" else "#104E8B" if phase == "dad" else "firebrick",
        hatch="//" if transmission == "Y" else None,
        label=label,
    )
    

    bottom += vals

# ax.set_xlim(0, 300)
ax.set_yticks(ind)
ax.set_yticklabels(phase_counts.drop_duplicates("alt_sample_id").sort_values("order_idx")["alt_sample_id"].unique())
ax.set_ylabel("Sample ID")
ax.set_xlabel("Number of phased TR DNMs")
ax.legend(title="Parent-of-origin", fancybox=True, shadow=True, fontsize=12)#, loc="center left", bbox_to_anchor=(0.7, 0.4))
ax.tick_params(axis='y', which='both',length=0)
sns.despine(ax=ax, left=True)


f.tight_layout()
f.savefig(f"phase.{ASSEMBLY}.{'only_phased' if REMOVE_UNPHASED else None}.png", dpi=200)
# f.savefig(f"phase.{ASSEMBLY}.{'only_phased' if REMOVE_UNPHASED else None}.pdf")
