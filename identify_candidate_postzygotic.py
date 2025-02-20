import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
import numpy as np

ASSEMBLY = "CHM13v2"
# ASSEMBLY = "GRCh38"

res = []

# for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
for fh in glob.glob(f"tr_validation/csv/phased/*.{ASSEMBLY}.phased.3gen.tsv"):   
    mutations = pd.read_csv(fh, sep="\t", dtype={"paternal_id": str})
    res.append(mutations)

res_df = pd.concat(res)

res_df = res_df[res_df["sample_id"].isin([2189, 2216])]

print(
    res_df[(res_df["most_common_freq"] >= 0.75) & (res_df["child_ratio"] < 0.3) & (res_df["sample_id"] == 2189)][
        [
            "#chrom",
            "start",
            "end",
            "denovo_al",
            "child_AL",
            "father_AL",
            "mother_AL",
            "child_ratio",
            "n_inf",
            "phase",
            "children_with_denovo_allele",
        ]
    ]
)

f, ax = plt.subplots()
ax.hist(res_df["most_common_freq"].values, bins=np.arange(0, 1, 0.05))
f.savefig("most_common_freq.png")


print (res_df.groupby("#chrom").size())
res_df = res_df[res_df["#chrom"] != "chrX"]

res_df["phase_2gen"] = res_df["phase_summary"].apply(lambda p: p.split(":")[0])
res_df = res_df[res_df["n_inf"] > 1]
res_df = res_df[(res_df["phase"] != "unknown")]
res_df = res_df[res_df["most_common_freq"] > 0.5]

print (res_df.groupby(["candidate_postzygotic", "phase"]).size())

res_df["generation"] = res_df["paternal_id"].apply(
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

moore_lm = ols('child_ratio ~ C(phase) + C(candidate_postzygotic)',
                    data=res_df[res_df["generation"].isin(["G3", "G2A", "G2B"])]).fit()

table = sm.stats.anova_lm(moore_lm, typ=2)
print (table)

g = sns.FacetGrid(data=res_df, row="sample_id", sharex=False, aspect=1.5)
# g.map(sns.boxplot,
#     "candidate_postzygotic",
#     "child_ratio",
#     "phase",
#     color="white",
#     fliersize=0,
#     hue_order=["dad", "mom"]
# )
g.map(sns.stripplot,
    "candidate_postzygotic",
    "child_ratio",
    "phase",
    dodge=True,
        hue_order=["dad", "mom"]

)
g.add_legend()
# sns.despine(ax=ax)
# ax.set_ylim(0, 1)
# f.tight_layout()
g.savefig("o.png")

# f, ax = plt.subplots()
# sns.boxplot(
#     data=res_df,
#     x="candidate_postzygotic",
#     y="child_ratio",
#     hue="phase",
#     palette="colorblind",
# boxprops=dict(alpha=0.5),    fliersize=0,
# )
# sns.stripplot(
#     data=res_df,
#     x="candidate_postzygotic",
#     y="child_ratio",
#     hue="phase",
#     dodge=True,
#     palette="colorblind",
# )
# # ax.set_ylabel("Fraction of reads supporting\nthe alternate allele")
# # ax.set_xlabel("Sample ID")
# handles, labels = ax.get_legend_handles_labels()

# # When creating the legend, only use the first two elements
# # to effectively remove the last two.
# # l = plt.legend(handles[0:2], labels[0:2], title="Candidate post-zygotic?", loc="upper right",shadow=True)
# sns.despine(ax=ax)
# f.tight_layout()
# f.savefig("o.png", dpi=200)

# print (res_df.groupby(["sample_id", "candidate_postzygotic", "phase"]).size())
