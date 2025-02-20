import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
import numpy as np

ASSEMBLY = "CHM13v2"

res = []

for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    mutations = pd.read_csv(fh, sep="\t", dtype={"paternal_id": str})
    res.append(mutations)

res_df = pd.concat(res)

print (res_df.groupby("sample_id").size())

res_df["phase"] = res_df["phase_summary"].apply(lambda p: p.split(":")[0] if (p != "unknown" and int(p.split(":")[1]) > 1) else "unknown")
res_df = res_df[(res_df["phase"] != "unknown")]# & (res_df["most_common_freq"] >= 0.8)]

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
res_df = res_df[res_df["generation"] == "G3"]

moore_lm = ols('child_ratio ~ C(phase)',
                    data=res_df[res_df["generation"] == "G3"]).fit()

table = sm.stats.anova_lm(moore_lm, typ=2)
print (table)

g = sns.FacetGrid(data=res_df, row="generation", sharex=False, aspect=1.5)
g.map(sns.boxplot,
    "sample_id",
    "child_ratio",
    "phase",
    fliersize=0,
    hue_order=["dad", "mom"]
)
g.map(sns.stripplot,
    "sample_id",
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
