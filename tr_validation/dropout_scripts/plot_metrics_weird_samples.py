import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob


ASSEMBLY = "GRCh38"

dfs = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.{ASSEMBLY}.tsv"):
    df = pd.read_csv(fh, sep="\t")
    dfs.append(df)

dfs = pd.concat(dfs)

dfs = dfs[dfs["paternal_id"].isin([200080, 2189])]

poi_bias_col = "Subfamily POI bias"
dfs[poi_bias_col] = dfs["paternal_id"].apply(lambda p: "maternal" if p == 200080 else "paternal")

dfs = dfs[np.abs(dfs["likely_denovo_size"]) < 100]
dfs["chrom"] = dfs["trid"].apply(lambda t: t.split("_")[0])

dfs["POI"] = dfs["phase_summary"].apply(lambda p: p.split(":")[0])
dfs = dfs[dfs["POI"] != "unknown"]
dfs = dfs[dfs["trid"].str.startswith("chr")]
dfs["dad_overlap"] = dfs["father_overlap_coverage"].apply(lambda c: sum(list(map(int, c.split(",")))))
dfs["mom_overlap"] = dfs["mother_overlap_coverage"].apply(lambda c: sum(list(map(int, c.split(",")))))
dfs["dad_overlap_frac"] = dfs["dad_overlap"] / dfs["dad_dp"]
dfs["mom_overlap_frac"] = dfs["mom_overlap"] / dfs["mom_dp"]
FEATURES = [
    "mom_dp",
    "dad_dp",
    "child_coverage",
    "denovo_al",
    "likely_denovo_size",
    "child_ratio",
    "allele_ratio",
    "denovo_coverage",
    "dad_overlap_frac",
    "mom_overlap_frac",
]
dfs_sub = dfs[[poi_bias_col, "POI"] + FEATURES]

print (dfs.groupby([poi_bias_col, "simple_motif_size"]).size())

dfs["phase_int"] = dfs["POI"].apply(lambda p: 1 if p == "dad" else 0)

f, ax = plt.subplots(figsize=(12, 4))
sns.pointplot(
    data=dfs, x="chrom", y="phase_int", hue=poi_bias_col, linestyle="none", ax=ax
)
f.savefig("phase_by_chrom.png")


dfs_sub_tidy = dfs_sub.melt(value_vars=FEATURES, id_vars=[poi_bias_col, "POI"], value_name="value", var_name="feature")

g = sns.FacetGrid(data=dfs_sub_tidy, col="feature", sharey=False, aspect=2, col_wrap=2)
g.map(sns.boxplot, "POI", "value", poi_bias_col, dodge=True, palette="colorblind", fliersize=0)
g.map(sns.stripplot, "POI", "value", poi_bias_col, palette="colorblind", dodge=True)
g.add_legend(title=poi_bias_col)
g.savefig("o.png", dpi=200)
