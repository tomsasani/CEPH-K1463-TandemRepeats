import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
from decimal import Decimal

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def is_parent_het(row: pd.Series, with_dnm: bool = True):
    poi = row["phase"]
    if poi == "dad":
        col = "father_AL" if with_dnm else "mother_AL"
        return False if len(set(row[col].split(","))) == 1 else True
    elif poi == "mom":
        col = "mother_AL" if with_dnm else "father_AL"
        return False if len(set(row[col].split(","))) == 1 else True


rng = np.random.default_rng(42)

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"sample_id": str})

mutations = mutations[mutations["phase"] != "unknown"]

# filter out complex
mutations = mutations[mutations["n_motifs"] == 1]

mutations["dnm"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=True), axis=1)
mutations["wt"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=False), axis=1)

tidy = mutations[["trid", "sample_id", "simple_motif_size", "dnm", "wt"]].melt(id_vars=["trid", "simple_motif_size"], value_vars=["dnm", "wt"])

tidy = tidy.groupby(["simple_motif_size", "variable"]).agg(count = ("value", "sum")).reset_index()
motif_counts = mutations.groupby("simple_motif_size").size().reset_index().rename(columns={0: "total"})

tidy = tidy.merge(motif_counts)
tidy["Fraction of DNMs"] = tidy["count"] / tidy["total"]

tidy["non_count"] = tidy["total"] - tidy["count"] 

tidy["variable"] = tidy["variable"].apply(lambda v: "Transmitted DNM" if v == "dnm" else "Did not transmit DNM")

f, axarr = plt.subplots(1, 3, figsize=(8, 6), sharex=True, sharey=True)

for mi, (motif, motif_df) in enumerate(tidy.groupby("simple_motif_size")):
    motif_df = motif_df.sort_values("variable", ascending=False)

    contingency = motif_df[["count", "non_count"]].values
    
    res = ss.chi2_contingency(contingency)
    pval = "{:.1e}".format(Decimal(str(res.pvalue))) if res.pvalue != 1 else "1"

    vals = motif_df["Fraction of DNMs"]
    top = 1 - vals
    axarr[mi].bar(motif_df["variable"], vals, ec="w", lw=1, label="Heterozygous", color="cornflowerblue")
    axarr[mi].bar(motif_df["variable"], top, bottom=vals, ec="w", lw=1, label="Homozygous", color="gainsboro")

    for vi, v in enumerate(vals):
        axarr[mi].text(vi - 0.1, v / 2, f"{contingency[vi, 0]}", family="monospace", color="k")
        axarr[mi].text(vi - 0.1, v + ((1 - v) / 2.2), f"{contingency[vi, 1]}", family="monospace", color="k")
    axarr[mi].set_xticklabels(["Y", "N"])
    axarr[mi].set_title(f"{motif} loci\n" + r"$\chi^{2}$" + f" p-value = {pval}")
    sns.despine(ax=axarr[mi])
    if mi == 0:
        axarr[mi].set_ylabel("Fraction of loci")
    if mi == 1:
        # axarr[mi].legend(title="Parental genotype", shadow=True)
        axarr[mi].set_xlabel("DNM observed in parental germline?")


f.tight_layout()
f.savefig("heterozygote_effect.png", dpi=200)