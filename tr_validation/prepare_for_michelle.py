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

t2t = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.CHM13v2.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    t2t.append(df)
t2t = pd.concat(t2t).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
t2t["assembly"] = "T2T"

hg38 = []
for fh in glob.glob(f"tr_validation/csv/filtered_and_merged/*.GRCh38.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    hg38.append(df)
hg38 = pd.concat(hg38).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
hg38["assembly"] = "hg38"

mutations = pd.concat([hg38, t2t])

mutations = mutations[mutations["paternal_id"] == 2209]

mutations["complex"] = mutations["child_MC"].str.contains("_")
mutations["motif_size"] = mutations.apply(lambda row: len(row["motif"]) if row["complex"] is False else -1, axis=1)
mutations["likely_denovo_size"] = mutations.apply(lambda row: utils.add_likely_de_novo_size(row), axis=1)


REMOVE_COMPLEX = False
FILTER_TRANSMITTED = False

## filter mutations
mutations = utils.filter_mutation_dataframe(
    mutations,
    remove_complex=REMOVE_COMPLEX,
    remove_duplicates=False,
    remove_gp_ev=False,
    remove_inherited=True,
    parental_overlap_frac_max=0.05,
    denovo_coverage_min=2,
    depth_min=10,
    child_ratio_min=0.2,
)


denovo_sizes = mutations.dropna(subset=["likely_denovo_size"]).groupby("likely_denovo_size").size().reset_index().rename(columns={0: "count"})
denovo_sizes["total"] = denovo_sizes["count"].sum()
denovo_sizes["frac"] = denovo_sizes["count"] / denovo_sizes["total"]
denovo_sizes = denovo_sizes.query("likely_denovo_size >= -20 and likely_denovo_size <= 20")

f, ax = plt.subplots()
sns.histplot(
    data=mutations.dropna(subset=["likely_denovo_size"]).query(
        "likely_denovo_size >= -20 and likely_denovo_size <= 20 and assembly == 'T2T'"
    ),
    x="likely_denovo_size",
    ax=ax,
    bins=np.arange(-20, 20, 2),
)

# ax.bar(np.arange(denovo_sizes.shape[0]), denovo_sizes["frac"], 1)
f.savefig("sizes.png")

# print (mutations.query("likely_denovo_size == 0 and assembly == 'hg38'")[["sample_id", "trid", "denovo_status", "genotype", "index", "father_AL", "mother_AL", "child_AL", "phase_summary"]])

# filter to sites that have transmission information if desired
if FILTER_TRANSMITTED:
    mutations = mutations[(mutations["has_transmission"] == False) | (mutations["children_with_denovo_allele"] != "unknown")]

mutations["chrom"] = mutations["trid"].apply(lambda t: t.split("_")[0])
mutations["start"] = mutations["trid"].apply(lambda t: t.split("_")[1])
mutations["end"] = mutations["trid"].apply(lambda t: t.split("_")[2])

# file_mapping = pd.read_csv("tr_validation/data/file_mapping.csv", sep=",", dtype={"sample_id": str})

# mutations = mutations.merge(file_mapping, on="sample_id")
# print (mutations)
mutations["variant_type"] = "de novo"
mutations["inheritance"] = mutations["phase_summary"].apply(
    lambda p: (
        "paternal"
        if p.split(":")[0] == "dad" and int(p.split(":")[1]) >= 1
        else (
            "maternal"
            if p.split(":")[0] == "mom" and int(p.split(":")[1]) >= 1
            else "cannot_determine"
        )
    )
)

print (mutations.groupby(["assembly", "sample_id"]).size())

for assembly, assembly_df in mutations.groupby("assembly"):
    assembly_df[
    [
        "alt_sample_id",
        "trid",
        "chrom",
        "start",
        "end",
        "motif",
        "likely_denovo_size",
        "variant_type",
        "inheritance",
        "children_with_denovo_allele",
    ]
].fillna({"likely_denovo_size": "?"}).rename(columns={"alt_sample_id": "sample_id"}).to_csv(f"{assembly}.michelle.tsv", sep="\t", index=False)

print (mutations.groupby(["assembly", "inheritance"]).size().sort_values())
