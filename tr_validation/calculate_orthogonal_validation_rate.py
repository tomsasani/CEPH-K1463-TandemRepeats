import pandas as pd
from utils import filter_mutation_dataframe
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance
import utils
import numpy as np

plt.rc("font", size=14)

ASSEMBLY = "GRCh38"
# ASSEMBLY = "CHM13v2"

# read in denominators
# denom = []
# for fh in glob.glob(f"tr_validation/data/all_loci/orthogonal_support/*.{ASSEMBLY}.element.read_support.csv"):
#     df = pd.read_csv(
#         fh,
#         dtype={"sample_id": str},
#     )
#     df = df[df["motif"].str.len() == 1]
#     denom.append(df)

# denom = pd.concat(denom)
# print (denom.groupby("sample_id").size())

# get mutations
mutations = []
for fh in glob.glob(f"tr_validation/data/denovos/orthogonal_support/*.{ASSEMBLY}.element.read_support.csv"):
    df = pd.read_csv(
        fh,
        dtype={"sample_id": str},
    )
    # df = df[df["motif"].str.len() == 1]

    mutations.append(df)

mutations = pd.concat(mutations)

mutations = utils.filter_mutation_dataframe(
    mutations,
    remove_complex=False,
    remove_duplicates=False,
    remove_gp_ev=False,
    remove_inherited=True,
    parental_overlap_frac_max=0.1,
    denovo_coverage_min=2,
    child_ratio_min=0.2,
)

print (mutations.groupby("sample_id").size())


# annotate mutations with validation status
res = []
for i, row in mutations.iterrows():
    row_dict = row.to_dict()
    parental_overlap = annotate_with_concordance(row)
    row_dict.update({"validation_status": parental_overlap})
    res.append(row_dict)

res = pd.DataFrame(res)
print (res[["sample_id", "region", "exp_allele_diff_denovo", "exp_allele_diff_non_denovo", "validation_status", "kid_evidence"]])#, "dad_evidence", "mom_evidence", "kid_evidence"]])



f, ax = plt.subplots()
bins = np.arange(1, 50)
for val, val_df in res.groupby("validation_status"):
    ax.hist(val_df["denovo_coverage"], bins=bins, label=val)
ax.legend()
f.savefig('o.png')

# # annotate mutations with motif size
# mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

# # apply baseline filters to DNMs
# mutations_filtered = filter_mutation_dataframe(
#     mutations,
#     remove_complex=False,
#     remove_duplicates=False,
#     remove_gp_ev=True,
#     remove_inherited=True,
#     parental_overlap_frac_max=1,
#     denovo_coverage_min=2,
# )

# # output to CSV
# mutations_filtered[
#     (mutations_filtered["sample_id"] == "2189")
#     & (mutations_filtered["motif_size"] == 1)
# ].sort_values("validation_status").to_csv("2211.homopolymers.csv", index=False)


# mutations_filtered_counts = (
#     mutations_filtered
#     .groupby(["motif_size", "validation_status"])
#     .size()
#     .reset_index()
#     .rename(columns={0: "count"})
# )


# f, ax = plt.subplots(figsize=(12, 5))
# sns.barplot(data=mutations_filtered_counts.query("motif_size <= 20"),
#     x="motif_size",
#     y="count",
#     hue="validation_status",
#     errorbar=None,
#     ax=ax,
#     ec="w",
#     lw=1,
# )
# ax.legend(frameon=False)
# f.savefig("validation.png", dpi=200)
# f.tight_layout()
