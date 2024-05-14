import pandas as pd
from utils import filter_mutation_dataframe
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance


plt.rc("font", size=14)

ASSEMBLY = "GRCh38"
# ASSEMBLY = "CHM13v2"

# read in denominators
denom = []
for fh in glob.glob(f"tr_validation/data/all_loci/orthogonal_support/*.{ASSEMBLY}.element.read_support.csv"):
    df = pd.read_csv(
        fh,
        dtype={"sample_id": str},
    )
    denom.append(df)

denom = pd.concat(denom)
print (denom.groupby("sample_id").size())

# get mutations
mutations = []
for fh in glob.glob(f"tr_validation/data/denovos/orthogonal_support/*.{ASSEMBLY}.element.read_support.csv"):
    df = pd.read_csv(
        fh,
        dtype={"sample_id": str},
    )
    mutations.append(df)

mutations = pd.concat(mutations)
mutations = mutations[mutations["motif"].str.len() == 1]
print (mutations.groupby("sample_id").size())

# annotate mutations with validation status
res = []
for i, row in mutations.iterrows():
    row_dict = row.to_dict()
    parental_overlap = annotate_with_concordance(row)
    row_dict.update({"validation_status": parental_overlap})
    res.append(row_dict)

res = pd.DataFrame(res)
print (res.groupby(["sample_id", "validation_status"]).size())

f, ax = plt.subplots()
for val, val_df in res.groupby("validation_status"):
    ax.hist(val_df["denovo_coverage"], label=val)
f.savefig('o.png')

print (res[(res["sample_id"] == "2216") & (res["validation_status"] == "pass")][["trid", "motif", "exp_allele_diff_denovo", "exp_allele_diff_non_denovo", "mom_evidence", "dad_evidence", "kid_evidence"]])
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
