import pandas as pd
from utils import filter_mutation_dataframe
import matplotlib.pyplot as plt
import seaborn as sns
import glob

plt.rc("font", size=14)

# read in mutations
mutations = []
for fh in glob.glob("tr_validation/data/orthogonal_support/*.element.tsv"):
    df = pd.read_csv(
        fh,
        sep="\t",
        dtype={"sample_id": str},
    )
    mutations.append(df)

mutations = pd.concat(mutations).dropna(subset=["motif"])
#print (mutations.groupby("motif").size().sort_values())

# annotate mutations with motif size
mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

# apply baseline filters to DNMs
mutations_filtered = filter_mutation_dataframe(
    mutations,
    remove_complex=False,
    remove_duplicates=False,
    remove_gp_ev=True,
    remove_inherited=True,
    parental_overlap_frac_max=1,
    denovo_coverage_min=2,
)

# output to CSV
mutations_filtered[
    (mutations_filtered["sample_id"] == "2189")
    & (mutations_filtered["motif_size"] == 1)
].sort_values("validation_status").to_csv("2211.homopolymers.csv", index=False)


mutations_filtered_counts = (
    mutations_filtered
    .groupby(["motif_size", "validation_status"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)


f, ax = plt.subplots(figsize=(12, 5))
sns.barplot(data=mutations_filtered_counts.query("motif_size <= 20"),
    x="motif_size",
    y="count",
    hue="validation_status",
    errorbar=None,
    ax=ax,
    ec="w",
    lw=1,
)
ax.legend(frameon=False)
f.savefig("validation.png", dpi=200)
f.tight_layout()
