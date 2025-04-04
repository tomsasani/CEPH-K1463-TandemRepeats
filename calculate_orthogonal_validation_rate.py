import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
import scipy.stats as ss


plt.rc("font", size=14)

TECH = "element"
ASSEMBLY = "CHM13v2"
MIN_SIZE = 1

# get mutations
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{TECH}.tsv"):
    df = pd.read_csv(
        fh,
        sep="\t",
        dtype={"sample_id": str},
    )
    df = df[df["simple_motif_size"] == "STR"]

    mutations.append(df)

mutations = pd.concat(mutations)

mutations = mutations[mutations["paternal_id"] == 2209]
mutations = mutations[np.abs(mutations["likely_denovo_size"]) >= MIN_SIZE]
mutations["is_phased"] = mutations["phase_consensus"].apply(lambda p: "Y" if p != "unknown" else "N")
mutations["is_phased_int"] = mutations["is_phased"].apply(lambda p: 1 if p == "Y" else 0)

# FILTERING
mutations["max_al"] = mutations["child_AL"].apply(lambda a: max(map(int, a.split(","))))
mutations = mutations[mutations["max_al"] <= 120]
mutations = mutations[mutations["validation_status"] != "no_data"]
mutations["is_homopolymer"] = mutations["max_motiflen"] == 1

print (mutations.groupby(["is_homopolymer", "validation_status"]).size())

# fisher's exact test to see if phased stuff is more likely to be validated
a_fore = mutations.query("validation_status == 'pass' and is_phased == 'Y'").shape[0]
a_back = mutations.query("validation_status == 'fail' and is_phased == 'Y'").shape[0]
b_fore = mutations.query("validation_status == 'pass' and is_phased == 'N'").shape[0]
b_back = mutations.query("validation_status == 'fail' and is_phased == 'N'").shape[0]

print (a_fore, a_back, b_fore, b_back)
print (ss.fisher_exact([[a_fore, a_back], [b_fore, b_back]], alternative="greater"))

mutations["is_validated"] = mutations["validation_status"].apply(lambda v: 1 if v == "pass" else 0)

f, ax = plt.subplots()
sns.pointplot(
    data=mutations,
    x="is_phased",
    y="is_validated",
    hue="is_homopolymer",
    ax=ax,
    # dodge=True,
    linestyle="none",
)
ax.set_xlabel("DNM has a confident POI?")
ax.set_ylabel("Fraction of DNMs\nconsistent with Element data")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("ortho.png", dpi=200)