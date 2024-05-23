import pandas as pd
from utils import filter_mutation_dataframe
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance
import utils
import numpy as np

plt.rc("font", size=14)

TECH = "element"

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
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.element.tsv"):
    df = pd.read_csv(
        fh,
        sep="\t",
        dtype={"sample_id": str},
    )
    df = df[df["motif"].str.len() <= 6]
    mutations.append(df)

mutations = pd.concat(mutations)


mutations = mutations.query("likely_denovo_size != 0")
mutations["motif_size"] = mutations["motif"].str.len()
mutations["is_homopolymer"] = mutations["motif_size"] == 1

mutations = mutations[mutations["validation_status"] != "no_data"]
print (mutations.groupby(["motif_size", "validation_status"]).size())
print (mutations.groupby("validation_status").size())
print (mutations.groupby(["is_homopolymer", "validation_status"]).size())

print (mutations.head())