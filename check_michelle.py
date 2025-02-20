import pandas as pd

michelle = pd.read_csv("T2T.michelle.tsv", sep="\t")

michelle = michelle[michelle["pass_all_filters"] == "YES"]

print (michelle.groupby("simplified_motif_structure").size())

strs = michelle[michelle["simplified_motif_structure"] == "STR"]

print (strs[strs["max_motif_size_in_locus"] > 1].groupby("element_validation_status").size())

michelle["check"] = michelle.apply(lambda row: "STR" if (row["min_motif_size_in_locus"] <= 6 and row["max_motif_size_in_locus"] <= 6) else "VNTR" if (row["min_motif_size_in_locus"] > 6 and row["max_motif_size_in_locus"] > 6) else "complex", axis=1)
print (michelle.groupby(["check", "simplified_motif_structure"]).size())
