import pandas as pd
from collections import Counter

dfs = []
for fh in snakemake.input.fhs:
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    dfs.append(df)
dfs = pd.concat(dfs)

samples_with_denovo = dfs.groupby("trid").agg(samples_with_denovo=("sample_id", lambda s: ",".join(s))).reset_index()
dfs = dfs.merge(samples_with_denovo)

dfs["generation"] = dfs["sample_id"].apply(lambda s: "G4" if s.startswith("200") else "G2" if s in ("2209", "2188") else "G3")

# figure out which TRIDs are observed as DNs
# multiple times in a single generation, rather than multiple
# times across generations
generational = (
    dfs.groupby("trid")
    .agg(
        inter=("generation", lambda g: len(set(g))),
        intra=("generation", lambda g: Counter(g).most_common()[0][1]),
    )
    .reset_index()
)

intra_trids = generational[generational["intra"] > 1]["trid"].unique()
inter_trids = generational[generational["inter"] > 1]["trid"].unique()
combined_trids = set(intra_trids).union(set(inter_trids))

print ("INTRA-GENERATIONAL", len(intra_trids))
print ("INTER-GENERATIONAL", len(inter_trids))
# output a subsetted dataframe with the information necessary to do orthogonal validation
dfs = dfs[dfs["trid"].isin(combined_trids)].drop_duplicates("trid")
dfs = dfs.merge(generational)

dfs[
    [
        "trid",
        "#chrom",
        "start",
        "end",
        "index",
        "motifs",
        "child_AL",
        "samples_with_denovo",
    ]
].drop_duplicates("trid").to_csv(snakemake.output.out, index=False, sep="\t")
