import pandas as pd
import tqdm

res = []
for fh in tqdm.tqdm(snakemake.input.fhs):
    df = pd.read_csv(fh, sep="\t")
    # group
    df_grouped = df.groupby(["hprc_sample", "min_motiflen", "is_het", "is_denovo"]).size().reset_index().rename(columns={0: "count"})
    res.append(df_grouped)
res = pd.concat(res)
res.to_csv(snakemake.output.fh, sep="\t", index=False)
