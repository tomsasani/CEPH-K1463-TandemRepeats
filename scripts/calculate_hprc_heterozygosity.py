import glob
from cyvcf2 import VCF
import pandas as pd
import tqdm

if snakemake.wildcards.KIND == "denovo":
    mutations = pd.read_csv(snakemake.input.denovos, sep="\t")
else:
    denovos = pd.read_csv(snakemake.input.denovos, sep="\t")
    mutations = pd.read_csv(snakemake.input.catalog, sep="\t", names=["#chrom", "start", "end", "info"]).drop_duplicates().sample(frac=0.005, replace=False, random_state=42)
    mutations["min_motiflen"] = mutations["info"].apply(lambda info: min([len(x) for x in info.split(";")[1].split("=")[1].split(",")]))
    mutations["trid"] = mutations["info"].apply(lambda info: info.split(";")[0].split("=")[1])
    mutations = mutations[~mutations["trid"].isin(denovos["trid"].to_list())]

mutations = mutations[["#chrom", "start", "end", "min_motiflen"]]

vcf = VCF(snakemake.input.vcf)

res = []

for i,row in tqdm.tqdm(mutations.iterrows()):
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    row_dict = row.to_dict()
    for v in vcf(f"{chrom}:{start}-{end}"):
        try:
            hap_a, hap_b, is_phased = v.genotypes[0]
        except: continue
        is_het = hap_a != hap_b
        row_dict.update({"hprc_sample": snakemake.wildcards.HPRC_SAMPLE, "is_het": is_het})
        res.append(row_dict)

res = pd.DataFrame(res)
res["is_denovo"] = 1 if snakemake.wildcards.KIND == "denovo" else 0
res.to_csv(snakemake.output.mutations, sep="\t", index=False)
