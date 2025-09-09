import pandas as pd
from cyvcf2 import VCF


mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"sample_id": str})

vcf_iters = [VCF(fh) for fh in snakemake.input.other_vcfs]

gen_to_query = snakemake.params.generation_to_query

res = []
for i, row in mutations.iterrows():
    row_dict = row.to_dict()
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    denovo_al = int(row["child_AL"].split(",")[row["index"]])
    region = f"{chrom}:{start}-{end}"
    others_with_denovo_allele = []
    for vcf in vcf_iters:
        sample = vcf.samples[0]
        for v in vcf(region):
            trid = v.INFO.get("TRID")
            if trid != row["trid"]: continue
            allele_lengths = list(v.format("AL")[0])
            if denovo_al in allele_lengths:
                others_with_denovo_allele.append(sample)
    row_dict.update({f"{gen_to_query}_with_denovo_allele": ",".join(others_with_denovo_allele)})
    res.append(row_dict)

res = pd.DataFrame(res)
res.to_csv(snakemake.output.fh, sep="\t", index=False)


