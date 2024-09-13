import argparse
import pandas as pd
from cyvcf2 import VCF


def main(args):
    mutations = pd.read_csv(args.mutations, sep="\t")

    dad_vcf = VCF(args.dad_vcf, gts012=True)
    mom_vcf = VCF(args.mom_vcf, gts012=True)
    kid_vcf = VCF(args.kid_vcf, gts012=True)

    res = []
    for i, row in mutations.iterrows():
        trid = row["trid"]
        try:
            chrom, start, end, _ = trid.split("_")
        except: continue

        region = f"{chrom}:{start}-{end}"

        row_dict = row.to_dict()

        for vcf, colname in zip((dad_vcf, mom_vcf, kid_vcf), ("father_AL", "mother_AL", "child_AL")):
            for v in vcf(region):
                allele_lengths = v.format("AL")[0]
            row_dict.update({colname: ",".join(list(map(str, allele_lengths)))})    
        res.append(row_dict)

    res_df = pd.DataFrame(res)
    res_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")    
    p.add_argument("--dad_vcf")
    p.add_argument("--kid_vcf")    
    p.add_argument("--mom_vcf")
    p.add_argument("--output")
    args = p.parse_args()
    main(args)