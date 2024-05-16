import pandas as pd
import argparse
import cyvcf2
import numpy as np


def query_al_at_trid(trid: str, vcf: cyvcf2.VCF, min_depth: int = 10):
    chrom, start, end, _ = trid.split("_")
    region = f"{chrom}:{start}-{end}"

    allele_lengths = None
    for v in vcf(region):
        if trid != v.INFO.get("TRID"): continue
        spanning_reads = v.format("SD")
        if np.sum(spanning_reads) <= min_depth: continue
        
        allele_lengths = v.format("AL")
    return allele_lengths

def query_motif_at_trid(trid: str, vcf: cyvcf2.VCF):
    chrom, start, end, _ = trid.split("_")
    region = f"{chrom}:{start}-{end}"

    motif = None
    for v in vcf(region):
        if trid != v.INFO.get("TRID"): continue
        motif = v.INFO.get("MOTIFS")
    return motif


def annotate_with_gp(
    row: pd.Series,
    pgf: cyvcf2.VCF,
    pgm: cyvcf2.VCF,
    mgf: cyvcf2.VCF,
    mgm: cyvcf2.VCF,
):
    kid_als = list(map(int, row["child_AL"].split(",")))
    kid_idx = row["index"]
    denovo_al = kid_als[kid_idx]
    trid = row["trid"]
    chrom = trid.split("_")[0]

    res = []

    for l, v in zip(("pgf", "pgm", "mgf", "mgm"), (pgf, pgm, mgf, mgm)):
        # get allele lengths at this TRID in the VCF of interest
        # if we're on the X chromosome, we can require min depth of 10 in the
        # grandmothers and half that in the grandfathers
        v_als = query_al_at_trid(
            trid,
            v,
            min_depth=5 if l in ("pgf", "mgf") and chrom in ("chrX", "chrY") else 10,
        )
        if v_als is None:
            res.append(f"{l}_?")
        else:
            if denovo_al in v_als[0]:
                res.append(f"{l}_Y")
            else:
                res.append(f"{l}_N")
    return "|".join(res)


def get_suffix(ped: pd.DataFrame, sample_id: str):
    return ped[ped["sample_id"] == sample_id]["suffix"].values[0]


def main(args):

    # read in ped file
    ped = pd.read_csv(
        args.ped,
        sep=",",
        dtype={"paternal_id": str, "maternal_id": str, "sample_id": str},
    )
    # figure out the grandparents of this sample
    dad, mom = (
        ped[ped["sample_id"] == args.focal]["paternal_id"].values[0],
        ped[ped["sample_id"] == args.focal]["maternal_id"].values[0],
    )
    dad_dad, dad_mom = (
        ped[ped["sample_id"] == dad]["paternal_id"].values[0],
        ped[ped["sample_id"] == dad]["maternal_id"].values[0],
    )
    mom_dad, mom_mom = (
        ped[ped["sample_id"] == mom]["paternal_id"].values[0],
        ped[ped["sample_id"] == mom]["maternal_id"].values[0],
    )

    (
        dad_dad_suff,
        dad_mom_suff,
        mom_dad_suff,
        mom_mom_suff,
    ) = [get_suffix(ped, s) for s in (dad_dad, dad_mom, mom_dad, mom_mom)]

    pgf = cyvcf2.VCF(f"{args.vcf_pref}/{dad_dad}_{dad_dad_suff}_{args.vcf_suff}")
    pgm = cyvcf2.VCF(f"{args.vcf_pref}/{dad_mom}_{dad_mom_suff}_{args.vcf_suff}")
    mgf = cyvcf2.VCF(f"{args.vcf_pref}/{mom_dad}_{mom_dad_suff}_{args.vcf_suff}")
    mgm = cyvcf2.VCF(f"{args.vcf_pref}/{mom_mom}_{mom_mom_suff}_{args.vcf_suff}")

    # read in mutations
    mutations = pd.read_csv(args.mutations, sep="\t")
    
    mutations["motif"] = mutations.apply(lambda row: query_motif_at_trid(row["trid"], pgf), axis=1)
    mutations["gp_ev"] = mutations.apply(lambda row: annotate_with_gp(row, pgf, pgm, mgf, mgm), axis=1)

    mutations.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--ped")
    p.add_argument("--focal")
    p.add_argument("--out")
    p.add_argument("--vcf_pref")
    p.add_argument("--vcf_suff")
    args = p.parse_args()
    main(args)
