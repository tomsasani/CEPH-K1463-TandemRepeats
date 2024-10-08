import pandas as pd
import argparse
import cyvcf2
import numpy as np


def query_al_at_trid(row: pd.Series, vcf: cyvcf2.VCF, min_depth: int = 10,):
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    region = f"{chrom}:{start}-{end}"
    trid = row["trid"]

    allele_lengths = None
    allele_length_ranges = None
    for v in vcf(region):
        if trid != v.INFO.get("TRID"): continue
        spanning_reads = v.format("SD")
        if np.sum(spanning_reads) <= min_depth: continue
        
        allele_lengths = v.format("AL")
        allele_length_ranges = v.format("ALLR")
        
    if allele_length_ranges is None: return None
    else:
        allrs = []
        for allr in allele_length_ranges[0].split(","):
            allr_int = list(map(int, allr.split("-")))
            allrs.append(allr_int)
        return allrs

def query_motif_at_trid(row: str, vcf: cyvcf2.VCF):
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    trid = row["trid"]
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
    chrom = row["#chrom"]

    res = []

    for l, v in zip(("pgf", "pgm", "mgf", "mgm"), (pgf, pgm, mgf, mgm)):
        if v is None: continue
        # get allele lengths at this TRID in the VCF of interest
        # if we're on the X chromosome, we can require min depth of 10 in the
        # grandmothers and half that in the grandfathers
        v_allrs = query_al_at_trid(
            row,
            v,
            min_depth=5 if l in ("pgf", "mgf") and chrom in ("chrX", "chrY") else 5,
        )
        if v_allrs is None:
            res.append(f"{l}_?")
        else:
            in_allr = False
            for allr in v_allrs:
                allr_min, allr_max = allr
                if allr_min <= denovo_al <= allr_max:
                    in_allr = True
            if in_allr:
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

    dad_dad, dad_mom = None, None
    mom_dad, mom_mom = None, None

    if args.lineage in ("paternal", "both"):
        dad_dad, dad_mom = (
            ped[ped["sample_id"] == dad]["paternal_id"].values[0],
            ped[ped["sample_id"] == dad]["maternal_id"].values[0],
        )
    if args.lineage in ("maternal", "both"):
        mom_dad, mom_mom = (
            ped[ped["sample_id"] == mom]["paternal_id"].values[0],
            ped[ped["sample_id"] == mom]["maternal_id"].values[0],
        )

    (
        dad_dad_suff,
        dad_mom_suff,
        mom_dad_suff,
        mom_mom_suff,
    ) = [
        get_suffix(ped, s) if s is not None else None
        for s in (
            dad_dad,
            dad_mom,
            mom_dad,
            mom_mom,
        )
    ]

    pgf, pgm, mgf, mgm = None, None, None, None

    if args.lineage in ("paternal", "both"):
        pgf = cyvcf2.VCF(f"{args.vcf_pref}/{dad_dad}_{dad_dad_suff}_{args.vcf_suff}")
        pgm = cyvcf2.VCF(f"{args.vcf_pref}/{dad_mom}_{dad_mom_suff}_{args.vcf_suff}")
    if args.lineage in ("maternal", "both"):
        mgf = cyvcf2.VCF(f"{args.vcf_pref}/{mom_dad}_{mom_dad_suff}_{args.vcf_suff}")
        mgm = cyvcf2.VCF(f"{args.vcf_pref}/{mom_mom}_{mom_mom_suff}_{args.vcf_suff}")

    # read in mutations
    mutations = pd.read_csv(args.mutations, sep="\t")

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
    p.add_argument("--lineage")
    args = p.parse_args()
    main(args)
