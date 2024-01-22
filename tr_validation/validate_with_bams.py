from cyvcf2 import VCF
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import argparse
from collections import Counter

from utils import extract_diffs_from_bam
from schema import DeNovoSchema, InheritedSchema

plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_bam(fh: str):
    # read in BAM files for this sample, as well as parents if applicable
    bam = (
        pysam.AlignmentFile(fh, "rb" if fh.endswith("bam") else "rc")
        if fh != "None"
        else None
    )

    return bam


def main(args):
    # read in mutation file and validate columns in dataframe using pandera schema
    if args.variant_type == "inherited":
        mutations = pd.read_csv(
            args.mutations,
            sep="\t",
            names=["trid", "maternal_id", "paternal_id", "sample_id"],
            dtype=str,
        ).sample(1_000)
        InheritedSchema.validate(mutations)
    elif args.variant_type == "dnm":
        mutations = pd.read_parquet(args.mutations)
        DeNovoSchema.validate(mutations)
        # if we're looking at de novos, ensure that we filter
        # the DataFrame to exclude the non-denovo candidates
        mutations = mutations[mutations["denovo_status"].str.contains("Y")]

    # filter to the samples of interest
    mutations = mutations[mutations["sample_id"].astype(str) == args.sample_id]

    # remove sex chromosomes, and add chromosome column if not already provided
    if "chrom" in mutations.columns:
        mutations = mutations[mutations["chrom"] != "chrX"]
    else:
        mutations["chrom"] = mutations["trid"].apply(lambda t: t.split("_")[0])
        mutations = mutations[mutations["chrom"] != "chrX"]

    # read in the VCF file
    vcf = VCF(args.vcf, gts012=True)

    # read in BAM files for this sample, as well as parents if applicable
    kid_bam, mom_bam, dad_bam = (
        read_bam(fh) for fh in (args.kid_bam, args.mom_bam, args.dad_bam)
    )

    # store output
    res = []

    SKIPPED, TOTAL_SITES = [], 0

    # loop over mutations in the dataframe
    for _, row in tqdm.tqdm(mutations.iterrows()):
        # extract chrom, start and end
        trid = row["trid"]
        try:
            chrom, start, end, _ = trid.split("_")
            start, end = int(start) - 1, int(end)
        except ValueError:
            continue

        region = f"{chrom}:{start}-{end}"

        # loop over VCF, allowing for slop to ensure we hit the STR
        for var in vcf(region):
            # ensure the VCF TRID matches the TR TRID
            var_trid = var.INFO.get("TRID")
            assert var_trid == trid

            # concatenate the alleles at this locus
            alleles = [var.REF]
            alleles.extend(var.ALT)

            # figure out the motif(s) that's repeated in the STR
            tr_motif = var.INFO.get("MOTIFS").split(",")

            # loop over motifs
            for motif_i, motif in enumerate(tr_motif):
                # the difference between the AL fields at this locus represents the
                # expected expansion size of the STR.
                ref_al, alt_al = var.format("AL")[0]

                # if the total length of the allele is greater than the length of
                # a typical read, move on
                if alt_al > args.max_al:
                    continue

                # calculate expected diffs between alleles and the refernece genome
                exp_diff_alt = alt_al - len(var.REF)
                exp_diff_ref = ref_al - len(var.REF)

                # loop over reads in the BAM for this individual
                # loop over reads in the BAMs
                for bam, label in zip(
                    (mom_bam, dad_bam, kid_bam),
                    ("mom", "dad", "kid"),
                ):
                    diff_counts = extract_diffs_from_bam(
                        bam,
                        chrom,
                        start,
                        end,
                        exp_diff_ref,
                        exp_diff_alt,
                    )

                    for diff, count in diff_counts:
                        d = {
                            "region": region,
                            "sample": label,
                            "motif": motif,
                            "motif_number": motif_i + 1,
                            "diff": diff,
                            "Read support": count,
                            "exp_allele_diff_alt": exp_diff_alt,
                            "exp_allele_diff_ref": exp_diff_ref,
                            "n_with_denovo_allele_strict": row[
                                "n_with_denovo_allele_strict"
                            ],
                        }
                        res.append(d)

    print(
        "TOTAL sites: ",
        mutations.shape[0],
        TOTAL_SITES,
    )
    for reason, count in Counter(SKIPPED).most_common():
        print(f"{reason}: {count}")
    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        help="""Path to parquet/tsv file containing candidate mutations""",
    )
    p.add_argument(
        "--vcf",
        help="""Path to TRGT VCF file for the sample of interest""",
    )
    p.add_argument(
        "--kid_bam",
        help="""Path to BAM file with aligned Elements reads for the sample of interest""",
    )
    p.add_argument(
        "--out",
        help="""Name of output file to store read support evidence for each call.""",
    )
    p.add_argument(
        "--sample_id",
        help="""Name of sample ID you wish to query.""",
    )
    p.add_argument(
        "--variant_type",
        type=str,
        help="Options are dnm and inherited",
    )
    p.add_argument(
        "-mom_bam",
        help="""Path to BAM file with aligned Elements reads for the sample of interest""",
        default=None,
    )
    p.add_argument(
        "-dad_bam",
        help="""Path to BAM file with aligned Elements reads for the sample of interest""",
        default=None,
    )
    p.add_argument(
        "-max_al",
        type=int,
        default=120,
        help="""Maximum allele length to consider.""",
    )

    args = p.parse_args()

    # lp = LineProfiler()
    # lp.add_function(reformat_cigartuples)

    # lp_wrapper = lp(main)
    # lp_wrapper(args)
    # lp.print_stats()
    main(args)
