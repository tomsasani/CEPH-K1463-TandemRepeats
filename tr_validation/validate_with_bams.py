from cyvcf2 import VCF
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import argparse
from collections import Counter
import numpy as np
import csv
from utils import extract_diffs_from_bam, filter_mutation_dataframe
from schema import DeNovoSchema, InheritedSchema

plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_bam(fh: str):
    # read in BAM files for this sample, as well as parents if applicable
    bam = (
        pysam.AlignmentFile(fh, "rb" if fh.endswith("bam") else "rc")
        if fh not in ("None", None)
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
        mutations = pd.read_csv(args.mutations, sep="\t")

    # read in the VCF file
    vcf = VCF(args.vcf, gts012=True)

    # read in BAM files for this sample, as well as parents if applicable
    kid_bam, mom_bam, dad_bam = (
        read_bam(fh) for fh in (args.kid_bam, args.mom_bam, args.dad_bam)
    )

    # store output
    res = []

    # loop over mutations in the dataframe
    with open(args.mutations, "r") as infh:
        csvf = csv.reader(infh, delimiter="\t")
        header = None
        for i,l in tqdm.tqdm(enumerate(csvf)):
            if i == 0: 
                header = l
                continue
            # convert the pd.Series into a dict we can update
            # with addtl info later
            row_dict = dict(zip(header, l))

            # extract chrom, start and end
            trid = row_dict["trid"]
            chrom, start, end, _ = trid.split("_")
            start, end = int(start) - 1, int(end)
            region = f"{chrom}:{start}-{end}"

            # for sites at which we've detected a de novo,
            # figure out which allele lengths correspond to the
            # de novo and non de novo alleles.
            denovo_idx = int(row_dict["index"])
            assert denovo_idx in (0, 1)

            allele_lengths = list(map(int, row_dict["child_AL"].split(",")))
            denovo_al = allele_lengths[denovo_idx]
            non_denovo_al = allele_lengths[1 - denovo_idx]

            # loop over VCF
            for var in vcf(region):
                # ensure the VCF TRID matches the TR TRID
                var_trid = var.INFO.get("TRID")
                assert var_trid == trid

                # concatenate the alleles at this locus
                alleles = [var.REF]
                alleles.extend(var.ALT)

                # figure out the motif(s) that's repeated in the STR
                tr_motif = var.INFO.get("MOTIFS")

                # calculate expected diffs between alleles and the reference genome.
                exp_diff_denovo = denovo_al - len(var.REF)
                exp_diff_non_denovo = non_denovo_al - len(var.REF)

                # if the total length of either allele is greater than the length of
                # a typical read, move on
                if max([denovo_al, non_denovo_al]) > int(args.max_al):
                    continue

                row_dict.update(
                    {
                        "region": region,
                        "motif": tr_motif,
                        "exp_allele_diff_denovo": exp_diff_denovo,
                        "exp_allele_diff_non_denovo": exp_diff_non_denovo,
                    }
                )

                # loop over reads in the BAM for this individual
                bam_evidence = {
                    "mom_evidence": None,
                    "dad_evidence": None,
                    "kid_evidence": None,
                }

                for bam, label in zip(
                    (mom_bam, dad_bam, kid_bam),
                    ("mom", "dad", "kid"),
                ):
                    diff_counts = extract_diffs_from_bam(
                        bam,
                        chrom,
                        start,
                        end,
                        min_mapq=20,
                    )

                    # calculate the total number of queryable reads
                    # for this individual. if that's < 10, move on.
                    total_depth = sum([v for k, v in diff_counts])
                    if total_depth < 10: 
                        continue
                    else:
                        # otherwise, summarize the orthogonal evidence
                        evidence = {
                            f"{label}_evidence": "|".join(
                                [
                                    ":".join(list(map(str, [diff, count])))
                                    for diff, count in diff_counts
                                ]
                            )
                        }
                        bam_evidence.update(evidence)

                # if any members of the trio had fewer than 10 "good" reads,
                # skip this site
                if any([v is None for k, v in bam_evidence.items()]):
                    continue
                else:
                    row_dict.update(bam_evidence)
                    res.append(row_dict)

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

    main(args)
