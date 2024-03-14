from cyvcf2 import VCF
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import argparse
from collections import Counter
import numpy as np

from utils import extract_diffs_from_bam
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
        #DeNovoSchema.validate(mutations)
        # if we're looking at de novos, ensure that we filter
        # the DataFrame to exclude the non-denovo candidates
        mutations = mutations[mutations["denovo_coverage"] >= 1]

    # remove sex chromosomes
    mutations = mutations[~mutations["trid"].str.contains("X|Y|Un|random", regex=True)]

    # read in the VCF file
    vcf = VCF(args.vcf, gts012=True)

    print (args.mom_bam, args.dad_bam)

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
        denovo_gt = row["genotype"]

        # loop over VCF, allowing for slop to ensure we hit the STR
        for var in vcf(region):
            # ensure the VCF TRID matches the TR TRID
            var_trid = var.INFO.get("TRID")
            assert var_trid == trid

            # concatenate the alleles at this locus
            alleles = [var.REF]
            alleles.extend(var.ALT)

            # figure out the motif(s) that's repeated in the STR
            tr_motif = var.INFO.get("MOTIFS")

            # the difference between the AL fields at this locus represents the
            # expected expansion size of the STR.
            allele_lengths = var.format("AL")[0]
            allele_lengths = [len(al) for al in alleles]

            #if region != "chr1:171246527-171246564": continue

            # figure out which AL corresponds to the de novo genotype in the child.
            # to do this, we can simply calculate the lengths of the various alleles in the
            # union of the REF and ALTs. then, just index the ALs based on the genotype ID associated
            # with the denovo allele.
            denovo_al = allele_lengths[denovo_gt]
            # to figure out which AL corresponds to the "non denovo" AL, we have to figure
            # out what the non-denovo genotype ID is.

            genotypes = np.array(var.genotypes[0][:-1])
            non_denovo_gt_idx = np.where(genotypes != denovo_gt)[0]
            if non_denovo_gt_idx.shape[0] == 0: continue
            non_denovo_gt_idx = non_denovo_gt_idx[0]
            non_denovo_gt = genotypes[non_denovo_gt_idx]
            non_denovo_al = allele_lengths[non_denovo_gt]

            # if the total length of the allele is greater than the length of
            # a typical read, move on
            if (max([denovo_al, non_denovo_al]) * 2) > int(args.max_al):
                continue

            
            # calculate expected diffs between alleles and the refernece genome.
            exp_diff_denovo = denovo_al - len(var.REF)
            exp_diff_non_denovo = non_denovo_al - len(var.REF)

            # different slop for ins vs del so we capture alleels
            is_expansion = denovo_al > non_denovo_al

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
                    exp_diff_denovo,
                    exp_diff_non_denovo,
                    min_mapq=1,
                    is_expansion=is_expansion,
                )

                #print (label, diff_counts)
                row_dict = row.to_dict()

                for diff, count in diff_counts:
                    row_dict = row.to_dict()
                    row_dict.update({
                        "region": region,
                        "sample": label,
                        "motif": tr_motif,
                        "diff": diff,
                        "Read support": count,
                        "exp_allele_diff_denovo": exp_diff_denovo,
                        "exp_allele_diff_non_denovo": exp_diff_non_denovo,
                    })
                    res.append(row_dict)

    print(
        "TOTAL sites: ",
        mutations.shape[0],
        TOTAL_SITES,
    )
    for reason, count in Counter(SKIPPED).most_common():
        print(f"{reason}: {count}")
    res_df = pd.DataFrame(res)
    
    res_df[res_df["sample"] == "kid"].to_csv(args.out, index=False)


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
