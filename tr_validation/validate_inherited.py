from cyvcf2 import VCF
import pysam
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import argparse
from validate_dnms import get_read_diff



plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def main(args):
    # read in de novo parquet file and VCF
    mutations = pd.read_csv(
        args.mutations,
        names=["trid", "maternal_id", "paternal_id", "sample_id"],
        sep="\t",
    )

    # filter to proband of interest
    mutations = mutations[mutations["sample_id"] == int(args.sample_id)]

    mom_bam, dad_bam, kid_bam = (
        pysam.AlignmentFile(args.mom_bam, "rb"),
        pysam.AlignmentFile(args.dad_bam, "rb"),
        pysam.AlignmentFile(args.kid_bam, "rb"),
    )

    # randomly sample 10_000

    #mutations = mutations[mutations["trid"] == "chr9_83786544_83786584_trsolve"]
    mutations = mutations.sample(1_000, replace=False)
    #print(mutations)


    # read in VCF for one sample so we can access AL info
    vcf = VCF(args.vcf, gts012=True)

    # store output
    res = []
    SKIPPED = []

    sizes = []

    for _, row in tqdm.tqdm(mutations.iterrows()):
        # extract info about the DNM
        trid = row["trid"]

        try:
            chrom, start, end, _ = trid.split("_")
            start, end = int(start), int(end)
            sizes.append(end - start)
        except ValueError:
            continue

        region = f"{chrom}:{start}-{end}"

        tr_motif = (None,)
        ref_al, alt_al = None, None

        for var in vcf(region):
            var_trid = var.INFO.get("TRID")


            if var_trid != trid:
                continue

            # figure out the motif that's repeated in the STR
            tr_motif = var.INFO.get("MOTIFS").split(",")
            # for now, ignore sites with multiple TR motifs.
            # don't think there are any in the DNM file at the moment
            if len(tr_motif) > 1:
                SKIPPED.append("2+ TR motifs")
                continue

            # the difference between the AL fields at this locus represents the
            # expected expansion size of the STR.
            try:
                ref_al, alt_al = var.format("AL")[0]
            except ValueError:
                SKIPPED.append("invalid AL field")
                continue
            
            if ref_al > 150: 
                SKIPPED.append("AL > 150bp")
                continue

        if ref_al is None or alt_al is None:
            continue


        # calculate expected diffs between alleles and the refernece genome
        exp_diff_alt = alt_al - len(var.REF)
        exp_diff_ref = ref_al - len(var.REF)


        # loop over reads in the BAMs
        for bam, label in zip(
            (mom_bam, dad_bam, kid_bam),
            ("mom", "dad", "kid"),
        ):
            diffs = []
            for read in bam.fetch(chrom, start, end):
                diff = get_read_diff(
                    read,
                    start,
                    end,
                    slop=max([abs(exp_diff_ref), abs(exp_diff_alt)]),
                )

                if diff is None:
                    continue
                else:
                    diffs.append(diff)

            # count up all recorded diffs between reads and reference allele
            diff_counts = Counter(diffs).most_common()

            for diff, count in diff_counts:
                d = {
                    "trid": trid,
                    "diff": diff,
                    "Read support": count,
                    "motif": tr_motif[0],
                    "sample": label,
                    "ref_al": ref_al,
                    "alt_al": alt_al,
                    "exp_allele_diff_alt": exp_diff_alt,
                    "exp_allele_diff_ref": exp_diff_ref,
                }
                res.append(d)
    
    for reason, count in Counter(SKIPPED).most_common():
        print(f"{reason}: {count}")

    print (min(sizes), max(sizes))

    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        help="""Path to TSV file containing mutations""",
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
        "--mom_bam",
        help="""Path to BAM file with aligned Elements reads for the sample of interest""",
    )
    p.add_argument(
        "--dad_bam",
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

    args = p.parse_args()
    main(args)
