from cyvcf2 import VCF
import pysam
from collections import Counter, defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gzip
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import re
from typing import Tuple, List, Union
from Bio import Align
import matplotlib.patches as patches
import cyvcf2
import argparse
import doctest

plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_exclude(fh: str) -> IntervalTree:
    """
    Read in a BED file containing genomic regions from which we want
    to exclude potential variants. 

    Args:
        fh (str): Path to filename containing regions. Must be BED-formatted. Can be \
            uncompressed or gzipped.

    Returns:
        tree (Dict[IntervalTree]): Dictionary of IntervalTree objects, with one key per \
            chromosome. Each IntervalTree containing the BED regions from `fh`, on which we can \
            quickly perform binary searches later.
    """

    tree = defaultdict(IntervalTree)
    is_zipped = fh.endswith(".gz")

    print("BUILDING EXCLUDE TREE")

    with gzip.open(fh, "rt") if is_zipped else open(fh, "rt") as infh:
        csvf = csv.reader(infh, delimiter="\t")
        for l in tqdm.tqdm(csvf):
            if l[0].startswith("#") or l[0] == "chrom":
                continue
            chrom, start, end = l
            interval = Interval(int(start), int(end))
            tree[chrom].insert_interval(interval)

    return tree


def count_indel_in_read(ct: List[Tuple[int, int]]) -> int:
    """count up inserted and deleted sequence in a read
    using the pysam cigartuples object. a cigartuples object
    is a list of tuples -- each tuple stores a CIGAR operation as
    its first element and the number of bases attributed to that operation
    as its second element. exact match is (0, N), where N is the number
    of bases that match the reference, insertions are (1, N), and deletions
    are (2, N).

    Args:
        ct (List[Tuple[int, int]]): cigartuples object from pysam aligned segment

    Returns:
        int: total in/del sequence (negative if net deleted, positive if net inserted)

    >>> count_indel_in_read([(0, 50), (1, 4), (2, 3), (1, 1), (0, 90), (4, 2)])
    2
    >>> count_indel_in_read([(0, 50), (1, 4), (2, 3), (1, 7), (0, 10)])
    8
    >>> count_indel_in_read([(0, 150)])
    0
    """
    total_indel = 0
    for op, v in ct:
        if op == 1:
            total_indel += v
        elif op == 2:
            total_indel -= v
    return total_indel


def reformat_cigartuples(
    ct: Tuple[int, int],
    rs: int,
    re: int,
    vs: int,
    ve: int,
    slop: int = 20,
):
    """_summary_

    Args:
        ct (Tuple[int, int]): _description_
        rs (int): _description_
        re (int): _description_
        vs (int): _description_
        ve (int): _description_

    Returns:
        _type_: _description_
    """
    # expand the CIGAR operations
    cigar_arr = np.zeros(re - rs)
    cur_s = 0
    for k, v in ct:
        cigar_arr[cur_s : cur_s + v] = k
        cur_s += v

    # figure out the interval of overlap between the
    # read and the variant
    assert rs < re
    assert vs < ve

    read_arr = set(range(rs, re))
    var_arr = set(range(vs - slop, ve + slop))
    overlap_arr = list(read_arr.intersection(var_arr))

    indel_ops = 0
    for i in np.arange(cigar_arr.shape[0]):
        cur_pos = i + rs
        if cur_pos not in overlap_arr:
            continue
        op = cigar_arr[i]
        if op == 1:
            indel_ops += 1
        elif op == 2:
            indel_ops -= 1
        # elif op in (4, 5):
        #     indel_ops += 1
        else: continue

    return indel_ops


def get_read_diff(
    read: pysam.AlignedSegment,
    start: int,
    end: int,
    min_mapq: int = 20,
    slop: int = 0,
) -> Union[None, int]:
    """compare a single sequencing read and to the reference. then,
    count up the net inserted/deleted sequence in the read.

    Args:
        read (pysam.AlignedSegment): pysam read (aligned segment) object.
        start (int): start position of the variant as reported in the VCF.
        end (int): end position of the variant as reported in the VCF.
        min_mapq (int, optional): minimum mapping quality for a read to be considered. Defaults to 60.
        slop (int, optional): amount of slop around the start and end of the variant reported site. Defaults to 0.

    Returns:
        Union[None, int]: _description_
    """
    # initial filter on mapping quality
    if read.mapping_quality < min_mapq:
        return None
    
    qs, qe = read.reference_start, read.reference_end
    # ensure that this read completely overlaps the expected repeat expansion.
    # NOTE: we are *very* stringent when considering reads, and only include
    # reads that align completely across the reported variant location +/-
    # the expected size of the expanded allele. for example, if the location 
    # of the reference STR is reported in the VCF as chr1:50-60 and we expect
    # the expansion to be 10bp in size, we only consider reads that align
    # completely across the interval chr1:40-70.
    if qs > (start - slop) or qe < (end + slop):
        return None

    # query the CIGAR string in the read.
    ct = read.cigartuples
    # we can simply count up the indel content of the read (w/r/t what's in the
    # reference genome) as a proxy for expanded/deleted sequence.
    diff = reformat_cigartuples(ct, qs, qe, start, end, slop=slop)
    
    return diff



def main(args):
    # read in de novo parquet file and VCF
    dnms = pd.read_csv(args.dnms) if args.dnms.endswith("csv") else pd.read_parquet(args.dnms) 

    # dnms = dnms.sample(10_000, replace=False)

    # count DNMs observed in
    num_obs = (
        dnms.groupby("trid")
        .size()
        .sort_values()
        .reset_index()
        .rename(columns={0: "num_obs"})
    )
    dnms = dnms.merge(num_obs, on="trid")
    # dnms = dnms.query("num_obs == 1")

    vcf = VCF(args.vcf, gts012=True)

    # read in BAM files for this sample
    kid_bam = pysam.AlignmentFile(
        args.kid_bam,
        "rb" if args.kid_bam.endswith("bam") else "rc",
    )
    mom_bam = (
        pysam.AlignmentFile(
            args.mom_bam,
            "rb" if args.mom_bam.endswith("bam") else "rc",
        )
        if args.mom_bam != "None"
        else None
    )
    dad_bam = (
        pysam.AlignmentFile(
            args.dad_bam,
            "rb" if args.dad_bam.endswith("bam") else "rc",
        )
        if args.mom_bam != "None"
        else None
    )
    
    if "sample_id" in dnms.columns:
        dnms = dnms[
            (dnms["sample_id"].astype(str) == args.sample_id) & (dnms["genotype"] != 0)
        ]

    # store output
    res = []

    SKIPPED, TOTAL_SITES = [], 0
    for _, row in tqdm.tqdm(dnms.iterrows()):
        # extract info about the DNM
        trid = row["trid"]
        # chrom, start, end = row["chrom"], row["start"], row["end"]
        try:
            chrom, start, end, method = trid.split("_")
        except ValueError:
            continue

        region = f"{chrom}:{start}-{end}"

        # loop over VCF, allowing for slop to ensure we hit the STR
        for var in vcf(region):
            TOTAL_SITES += 1
            var_trid = var.INFO.get("TRID")

            assert var_trid == trid

            # concatenate the alleles at this locus
            alleles = [var.REF]
            alleles.extend(var.ALT)

            # figure out the motif that's repeated in the STR
            tr_motif = var.INFO.get("MOTIFS").split(",")
            # for now, ignore sites with multiple TR motifs.
            # don't think there are any in the DNM file at the moment
            if len(tr_motif) > 1:
                SKIPPED.append("2+ TR motifs")
                continue
            tr_motif = tr_motif[0]

            # figure out the genotype of this sample
            gt = var.genotypes
            if gt[0] == 0:
                SKIPPED.append("Sample doesn't have an ALT allele")
                continue

            # the difference between the AL fields at this locus represents the
            # expected expansion size of the STR.
            try:
                ref_al, alt_al = var.format("AL")[0]
            except ValueError:
                SKIPPED.append("invalid AL field")
                continue

            # if the total length of the allele is greater than the length of
            # a typical element read, move on
            if alt_al > 150:
                SKIPPED.append("ALT allele length > 150 bp")
                continue

            if alt_al - ref_al == 0:
                SKIPPED.append("Allele lengths are identical")
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
                diffs = []

                if bam is None:
                    diffs.append(0)
                else:
                    for read in bam.fetch(chrom, int(start), int(end)):
                        diff = get_read_diff(
                            read,
                            int(start),
                            int(end),
                            slop=max([abs(exp_diff_ref), abs(exp_diff_alt)]),
                        )
                        if diff is None:
                            continue
                        else:
                            diffs.append(diff)
                # count up all recorded diffs between reads and reference allele
                diff_counts = Counter(diffs).most_common()

                # total_reads = sum([v for k,v in diff_counts])
                # if total_reads < 10: continue

                for diff, count in diff_counts:
                    d = {
                        "region": region,
                        "sample": label,
                        "motif": tr_motif,
                        "diff": diff,
                        "Read support": count,
                        "exp_allele_diff_alt": exp_diff_alt,
                        "exp_allele_diff_ref": exp_diff_ref,
                        "n_with_denovo_allele_strict": row[
                            "n_with_denovo_allele_strict"
                        ] if "n_with_denovo_allele_strict" in row else 0,
                    }
                    res.append(d)

    # print(
    #     "TOTAL sites: ",
    #     dnms[dnms["sample_id"] == args.sample_id].shape[0],
    #     TOTAL_SITES,
    # )
    for reason, count in Counter(SKIPPED).most_common():
        print(f"{reason}: {count}")
    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--dnms",
        help="""Path to parquet file containing candidate DNMs""",
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
        "--out",
        help="""Name of output file to store read support evidence for each call.""",
    )
    p.add_argument(
        "--sample_id",
        help="""Name of sample ID you wish to query.""",
    )

    args = p.parse_args()
    doctest.testmod()
    main(args)
