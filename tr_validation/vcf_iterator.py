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

plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
BAM_FH = "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/complete/"


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


def revcomp(s: str):
    return "".join([RC[k] for k in s[::-1]])


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """take a list of (a, b) interval locations and merge into the
    simplest representation of those intervals.

    Args:
        intervals (List[Tuple[int, int]]): _description_

    Returns:
        List[Tuple[int, int]]: _description_
    """
    merged = []
    merged.append(intervals[0])
    for a, b in intervals[1:]:
        prev_a, prev_b = merged[-1]
        # if the right side of the previous interval
        # is the same as the left side of the current interval...
        if prev_b == a:
            # ...create a new merged interval
            new = (prev_a, b)
            merged[-1] = new
        else:
            merged.append((a, b))
    return merged


def count_motif_regex(read: str, motif: str) -> Tuple[int, int, int]:
    """count TR motifs in a string using basic regex.

    Args:
        read (str): target sequence
        motif (str): query sequence

    Returns:
        Tuple[int, int, int]: tuple of three integers
    """
    locations = []
    cur_start, cur_end = None, None
    for m in re.finditer(motif, read):
        ms, me = m.start(), m.end()
        if cur_start is None:
            locations.append((ms, me))
        else:
            if ms == cur_end:
                locations.append((ms, me))
            else:
                pass
        cur_start = ms
        cur_end = me

    # merge the intervals
    if len(locations) == 0:
        return []
    else:
        merged_intervals = merge_intervals(locations)
        return merged_intervals


def find_longest_interval(intervals: List[Tuple[int, int]]) -> Tuple[int, int]:
    ilens = [i[-1] - i[0] for i in intervals]
    max_i = np.argmax(ilens)
    return intervals[max_i]


def count_indel_in_read(ct: List[Tuple[int, int]]) -> int:
    """count up inserted and deleted sequence in a read
    using the pysam cigartuples object.

    Args:
        ct (List[Tuple[int, int]]): cigartuples object from pysam aligned segment

    Returns:
        int: total in/del sequence (negative if net deleted, positive if net inserted)
    """
    total_indel = 0
    for op, v in ct:
        if op == 1:
            total_indel += v
        elif op == 2:
            total_indel -= v
    return total_indel


def get_read_diff(
    read: pysam.AlignedSegment,
    ref_al: int,
    alt_al: int,
    var: cyvcf2.Variant,
    min_mapq: int = 10,
) -> Union[None, int]:
    """compare a single sequencing read and compare it to the reference.
    count up the net inserted/deleted sequence in the read, and figure out 
    if that amount of indel sequence matches what we'd expect for a read 
    supporting the reference allele at this locus. if not, record the amount
    of indel sequence in the read.

    Args:
        read (pysam.AlignedSegment): _description_
        ref_al (int): _description_
        alt_al (int): _description_
        var (cyvcf2.Variant): _description_
        min_mapq (int, optional): _description_. Defaults to 10.

    Returns:
        Union[None, int]: _description_
    """
    
    # initial filter on mapping quality
    if read.mapping_quality < min_mapq:
        return None
    
    # NOTE: need to filter reads so that we only deal with reads
    # that completely overlap the expanded locus
    qs, qe = read.reference_start, read.reference_end
    # ensure that this read completely overlaps the expected repeat expansion.
    # NOTE: this is hacky and doesn't account for flanking sequence
    if qs > var.start or qe < var.end:
        return None

    # we can simply count up the indel content of the read (w/r/t what's in the
    # reference genome) as a proxy for expanded/deleted sequence.

    # query the CIGAR string in the read.
    ct = read.cigartuples
    # figure out the difference in allele length between the read and the reference
    diff = count_indel_in_read(ct)

    # if this difference matches what we'd expect for a read supporting the 
    # reference allele, move on.
    # NOTE: should probably count this kind of read support in some way.
    if diff == ref_al - len(var.REF):
        return None
    
    return diff


def main(args):
    # read in de novo parquet file and VCF
    dnms = pd.read_parquet(args.dnms)
    vcf = VCF(args.vcf, gts012=True)

    # store output
    res = []

    SKIPPED = 0
    for _, row in dnms.iterrows():
        trid = row["trid"]
        chrom, start, end = row["chrom"], row["start"], row["end"]
        region = f"{chrom}:{start}-{end}"

        if row["genotype"] == 0:
            continue
        if row["sample_id"] != args.sample_id:
            continue

        # read in BAM file for this sample
        relevant_bam = BAM_FH + "{}_1_merged_sort.bam".format(args.sample_id)
        bam = pysam.AlignmentFile(relevant_bam, "rb")

        # loop over VCF, allowing for slop to ensure we hit the STR
        for var in vcf(region):
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
                SKIPPED += 1
                continue
            tr_motif = tr_motif[0]

            # figure out the genotype of this sample
            gt = var.genotypes
            # the sample might be genotype 0/1 or 1/2
            ref_idx, alt_idx = gt[0][:2]
            ref_allele, alt_allele = alleles[ref_idx], alleles[alt_idx]

            # the difference between the AL fields at this locus represents the
            # expected expansion size of the STR.
            try:
                ref_al, alt_al = var.format("AL")[0]
            except ValueError:
                SKIPPED += 1
                continue

            # if the length of the repeat expansion is greater than the length
            # of a typical element read, move on
            allele_length_diff = alt_al - ref_al
            if abs(allele_length_diff) > 100 or allele_length_diff == 0:
                SKIPPED += 1
                continue

            # loop over reads in the BAM for this individual
            diffs = []
            for ri, read in enumerate(bam.fetch(var.CHROM, var.start, var.end)):
                diff = get_read_diff(read, ref_al, alt_al, var)
                if diff is None: continue
                else: diffs.append(diff)

            # count up all recorded diffs between reads and reference allele
            diff_counts = Counter(diffs).most_common()

            for diff, count in diff_counts[:2]:
                res.append(
                    {
                        "region": region,
                        "motif": tr_motif,
                        "diff": diff,
                        "count": count,
                        "exp_allele_diff": alt_al - len(var.REF),
                        "n_with_denovo_allele_strict": row[
                            "n_with_denovo_allele_strict"
                        ],
                    }
                )

    print("SKIPPED", SKIPPED)

    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dnms")
    p.add_argument("--vcf")
    p.add_argument("--out")
    p.add_argument("--sample_id")
    args = p.parse_args()
    main(args)

