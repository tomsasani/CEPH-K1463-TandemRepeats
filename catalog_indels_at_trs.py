import pysam
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import argparse
from collections import Counter
from typing import List, Tuple, Union
import math
from cyvcf2 import VCF

MATCH, INS, DEL = range(3)
OP2DIFF = {MATCH: 0, INS: 1, DEL: -1}
# from SAM spec (https://samtools.github.io/hts-specs/SAMv1.pdf)
# page 8 -- 1 if the operation consumes reference sequence
CONSUMES_REF = dict(zip(range(9), [1, 0, 1, 1, 0, 0, 0, 1, 1]))


def bp_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
    """simple utility to determine the amount of
    overlap between two interval coordinates

    Args:
        s1 (int): start of interval 1
        e1 (int): end of interval 1
        s2 (int): start of interval 2
        e2 (int): end of interval 2

    Returns:
        int: size of overlap
    """
    return max(
        max((e2 - s1), 0) - max((e2 - e1), 0) - max((s2 - s1), 0),
        0,
    )


def count_indel_in_read(
    ct: List[Tuple[int, int]],
    rs: int,
    vs: int,
    ve: int,
    slop: int = 1,
) -> int:
    """count up inserted and deleted sequence in a read
    using the pysam cigartuples object. a cigartuples object
    is a list of tuples -- each tuple stores a CIGAR operation as
    its first element and the number of bases attributed to that operation
    as its second element. exact match is (0, N), where N is the number
    of bases that match the reference, insertions are (1, N), and deletions
    are (2, N). we only count CIGAR operations that completely overlap
    the expected STR locus interval.

    Args:
        ct (Tuple[int, int]): pysam cigartuples object
        rs (int): start of read w/r/t reference
        vs (int): start of TR locus in reference
        ve (int): end of TR locus in reference

    Returns:
        int: net ins/del sequence in read w/r/t reference
    """

    # keep a running total of net ins/del operations in the read.
    cigar_op_total = 0
    # count from the first position in the read w/r/t the reference
    cur_pos = rs

    # loop over the cigartuples object, operation by operation
    for op, bp in ct:
        # check if this operation type consumes reference sequence.
        # if so, we'll keep track of the bp associated with the op.
        # if not, 
        if bool(CONSUMES_REF[op]):
            op_length = bp
        else:
            op_length = 0
        # if the operation isn't an INS, DEL, or MATCH, we can
        # just increment the running start and move on. e.g., if
        # the operation is a mismatch or a soft clip, we're not interested
        # in incrementing our net CIGAR totals by the number of bp affected
        # by that operation. however, we *do* need to increment our current
        # position.
        if op not in (INS, DEL, MATCH):
            cur_pos += op_length
            continue
        else:
            # keep track of the *relative* start and end of the operation,
            # given the number of consumable base pairs that have been
            # encountered in the iteration so far.
            op_s, op_e = cur_pos, cur_pos + max([1, op_length])
            # figure out the amount of overlap between this opeartion and
            # our TR locus. increment our counter of net CIGAR operations by this overlap.
            overlapping_bp = bp_overlap(op_s, op_e, vs - slop, ve + slop)
            if overlapping_bp > 0:
                cigar_op_total += (bp * OP2DIFF[op])
            # increment our current position counter regardless of whether the
            # operation overlaps our STR locus of interest
            cur_pos += op_length

    return cigar_op_total


def get_read_diff(
    read: pysam.AlignedSegment,
    start: int,
    end: int,
    min_mapq: int = 60,
    slop: int = 1,
) -> Union[None, int]:
    """compare a single sequencing read and to the reference. then,
    count up the net inserted/deleted sequence in the read.

    Args:
        read (pysam.AlignedSegment): pysam read (aligned segment) object.
        start (int): start position of the TR locus in the reference.
        end (int): end position of the TR locus.
        min_mapq (int, optional): minimum mapping quality for a read to be considered. Defaults to 60.
        slop (int, optional): amount of slop around the start and end of the variant reported site. Defaults to 1.

    Returns:
        Union[None, int]: either None (if the read fails basic checks) or the net ins/del in read w/r/t reference
    """
    # initial filter on mapping quality
    if read.mapping_quality < min_mapq:
        return None
    
    # get the start and end positions of the read w/r/t the reference
    qs, qe = read.reference_start, read.reference_end

    # ensure that this read completely overlaps the TR locus along with
    # slop. if the read starts after the adjusted/slopped start
    # or if it ends before the adjusted/slopped end, skip the read
    adj_start, adj_end = start - slop, end + slop
    if qs > adj_start or qe < adj_end:
        return None

    # query the CIGAR string in the read.
    diff = count_indel_in_read(
        read.cigartuples,
        qs,
        start,
        end,
        slop=slop,
    )

    return diff


def extract_diffs_from_bam(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    min_mapq: int = 60,
) -> List[Tuple[int, int]]:
    """gather information from all reads aligned to the specified region.

    Args:
        bam (pysam.AlignmentFile): pysam.AlignmentFile object
        chrom (str): chromosome we want to query
        start (int): start of TR locus
        end (int): end of TR locus
        min_mapq (int, optional): minimum MAPQ required for reads. Defaults to 60.

    Returns:
        List[Tuple[int, int]]: List of (X, Y) tuples where X is the read "diff" and Y is the
         number of reads with that "diff" 
    """
    diffs = []

    if bam is None:
        diffs.append(0)
    else:
        for read in bam.fetch(chrom, start, end):
            diff = get_read_diff(
                read,
                start,
                end,
                slop=math.floor(0.1 * (end - start)),
                min_mapq=min_mapq,
            )

            if diff is None:
                continue
            else:
                diffs.append(diff)

    # count up all recorded diffs between reads and reference allele
    diff_counts = Counter(diffs).most_common()
    return diff_counts


plt.rc("font", size=14)

RC = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_bam(fh: str):
    # read in BAM files for this sample, as well as parents if applicable
    bam = (
        pysam.AlignmentFile(
            fh,
            "rb" if fh.endswith("bam") else "rc",
            #reference_filename=ref_fh if fh.endswith("cram") else None,
        )
        if fh not in ("None", None)
        else None
    )

    return bam


def main(args):

    print (args.str_vcf)

    vcf = VCF(args.str_vcf, gts012=True)

    trids = pd.read_csv(args.trids, sep="\t")

    # read in BAM files for this sample, as well as parents if applicable
    bam = read_bam(args.bam)

    # store output
    res = []

    for i, row in tqdm.tqdm(trids.iterrows()):
        chrom, start, end = row["#chrom"], row["start"], row["end"]
        region = f"{chrom}:{start}-{end}"
        orig_trid = row["TRid"]

        for v in vcf(region):
        
            trid = v.INFO.get("TRID")
            if trid != orig_trid: continue

            start, end = int(start), int(end)
            start -= 5_000
            end += 5_000

            motif = v.INFO.get("MOTIFS")

            if len(motif) != 1: continue
            gts = v.gt_types
            if gts[0] != 0: continue

            allele_lengths = v.format("AL")[0]
            assert len(set(allele_lengths)) == 1
            al = allele_lengths[0]

            exp_allele_diff = al - (end - start)

            diff_counts = extract_diffs_from_bam(
                bam,
                chrom,
                start,
                end,
                min_mapq=60,
            )

            if len(diff_counts) == 0: continue

            evidence = {
                        f"evidence": "|".join(
                            [
                                ":".join(list(map(str, [diff, count])))
                                for diff, count in diff_counts
                            ]
                        )
                    }
        
            row_dict = {"chrom": chrom, "start": start, "end": end, "trid": trid, "exp_allele_diff": exp_allele_diff, "motif": motif}

            row_dict.update(evidence)
            res.append(row_dict)
       
    res_df = pd.DataFrame(res)

    res_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--str_vcf",
    )
    p.add_argument(
        "--bam",
        help="""Path to BAM file with aligned Elements reads for the sample of interest""",
    )
    p.add_argument("--trids")
    p.add_argument(
        "--out",
        help="""Name of output file to store read support evidence for each call.""",
    )

    args = p.parse_args()

    main(args)
