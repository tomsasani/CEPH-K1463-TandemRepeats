from collections import Counter
from typing import Tuple, Union

import pysam


def bp_overlap(s1: int, e1: int, s2: int, e2: int):
    return max(
        max((e2 - s1), 0) - max((e2 - e1), 0) - max((s2 - s1), 0),
        0,
    )


def count_indel_in_read(
    ct: Tuple[int, int],
    rs: int,
    vs: int,
    ve: int,
    slop: int = 20,
):
    """count up inserted and deleted sequence in a read
    using the pysam cigartuples object. a cigartuples object
    is a list of tuples -- each tuple stores a CIGAR operation as
    its first element and the number of bases attributed to that operation
    as its second element. exact match is (0, N), where N is the number
    of bases that match the reference, insertions are (1, N), and deletions
    are (2, N). we only count CIGAR operations that completely overlap
    the expected STR locus interval (+/- the size of the STR expansion).

    Args:
        ct (Tuple[int, int]): _description_
        rs (int): _description_
        re (int): _description_
        vs (int): _description_
        ve (int): _description_

    Returns:
        _type_: _description_
    """

    MATCH, INS, DEL = range(3)

    OP2DIFF = {MATCH: 0, INS: 1, DEL: -1}

    # figure out start and end positions of the variant, adding
    # some slop if desired to either side
    var_s, var_e = vs - slop, ve + slop

    # loop over cigar operations in the read. if a cigar operation
    # is within the variant locus, add it to the running total of
    # ins/del operations in the read.
    cigar_op_total = 0
    cur_pos = rs
    for op, bp in ct:
        if op not in (MATCH, INS, DEL):
            continue
        # figure out intersection between current cigar operation span
        # and the variant locus span
        op_s, op_e = cur_pos, cur_pos + bp
        overlapping_bp = bp_overlap(op_s, op_e, var_s, var_e)
        # increment cigar ops by the total base pairs of this op
        cigar_op_total += overlapping_bp * OP2DIFF[op]
        cur_pos += bp
    return cigar_op_total


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

    diff = count_indel_in_read(ct, qs, start, end, slop=slop)

    return diff


def extract_diffs_from_bam(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    ref_al_diff: int,
    alt_al_diff: int,
    min_mapq: int = 60,
):
    diffs = []

    if bam is None:
        diffs.append(0)
    else:
        for read in bam.fetch(chrom, start, end):
            diff = get_read_diff(
                read,
                start,
                end,
                slop=max([abs(ref_al_diff), abs(alt_al_diff)]),
                min_mapq=min_mapq,
            )
            if diff is None:
                continue
            else:
                diffs.append(diff)

    # count up all recorded diffs between reads and reference allele
    diff_counts = Counter(diffs).most_common()
    return diff_counts
