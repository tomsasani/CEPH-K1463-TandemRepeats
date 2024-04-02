from collections import Counter
from typing import Tuple, Union
import pandas as pd

import pysam
import numpy as np

def filter_mutation_dataframe(
    mutations: pd.DataFrame,
    remove_complex: bool = True,
    remove_duplicates: bool = True,
    mc_column: str = "child_MC",
    denovo_coverage_min: int = 1,
):

    mutations = mutations[mutations["denovo_coverage"] >= denovo_coverage_min]
    # only look at autosomes
    mutations = mutations[~mutations["trid"].str.contains("X|Y|Un|random", regex=True)]
    # remove pathogenics
    mutations = mutations[~mutations["trid"].str.contains("pathogenic")]
    
    if remove_complex:
        mutations = mutations[~mutations[mc_column].str.contains("_")]
    if remove_duplicates:
        trid_counts = mutations.groupby("trid").size().to_dict()
        duplicate_trids = [k for k,v in trid_counts.items() if v > 1]
        mutations = mutations[~mutations["trid"].isin(duplicate_trids)]

    return mutations


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

    # from SAM spec (https://samtools.github.io/hts-specs/SAMv1.pdf)
    # page 8 -- 1 if the operation consumes reference sequence
    CONSUMES_REF = dict(zip(range(9), [1, 0, 1, 1, 0, 0, 0, 1, 1]))

    # figure out start and end positions of the variant, adding
    # some slop if desired to either side
    var_s, var_e = vs - slop, ve + slop

    # loop over cigar operations in the read. if a cigar operation
    # is within the variant locus, add it to the running total of
    # ins/del operations in the read.
    cigar_op_total = 0
    cur_pos = rs
    for op, bp in ct:

        # check if this operation type consumes reference sequence.
        # if so, we'll keep track of the bp associated with the op.
        if bool(CONSUMES_REF[op]):
            op_length = bp
        else:
            op_length = 0

        # keep track of the *relative* start and end of the operation,
        # given the number of consumable base pairs that have been
        # encountered in the iteration so far.
        op_s, op_e = cur_pos, cur_pos + max([1, op_length])

        # if this operation overlaps our STR locus, increment our counter
        # of net CIGAR operations
        overlapping_bp = bp_overlap(op_s, op_e, var_s, var_e)

        if overlapping_bp > 0:
            if op == INS:
                # print(
                #     f"overlap: {overlapping_bp}, operation: {op}, size: {bp}, true_size: {op_length}, true_start: {op_s}, true_end: {op_e}"
                # )
                cigar_op_total += (bp * OP2DIFF[op])
            elif op in (MATCH, DEL):
                cigar_op_total += (overlapping_bp * OP2DIFF[op])
        # increment our current position counter regardless of whether the
        # operation overlaps our STR locus of interest

        cur_pos += op_length

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
    is_expansion: bool = False,
):
    diffs = []

    # if we're dealing with an expansion, we expect CIGAR
    # operations to be insertions relative to the reference.
    # thus, we can get away with requiring less slop, as the 
    # read only needs to map across the full length of the reference
    # allele.
    slop = max([abs(ref_al_diff), abs(alt_al_diff)])
    if is_expansion:
        slop = 1


    if bam is None:
        diffs.append(0)
    else:
        for read in bam.fetch(chrom, start, end):
            diff = get_read_diff(
                read,
                start,
                end,
                slop=slop,
                min_mapq=min_mapq,
            )

            if diff is None:
                continue
            else:
                diffs.append(diff)

    # count up all recorded diffs between reads and reference allele
    diff_counts = Counter(diffs).most_common()
    return diff_counts
