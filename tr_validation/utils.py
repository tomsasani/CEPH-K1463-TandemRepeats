from collections import Counter
from typing import Tuple, Union, List
import pandas as pd

import pysam
import numpy as np


MATCH, INS, DEL = range(3)
OP2DIFF = {MATCH: 0, INS: 1, DEL: -1}
# from SAM spec (https://samtools.github.io/hts-specs/SAMv1.pdf)
# page 8 -- 1 if the operation consumes reference sequence
CONSUMES_REF = dict(zip(range(9), [1, 0, 1, 1, 0, 0, 0, 1, 1]))

def al_in_parents(row: pd.Series):
    child_als = row["child_AL"].split(",")
    # denovo idx
    denovo_idx = int(row["index"])

    mom_als, dad_als = row["mother_AL"].split(","), row["father_AL"].split(",")

    # if we're on a male X, there's only one haplotype so the "non-denovo idx"
    # is not applicable. we only have to ask if the denovo AL is among the mom's ALs
    if row["suffix"].startswith("S") and row["trid"].split("_")[0] in ("chrX", "chrY"):
        denovo_al = child_als[denovo_idx]
        if denovo_al in mom_als: return True
        else: return False

    else:
        non_denovo_idx = 1 - denovo_idx
        try:
            denovo_al, non_denovo_al = child_als[denovo_idx], child_als[non_denovo_idx]
        except IndexError: print (child_als, row["trid"], row["suffix"])
        if (denovo_al in mom_als) and (non_denovo_al in dad_als):
            return True
        elif (denovo_al in dad_als) and (non_denovo_al in mom_als):
            return True
        else: return False


def add_likely_de_novo_size(row: pd.Series):

    child_als = list(map(int, row["child_AL"].split(",")))
    # denovo idx
    denovo_idx = int(row["index"])
    denovo_al = child_als[denovo_idx]

    if not row["trid"].startswith("chrX"):
        non_denovo_idx = 1 - denovo_idx
        non_denovo_al = child_als[non_denovo_idx]

    # figure out the inferred phase of the site, if not unknown
    phase = row["phase_summary"].split(":")[0]

    if phase == "unknown":
        return np.nan

    else:
        if phase == "dad":
            parent_als = list(map(int, row["father_AL"].split(",")))
        elif phase == "mom":
            parent_als = list(map(int, row["mother_AL"].split(",")))

        abs_denovo_diffs = np.array([abs(denovo_al - al) for al in parent_als])
        denovo_diffs = np.array([denovo_al - al for al in parent_als])
        min_diff_idx = np.argmin(abs_denovo_diffs)
        # likely original AL
        return denovo_diffs[min_diff_idx]

def filter_mutation_dataframe(
    mutations: pd.DataFrame,
    remove_complex: bool = False,
    remove_duplicates: bool = False,
    remove_gp_ev: bool = False,
    remove_inherited: bool = False,
    parental_overlap_frac_max: float = 1,
    denovo_coverage_min: int = 1,
    child_ratio_min: float = 0.,
    depth_min: int = 10,
):

    if denovo_coverage_min > 0:
        mutations = mutations[
            (mutations["denovo_coverage"] >= denovo_coverage_min)
           & (mutations["allele_ratio"] > 0)
           & (mutations["child_ratio"] > 0)
        ]
    else:
        mutations = mutations[
            (mutations["denovo_coverage"] >= denovo_coverage_min)
        ]

    mutations = mutations[mutations["child_ratio"] >= child_ratio_min]
    # only look at autosomes + X
    mutations = mutations[~mutations["trid"].str.contains("Un|random", regex=True)]
    # remove pathogenics
    mutations = mutations[~mutations["trid"].str.contains("pathogenic")]
    # remove non-TRF/TRSOLVE
    mutations = mutations[mutations["trid"].str.contains("trsolve|TRF", regex=True)]


    # remove sites where de novo AL is observed in a parent if we're interested
    # in filtering DNMs
    if denovo_coverage_min > 0:
        mutations["denovo_al_in_parents"] = mutations.apply(lambda row: al_in_parents(row), axis=1)
        if remove_inherited:
            mutations = mutations[~mutations["denovo_al_in_parents"]]

        # no parental dropout, only when we filter DNMs
        mutations = mutations[(mutations["father_dropout"] == "N") & (mutations["mother_dropout"] == "N")]

    # annotate with depth and filter on depth
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(lambda a: sum(list(map(int, a.split(",")))))
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(lambda a: sum(list(map(int, a.split(",")))))
    mutations = mutations.query(f"child_coverage >= {depth_min} and mom_dp >= {depth_min} and dad_dp >= {depth_min}")

    if remove_complex:
        mutations = mutations[~mutations["child_MC"].str.contains("_")]
    if remove_duplicates:
        dup_cols = ["trid"]
        if "sample_id" in mutations.columns:
            dup_cols.append("sample_id")
        mutations = mutations.drop_duplicates(dup_cols)

    mutations["dad_ev"] = mutations["father_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))
    mutations["mom_ev"] = mutations["mother_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))

    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(lambda o: sum(list(map(int, o.split(",")))))
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(lambda o: sum(list(map(int, o.split(",")))))

    # calculate fraction of reads supporting the denovo in the parents
    mutations["parental_overlap_coverage_total"] = mutations["dad_ev"] + mutations["mom_ev"]
    mutations["parental_coverage_total"] = mutations["mom_dp"] + mutations["dad_dp"]
    mutations["parental_overlap_coverage_frac"] = mutations["parental_overlap_coverage_total"] / mutations["parental_coverage_total"]
    mutations = mutations.query(f"parental_overlap_coverage_frac <= {parental_overlap_frac_max}")

    # if we want to remove grandparental evidence
    if remove_gp_ev and "gp_ev" in mutations.columns:
        mutations = mutations[~mutations["gp_ev"].str.contains("Y")]

    # mutations["denovo_size"] = mutations.apply(lambda row: add_likely_de_novo_size(row), axis=1)

    return mutations


# @numba.njit
def bp_overlap(s1: int, e1: int, s2: int, e2: int):
    return max(
        max((e2 - s1), 0) - max((e2 - e1), 0) - max((s2 - s1), 0),
        0,
    )


def count_indel_in_read(
    ct: List[Tuple[int, int]], rs: int, vs: int, ve: int, slop: int = 1
) -> int:
    """count up inserted and deleted sequence in a read
    using the pysam cigartuples object. a cigartuples object
    is a list of tuples -- each tuple stores a CIGAR operation as
    its first element and the number of bases attributed to that operation
    as its second element. exact match is (0, N), where N is the number
    of bases that match the reference, insertions are (1, N), and deletions
    are (2, N). we only count CIGAR operations that completely overlap
    the expected STR locus interval (+/- the size of the STR expansion).

    Args:
        ct (Tuple[int, int]): cigartuples object
        rs (int): start of read w/r/t reference
        re (int): end of read w/r/t reference
        vs (int): start of TR locus in reference
        ve (int): end of TR locus in reference

    Returns:
        int: net ins/del sequence in read w/r/t reference
    """

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
        # if the operation isn't an INS, DEL, or MATCH, we can
        # just increment the running start and move on
        if op not in (INS, DEL, MATCH):
            cur_pos += op_length
            continue
        else:
            # keep track of the *relative* start and end of the operation,
            # given the number of consumable base pairs that have been
            # encountered in the iteration so far.
            op_s, op_e = cur_pos, cur_pos + max([1, op_length])
            # if this operation overlaps our STR locus, increment our counter
            # of net CIGAR operations.
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
    min_mapq: int = 20,
    slop: int = 1,
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
        Union[None, int]: either None (if the read fails basic checks) or the net ins/del in read w/r/t reference
    """

    # initial filter on mapping quality
    if read.mapping_quality < min_mapq:
        return None

    qs, qe = read.reference_start, read.reference_end

    # ensure that this read completely overlaps the TR locus
    if qs > (start - slop) or qe < (end + slop):
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
                slop=0.5 * (end - start),
                min_mapq=min_mapq,
            )

            if diff is None:
                continue
            else:
                diffs.append(diff)

    # count up all recorded diffs between reads and reference allele
    diff_counts = Counter(diffs).most_common()
    return diff_counts
