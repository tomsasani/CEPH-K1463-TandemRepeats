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
from typing import Tuple, List
from Bio import Align
import matplotlib.patches as patches

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
    # find the longest continuous interval
    ilens = [i[-1] - i[0] for i in intervals]
    max_i = np.argmax(ilens)
    return intervals[max_i]


def count_indel_in_read(ct: List[Tuple[int, int]]) -> int:
    total_indel = 0
    for op, v in ct:
        if op in (1, 2):
            total_indel += v
    return total_indel


# read in de novos
dnms = pd.read_parquet("tr_validation/data/df_transmission_C5.parquet.gzip")

VCF_FH = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/f135436/2188_DM_adotto_v02.sorted.vcf.gz"
BAM_FH = "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/merged/"

FATHER, MOTHER = "2214", "2213"

vcf = VCF(VCF_FH, gts012=True)

aligner = Align.PairwiseAligner()

res = []

for i, row in dnms.iterrows():
    trid = row["trid"]
    chrom, start, end = row["chrom"], row["start"], row["end"]
    region = f"{chrom}:{start}-{end}"

    # print (row)
    # break

    # read in BAM file for this sample
    relevant_bam = BAM_FH + "{}_1_merged_sort.bam".format(row["sample_id"])
    bam = pysam.AlignmentFile(relevant_bam, "rb")

    # loop over VCF, allowing for slop to ensure we hit the STR
    for var in vcf(region):
        var_trid = var.INFO.get("TRID")
        assert var_trid == trid

        tr_motif = var.INFO.get("MOTIFS").split(",")

        # for now, ignore sites with multiple TR motifs.
        # don't think there are any in the DNM file at the moment
        if len(tr_motif) > 1:
            continue
        tr_motif = tr_motif[0]
        
        # if the repeated motif is longer than an element read, move on
        if len(tr_motif) > 150: continue
        if len(var.ALT) != 1: continue

        # perform a local alignment of the REF and ALT alleles
        # at this site to figure out how much extra sequence
        # is added in the ALT w/r/t to the REF
        # alignments = aligner.align(var.REF, var.ALT[0])

        # perform a simple regex query to figure out where in the
        # REF/ALT alleles the TR motif begins. we need to do this
        # in order to know how much padding is on either side of the
        # actual TR motif in the reported REF/ALT.
        motif_locs = count_motif_regex(var.REF, tr_motif)
        if len(motif_locs) == 0:
            print (var.format("AL"))
            print (tr_motif)
            longest_motif_interval = (0, 0)
        else:
            # get the index of the longest contiguous chunk of repeated motifs
            longest_motif_interval = find_longest_interval(motif_locs)
        # figure out the position at which the motif repeats start and end. we'll
        # use this to ensure that reads overlap the full repeat sequence.
        motif_s, motif_e = longest_motif_interval

        # figure out how long the reference and alternate alleles are
        ref_al, alt_al = var.format("AL")[0]
        # if the length of the repeat expansion is greater than the length
        # of a typical element read, move on
        if abs(alt_al - ref_al) > 100: continue

        gts = var.gt_types

        # loop over reads in the BAM for this individual
        motif_counts = []
        for ri, read in enumerate(bam.fetch(var.CHROM, var.start, var.end)):
            # NOTE: need to filter reads so that we only deal with reads
            # that completely overlap the expanded locus
            qs, qe = read.reference_start, read.reference_end
            # ensure that this read completely overlaps the expected repeat expansion.
            # NOTE: the repeat expansion begins at the start of the variant (var.start) +
            # the index in the REF allele at which the motif expansion starts (longest_motif_interval[0]).
            # the repeat expansion ends at the end of the variant (var.end) - (longest_motif_interval[1]).
            if qs > (var.start + motif_s) or qe < (var.end - motif_e):
                continue

            # query the CIGAR string in the read.
            ct = read.cigartuples
            total_indel_in_read = count_indel_in_read(ct)
            # motif_counts.append(total_indel_in_read)
            # for k, v in d.items():
            res.append(
                {
                    "region": region,
                    "motif": tr_motif,
                    "indel_in_read": total_indel_in_read,
                    "ref_al": ref_al,
                    "alt_al": alt_al,
                    "read_i": ri,
                }
            )


res_df = pd.DataFrame(res)
res_df["allele_length_diff"] = abs(res_df["ref_al"] - res_df["alt_al"])

res_df.to_csv("out.csv", index=False)

print (res_df.groupby("region").size().shape)

res_df = res_df.query("allele_length_diff > 0 and indel_in_read > 0")
f, ax = plt.subplots(figsize=(8, 6))

x = res_df["allele_length_diff"].values
y = res_df["indel_in_read"].values

jitter = np.random.normal(loc=0, scale=0.05, size=x.shape[0])

ax.scatter(x + jitter, y, alpha=0.25, c="cornflowerblue", s=100)
ax.axline((0, 0), slope=1, ls=":", c="k", lw=2, zorder=-1)
if ax.get_legend() is not None:
    ax.get_legend().remove()
sns.despine(ax=ax, top=True, right=True)
ax.set_xlabel("Sum of indel CIGAR operations in read")
ax.set_ylabel("Allele length difference b/w REF and ALT")
f.tight_layout()
f.savefig("o.png", dpi=200)

# region = res_df.iloc[0]["region"]
# res_df_sub = res_df[res_df["region"] == region]
# motif = res_df_sub["motif"].unique()[0]

# f, ax = plt.subplots()
# ax.bar(
#     res_df_sub["indel_in_read"],
#     res_df_sub["num_reads"],
#     1,
#     ec="k",
#     lw=1,
#     color="cornflowerblue",
# )
# ax.axvline(res_df_sub["allele_length_diff"].values[0])
# ax.set_xlabel("Sum of indel CIGAR operations in read")
# ax.set_ylabel("Number of reads")
# f.tight_layout()
# f.savefig("o.png", dpi=200)
