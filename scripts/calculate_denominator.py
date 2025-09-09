import tqdm 
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import utils
import numpy as np
from snakemake.script import snakemake
from collections import defaultdict
from bx.intervals.intersection import Interval, IntervalTree
import csv

def annotate_with_censat(row: pd.Series, censat):
    overlaps = censat[row["#chrom"]].find(row["start"], row["end"])
    if len(overlaps) == 0: 
        return "no"
    else:
        return overlaps[0].value["kind"].split("_")[0]

censat = defaultdict(IntervalTree)
with open(snakemake.input.censat, "r") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l[:3]
        censat[chrom].insert_interval(Interval(int(start), int(end), value={"kind": l[3]}))


# read in raw mutations to define a denominator
raw_loci = pd.read_csv(snakemake.input.loci, sep="\t")
raw_loci["overlaps_censat"] = raw_loci.apply(lambda row: annotate_with_censat(row, censat), axis=1)

# add column with binned reference length
raw_loci["reflen"] = raw_loci["end"] - raw_loci["start"]
min_reflen, max_reflen = raw_loci["reflen"].min(), raw_loci["reflen"].max()

reflen_bins = np.arange(min_reflen, max_reflen, 10)
bins, labels = pd.cut(raw_loci["reflen"].values, bins=reflen_bins, retbins=True)

raw_loci["reflen_bin"] = bins

# calculate total number of base pairs contained within each locus
raw_loci["aff_bp"] = raw_loci["end"] - raw_loci["start"]
total_bp = raw_loci.groupby("simple_motif_size").agg(total_aff_bp = ("aff_bp", "sum")).reset_index()

GROUP_COLS = ["sample_id", "motif_size", "simple_motif_size", "overlaps_censat", "#chrom"]

res = []
# loop over samples grouped by motif size to start
for (
    sample_id,
    motif_size,
    simple_motif_size,
    in_censat,
    chrom,
), sub_df in tqdm.tqdm(raw_loci.groupby(GROUP_COLS)):
    # denominator is the number of annotated loci of the specified motif size
    # in the specified sample
    denom = sub_df.shape[0]

    res.append(
        {
            "sample_id": sample_id,
            "motif_size": motif_size,
            "simple_motif_size": simple_motif_size,
            "overlaps_censat": in_censat,
            "chrom": chrom,
            "denominator": denom,
        }
    )

res_df = pd.DataFrame(res).merge(total_bp)
res_df.to_csv(snakemake.output.out, sep="\t", index=False)
