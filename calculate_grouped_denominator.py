import tqdm 
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import utils
import numpy as np

plt.rc("font", size=14)

def main(args):

    # read in raw mutations to define a denominator
    raw_loci = pd.read_csv(args.loci, sep="\t")

    # add column with binned reference length
    raw_loci["reflen"] = raw_loci["end"] - raw_loci["start"]
    min_reflen, max_reflen = raw_loci["reflen"].min(), raw_loci["reflen"].max()

    reflen_bins = np.arange(min_reflen, max_reflen, 10)
    bins, labels = pd.cut(raw_loci["reflen"].values, bins=reflen_bins, retbins=True)

    raw_loci["reflen_bin"] = bins

    # calculate total number of base pairs
    raw_loci["aff_bp"] = raw_loci["end"] - raw_loci["start"]
    total_bp = raw_loci["aff_bp"].sum()

    GROUP_COLS = ["sample_id", "motif_size", "simple_motif_size", "reflen_bin", "#chrom"]

    res = []
    # loop over samples grouped by motif size to start
    for (sample_id, motif_size, simple_motif_size, reflen_bin, chrom), sub_df in tqdm.tqdm(raw_loci.groupby(GROUP_COLS)):
        # denominator is the number of annotated loci of the specified motif size
        # in the specified sample
        denom = sub_df.shape[0]

        res.append(
            {
                "sample_id": sample_id,
                "motif_size": motif_size,
                "simple_motif_size": simple_motif_size,
                "reflen_bin": reflen_bin,
                "chrom": chrom,
                "denominator": denom,
            }
        )

    res_df = pd.DataFrame(res)
    res_df["total_bp"] = total_bp

    res_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--loci")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
