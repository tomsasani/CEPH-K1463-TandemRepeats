import tqdm 
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import utils
import numpy as np

plt.rc("font", size=14)

def main(args):

    DTYPES = {
            "trid": str,
            "index": int,
            "child_AL": str,
            "mother_AL": str, 
            "father_AL": str,
            "child_MC": str,
            "per_allele_reads_father": str,
            "per_allele_reads_mother": str,
            "child_coverage": int,
            "denovo_coverage": int,
            "child_ratio": float,
            "father_overlap_coverage": str,
            "mother_overlap_coverage": str,
        }

    # read in raw mutations to define a denominator
    raw_loci = pd.read_csv(
        args.denominator,
        sep="\t",
        usecols=DTYPES.keys(),
        dtype=DTYPES,
    )

    raw_loci = utils.filter_mutation_dataframe(
        raw_loci,
        remove_complex=False,
        remove_duplicates=False,
        remove_gp_ev=False,
        remove_inherited=False,
        parental_overlap_frac_max=1,
        denovo_coverage_min=0,
        depth_min=10,
    )

    if args.trids != "None":
        trids_to_use = []
        with open(args.trids, "r") as infh:
            for l in infh:
                trids_to_use.append(l.strip())
        raw_loci = raw_loci[raw_loci["trid"].isin(trids_to_use)]

    annotations = pd.read_csv(args.annotations, sep="\t")

    raw_loci = raw_loci.merge(annotations, how="left", left_on="trid", right_on="TRid")

    # add some columns we'll use to group by
    raw_loci["motif_size"] = raw_loci.apply(lambda row: utils.determine_motif_size(row), axis=1)
    raw_loci["simple_motif_size"] = raw_loci.apply(lambda row: utils.determine_simplified_motif_size(row), axis=1)

    # add column with binned reference length
    raw_loci["reflen"] = raw_loci["end"] - raw_loci["start"]
    min_reflen, max_reflen = raw_loci["reflen"].min(), raw_loci["reflen"].max()

    reflen_bins = np.arange(min_reflen, max_reflen, 10)
    bins, labels = pd.cut(raw_loci["reflen"].values, bins=reflen_bins, retbins=True)

    raw_loci["reflen_bin"] = bins

    # calculate total number of base pairs
    raw_loci["aff_bp"] = raw_loci["end"] - raw_loci["start"]
    total_bp = raw_loci["aff_bp"].sum()

    # NOTE: if we are using transmission, we need to adjust the denominator to
    # only include sites with >= 10 reads in the G3 individuals' children.
    insufficient_depth_sites = []
    for fh in args.children_sites:
        if fh == "None": continue
        df = pd.read_csv(fh, sep="\t")
        trids = df["trid"].to_list()
        insufficient_depth_sites.extend(trids)
    insufficient_depth_sites = list(set(insufficient_depth_sites))
    raw_loci = raw_loci[~raw_loci["trid"].isin(insufficient_depth_sites)]

    GROUP_COLS = ["motif_size", "simple_motif_size", "reflen_bin", "#chrom"]

    res = []
    # loop over samples grouped by motif size to start
    for (motif_size, simple_motif_size, reflen_bin, chrom), sub_df in tqdm.tqdm(raw_loci.groupby(GROUP_COLS)):
        # denominator is the number of annotated loci of the specified motif size
        # in the specified sample
        denom = sub_df.shape[0]
        
        res.append(
            {
                "motif_size": motif_size,
                "simple_motif_size": simple_motif_size,
                "reflen_bin": reflen_bin,
                "chrom": chrom,
                "denominator": denom,
            }
        )

    res_df = pd.DataFrame(res)
    res_df["sample_id"] = args.sample_id
    res_df["total_bp"] = total_bp
    
    res_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--denominator")
    p.add_argument("--annotations")
    p.add_argument("--out")
    p.add_argument("--sample_id")
    p.add_argument("-children_sites", nargs="*", required=False)
    p.add_argument("-grandparent_sites", nargs="*", required=False)
    p.add_argument("-trids")
    args = p.parse_args()
    main(args)
