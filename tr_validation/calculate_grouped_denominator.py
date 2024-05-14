import seaborn as sns
import numpy as np
import seaborn.objects as so
import tqdm 
import glob 
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import utils

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
        remove_duplicates=True,
        remove_gp_ev=False,
        remove_inherited=False,
        parental_overlap_frac_max=1,
        denovo_coverage_min=0,
        depth_min=10,
    )

    annotations = pd.read_csv(args.annotations, sep="\t")

    raw_loci = raw_loci.merge(annotations, left_on="trid", right_on="TRid")

    # add some columns we'll use to group by
    raw_loci["motif_size"] = raw_loci.apply(lambda row: len(row["motifs"]) if row["n_motifs"] == 1 else -1, axis=1)

    print (raw_loci.query("motif_size == 1").groupby("motifs").size())

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

    # NOTE: if we are filtering on grandparental evidence, we need to adjust the denominator to
    # only include sites with >= 10 reads in the G3 individuals' grandparents.
    insufficient_depth_sites = []
    for fh in args.grandparent_sites:
        df = pd.read_csv(fh, sep="\t")
        trids = df["trid"].to_list()
        insufficient_depth_sites.extend(trids)
    insufficient_depth_sites = list(set(insufficient_depth_sites))
    raw_loci = raw_loci[~raw_loci["trid"].isin(insufficient_depth_sites)]

    raw_loci["chrom"] = raw_loci["trid"].apply(lambda t: t.split("_")[0])

    GROUP_COLS = ["motif_size", "chrom"]#, "insufficient_depth_in_children", "insufficient_depth_in_grandparents"]

    res = []
    # loop over samples grouped by motif size to start
    for (motif_size, chrom), sub_df in tqdm.tqdm(raw_loci.groupby(GROUP_COLS)):
        # denominator is the number of annotated loci of the specified motif size
        # in the specified sample
        denom = sub_df.shape[0]
        
        res.append(
            {
                "motif_size": motif_size,
                "chrom": chrom,
                "denominator": denom,
            }
        )

    res_df = pd.DataFrame(res)
    res_df["sample_id"] = args.sample_id

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
    args = p.parse_args()
    main(args)
