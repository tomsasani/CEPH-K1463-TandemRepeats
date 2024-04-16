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
    raw_loci["motif_size"] = raw_loci["motifs"].apply(lambda m: len(m))
    raw_loci["is_complex"] = raw_loci["n_motifs"].apply(lambda n: n > 1)

    # bin by reference allele length
    # raw_loci["reference_allele_length"] = raw_loci["end"] - raw_loci["start"]
    # raw_loci["binned_reference_allele_length"] = raw_loci

    GROUP_COLS = ["motif_size"]

    res = []
    # loop over samples grouped by motif size to start
    for motif_size, sub_df in tqdm.tqdm(raw_loci.groupby(GROUP_COLS)):
        # denominator is the number of annotated loci of the specified motif size
        # in the specified sample
        denom = sub_df.shape[0]
        
        res.append(
            {
                "motif_size": motif_size[0],
                "denominator": denom,
            }
        )

    res_df = pd.DataFrame(res)
    res_df["sample_id"] = args.sample_id
    print (args.sample_id, res_df["denominator"].sum())

    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--denominator")
    p.add_argument("--annotations")
    p.add_argument("--out")
    p.add_argument("--sample_id")
    p.add_argument("-children_sites", nargs="*", required=False)
    args = p.parse_args()
    main(args)
