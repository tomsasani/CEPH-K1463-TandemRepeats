import pandas as pd
import numpy as np
import argparse
from collections import Counter
from typing import List
import matplotlib.pyplot as plt
from schema import HaplotypedSchema
import utils

def measure_consistency(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """given a dataframe that contains information about informative sites
    surrounding a candidate DNM, find the longest stretch of 'consistent'
    informative sites that all support the same haplotype assignment.

    Args:
        df (pd.DataFrame): pandas DataFrame object containing information about informative sites.
        columns (List[str]): list of columns to use in consistency checks

    Returns:
        pd.DataFrame: pandas DataFrame containing a subset of only the informative sites in the
            longest continuous stretch of consistent sites.
    """

    # sort the informative sites by absolute distance to the STR
    df_sorted = df.sort_values(
        "abs_diff_to_str",
        ascending=True,
    ).reset_index(drop=True)#.head(50)

    sorted_values = df_sorted[columns].values

    # NOTE: this is a new approach. if every one of these sites is in the same
    # phase block, why not just take a consensus approach?
    freqs = Counter(sorted_values).most_common()
    most_common_poi = freqs[0][0]
    poi_frac = freqs[0][1] / len(sorted_values)

    # figure out how many sites are consistent closest to the STR. we can do this
    # simply by figuring out the first index where they are *inconsistent.*
    # inconsistent_phases = np.where(sorted_values[1:] != sorted_values[:-1])[0]

    # # if we have more than 1 informative site and none of them are inconsistent,
    # # then all of them are consistent
    # if sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] == 0:
    #     df_sorted_consistent = df_sorted.iloc[:]

    # # if we have more than 1 informative site and some of them are inconsistent,
    # # return the consistent ones up until the first inconsistency.
    # elif sorted_values.shape[0] > 1 and inconsistent_phases.shape[0] > 0:
    #     df_sorted_consistent = df_sorted.iloc[:inconsistent_phases[0]]
    # # if we only have one informative site
    # else:
    #     df_sorted_consistent = df_sorted.iloc[:]

    return most_common_poi, poi_frac





def main(args):

    dtypes = utils.DTYPES.copy()
    dtypes.update({"sample_id": "string"})

    inf_sites = pd.read_csv(args.annotated_dnms, sep="\t", dtype=dtypes)
    HaplotypedSchema.validate(inf_sites)

    print (inf_sites["trid"].nunique())

    COLS = ["trid", "sample_id", "genotype", "index", "suffix"]

    res = []

    # get a dataframe with N rows for each DNM. each of the N rows
    # represents a single informative site
    for (
        trid,
        sample_id,
        genotype,
        index,
        suffix,
    ), df in inf_sites.groupby(COLS):
        
        # subset to only sites with informative info
        
        # if we're on a male sex chromosome, the phase is simple
        is_male_x = suffix.startswith("S") and trid.startswith("chrX")
        is_male_y = suffix.startswith("S") and trid.startswith("chrY")

        poi, poi_support = "unknown", 0
        hoi, hoi_support = "unknown", 0

        if is_male_x:
            poi = "mom"
            poi_support = 1
        elif is_male_y:
            poi = "dad"
            poi_support = 1
        elif df.shape[0] == 0:
            poi = "unknown"
            poi_support = 0

        # otherwise, we have to use informative sites
        else:
            # subset the dataframe to only include the longest stretch of consistent
            # STR parent of origin inferences
            # consistent_poi_sites = measure_consistency(df, "str_parent_of_origin")
            poi, poi_support = measure_consistency(df, "str_parent_of_origin")

            # if consistent_poi_sites.shape[0] == 0:
            #     consensus_phase, consensus_support = "unknown", 0
            # else:
            #     parent_of_origin_support = Counter(consistent_poi_sites["str_parent_of_origin"]).most_common()
            #     assert len(parent_of_origin_support) == 1
            #     consensus_phase, consensus_support = parent_of_origin_support[0]

            # subset the longest stretch of consistent STR parent of origin sites to the ones with
            # the longest contiguous stretch of consistent haplotype inferences
            # hap_sites = consistent_poi_sites[consistent_poi_sites["haplotype_in_parent"] != "UNK"]

            # these sites have to be the same as the ones for which the parent of origin is consistent
            hap_sites = df[(df["haplotype_in_parent"] != "unknown") & (df["str_parent_of_origin"] == poi)]
            if hap_sites.shape[0] == 0:
                hoi = "unknown"
            else:
            # consistent_hap_sites = measure_consistency(
            #     hap_sites,
            #     ["haplotype_in_parent"],
            # )
                hoi, hoi_support = measure_consistency(hap_sites, "haplotype_in_parent")

            # if consistent_hap_sites.shape[0] == 0:
            #     haplotype_in_parent = "UNK"
            # else:
            #     haplotype_of_origin_support = Counter(consistent_hap_sites["haplotype_in_parent"]).most_common()
            #     assert len(haplotype_of_origin_support) == 1
            #     haplotype_in_parent, _ = haplotype_of_origin_support[0]
        df["phase_consensus"] = f"{poi}:{poi_support}"#poi if poi_support >= 0.9 else "unknown"
        df["haplotype_in_parent_consensus"] = f"{hoi}:{hoi_support}"

        res.append(df)

    res_df = (
        pd.concat(res)
        .drop_duplicates(COLS)
        .drop(
            columns=[
                "abs_diff_to_str",
                "inf_chrom",
                "inf_pos",
                "dad_inf_gt",
                "mom_inf_gt",
                "str_parent_of_origin",
                "haplotype_in_parent",
            ]
        )
    )
    print (res_df.shape)

    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
