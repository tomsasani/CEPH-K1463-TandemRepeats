import pandas as pd
import numpy as np
import argparse
from collections import Counter

from schema import HaplotypedSchema
from annotate_with_informative_sites import measure_consistency

def main(args):

    inf_sites = pd.read_csv(args.annotated_dnms, sep="\t", dtype={"sample_id": str})
    HaplotypedSchema.validate(inf_sites)

    COLS = ["trid", "sample_id", "genotype", "index", "suffix"]

    res = []

    # loop over all sites
    for (trid, sample_id, genotype, index, suffix), df in inf_sites.groupby(COLS):

        # if we're on a male sex chromosome, the phase is simple
        is_male_x = suffix.startswith("S") and trid.startswith("chrX")
        is_male_y = suffix.startswith("S") and trid.startswith("chrY")

        consensus_phase, consensus_support = None, 0
        al_in_parent = -1

        if is_male_x:
            consensus_phase = "mom"
            consensus_support = 99
        elif is_male_y:
            consensus_phase = "dad"
            consensus_support = 99
        
        # otherwise, we have to use informative sites
        else:
            # measure the consistency of str parent of origin inference
            consistent_poi_sites = measure_consistency(df, "str_parent_of_origin")
            # measure consistency of the haplotype inferences, ignoring informative sites at which it's unknown
            hap_df = df[df["haplotype_in_parent"] != "UNK"]
            consistent_hap_sites = measure_consistency(hap_df, ["haplotype_in_parent", "allele_length_in_parent"])
            
            if consistent_poi_sites.shape[0] <= 1:
                consensus_phase, consensus_support = "unknown", 0
            else:
                parent_of_origin_support = Counter(consistent_poi_sites["str_parent_of_origin"]).most_common()
                assert len(parent_of_origin_support) == 1
                consensus_phase, consensus_support = parent_of_origin_support[0]

            if consistent_hap_sites.shape[0] <= 1:
                al_in_parent = -1
            else:
                haplotype_of_origin_support = Counter(consistent_hap_sites["allele_length_in_parent"]).most_common()
                assert len(haplotype_of_origin_support) == 1
                al_in_parent, _ = haplotype_of_origin_support[0]
        
        df["phase_consensus"] = f"{consensus_phase}:{consensus_support}"
        df["allele_length_consensus"] = al_in_parent

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
                "allele_length_in_parent",
            ]
        )
    )

    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
