import pandas as pd
import numpy as np
import argparse
from collections import Counter

def main(args):

    merged_dnms_inf = pd.read_csv(args.annotated_dnms, sep="\t")

    merged_dnms_inf["is_upstream"] = merged_dnms_inf["diff_to_str"] < 0
    merged_dnms_inf["is_downstream"] = merged_dnms_inf["diff_to_str"] > 0

    COLS = ["trid", "sample_id", "genotype"]

    merged_dnms_inf_totals = (
        merged_dnms_inf.groupby(COLS)
        .agg({"is_upstream": "sum", "is_downstream": "sum", "str_parent_of_origin": lambda p: ":".join(list(map(str, Counter(p).most_common()[0])))})
        .reset_index()
        .rename(
            columns={"is_upstream": "n_upstream", "is_downstream": "n_downstream", "str_parent_of_origin": "phase_summary"},
        )
    )

    print (merged_dnms_inf_totals.groupby(COLS).size().sort_values())

    # print (merged_dnms_inf_totals[merged_dnms_inf_totals["trid"] == "chr7_150033538_150042643_trsolve"])

    merged_dnms_inf = merged_dnms_inf.merge(merged_dnms_inf_totals).drop(
        columns=[
            "is_upstream",
            "is_downstream",
            "diff_to_str",
            "abs_diff_to_str",
            "inf_chrom",
            "inf_pos",
            "dad_inf_gt",
            "mom_inf_gt",
            "kid_inf_gt",
            "kid_inf_ps",
            "haplotype_A_origin",
            "haplotype_B_origin",
            "dad_haplotype_origin",
            "mom_haplotype_origin",
            "str_midpoint",
        ]
    )

    merged_dnms_inf.drop_duplicates(COLS).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
