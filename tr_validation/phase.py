import pandas as pd
import numpy as np
import argparse

def main(args):

    merged_dnms_inf = pd.read_csv(args.annotated_dnms, sep="\t")

    merged_dnms_inf["is_upstream"] = merged_dnms_inf["diff_to_str"] < 0
    merged_dnms_inf["is_downstream"] = merged_dnms_inf["diff_to_str"] > 0

    COLS = ["trid", "sample_id", "genotype", "str_parent_of_origin"]

    merged_dnms_inf_totals = (
        merged_dnms_inf.groupby(COLS)
        .agg({"is_upstream": "sum", "is_downstream": "sum"})
        .reset_index()
        .rename(
            columns={"is_upstream": "n_upstream", "is_downstream": "n_downstream"},
        )
    )

    merged_dnms_inf = merged_dnms_inf.merge(merged_dnms_inf_totals)

    merged_dnms_inf.drop_duplicates(COLS).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)

