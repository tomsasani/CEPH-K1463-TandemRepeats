import pandas as pd
import argparse


def main(args):

    mutations = pd.read_csv(
        args.mutations,
        sep="\t"
    )

    INH_COL = "children_with_denovo_allele"

    if args.transmission != "UNK":

        transmission = pd.read_csv(
            args.transmission,
            sep=";",
            dtype={INH_COL: "str"},
        )

        mutations = mutations.merge(transmission)

    mutations.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--transmission")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
