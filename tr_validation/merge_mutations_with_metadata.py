import pandas as pd
import argparse

import utils

def main(args):
    raw_denovos = pd.read_csv(
        args.raw_denovos,
        sep="\t",
        dtype={"sample_id": str},
    )

    grandparental_evidence = pd.read_csv(
        args.grandparental_evidence,
        sep="\t",
        dtype={"sample_id": str},
    )

    transmission_evidence = pd.read_csv(
        args.transmission_evidence,
        sep="\t",
        dtype={"sample_id": str},
    )

    # read in phasing information
    phasing = pd.read_csv(
        args.phasing,
        sep="\t",
        dtype={"sample_id": str},
    )


    raw_plus_transmission = raw_denovos.merge(transmission_evidence)
    raw_plus_transmission_plus_grandparental = raw_plus_transmission.merge(grandparental_evidence)


    final = raw_plus_transmission_plus_grandparental.merge(phasing, how="left").fillna({"phase_summary": "unknown"})

    final["likely_denovo_size"] = final.apply(lambda row: utils.add_likely_de_novo_size(row), axis=1)

    final.to_csv(args.out, index=False, sep="\t")
    print (raw_denovos.shape, transmission_evidence.shape, final.shape)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--raw_denovos")
    p.add_argument("--grandparental_evidence")
    p.add_argument("--transmission_evidence")

    p.add_argument("--phasing")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
