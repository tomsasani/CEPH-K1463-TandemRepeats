import pandas as pd
import argparse
import numpy as np

import utils

def main(args):
    raw_denovos = pd.read_csv(
        args.raw_denovos,
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

    multiples = raw_denovos.groupby(["trid", "genotype"]).size().to_dict()
    multiple_trids = [k[0] for k,v in multiples.items() if v > 1]
    for c in raw_denovos.columns:
        raw_vals, trans_vals = raw_denovos[c].values, transmission_evidence[c].values
        if not np.array_equal(raw_vals, trans_vals):
            for i in range(raw_vals.shape[0]):
                if raw_vals[i] != trans_vals[i]:
                    print (raw_vals[i], trans_vals[i])

    raw_plus_transmission = raw_denovos.merge(transmission_evidence)

    # for (trid, genotype), sub_df in raw_plus_transmission.groupby(["trid", "genotype"]):
    #     if sub_df.shape[0] > 1:
    #         for c in sub_df.columns:
    #             if np.unique(sub_df[c].values).shape[0] > 1:
    #                 print (c)

    final = raw_plus_transmission.merge(phasing, how="left").fillna({"phase_summary": "unknown"})

    final["likely_denovo_size"] = final.apply(lambda row: utils.add_likely_de_novo_size(row, use_phase=True), axis=1)

    # ensure we're only looking at sites where the de novo size at least as big as the motif size.
    # if a complex locus contains a homopolymer, allow for 1bp expansions/contractions
    final = final[
        (final["motif_size"].isin([-1, 1])) | (np.abs(final["likely_denovo_size"]) >= 2)
    ]
    final.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--raw_denovos")
    p.add_argument("--transmission_evidence")
    p.add_argument("--phasing")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
