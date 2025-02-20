import pandas as pd
import argparse
import numpy as np
import glob

from FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS

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
    
    raw_plus_transmission = raw_denovos.merge(transmission_evidence)

    mutations = raw_plus_transmission.merge(phasing, how="left").fillna({"phase_summary": "unknown"})
    mutations["likely_denovo_size"] = mutations.apply(lambda row: utils.add_likely_de_novo_size(row, use_phase=True), axis=1)


    # if we 're not at a homopolymer, ensure that the likely de novo size is at least 2 bp
    mutations = mutations[
        (mutations["motif_size"].isin([-1, 1])) | (np.abs(mutations["likely_denovo_size"]) >= 2)
    ]
    mutations = mutations[mutations["likely_denovo_size"] != 0]

    # if desired, filter SVs that failed manual inspection
    mutations = mutations[
        (
            (np.abs(mutations["likely_denovo_size"]) >= 50)
            & (~mutations["trid"].isin(FAIL_SVS))
        )
        | (np.abs(mutations["likely_denovo_size"]) < 50)
    ]

    # if desired, filter VNTRs that failed manual inspection
    mutations = mutations[
        (
            (mutations["simple_motif_size"] == "VNTR")
            & (~mutations["trid"].isin(FAIL_VNTRS))
        )
        | (mutations["simple_motif_size"] != "VNTR")
    ]

    mutations.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--raw_denovos")
    p.add_argument("--transmission_evidence")
    p.add_argument("--phasing")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
