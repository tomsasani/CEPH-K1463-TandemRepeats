import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS

import utils

from schema import TRGTDeNovoSchema, TransmissionSchema, PhasedSchema


def main(args):

    DTYPES = {
        "sample_id": str,
        "children_with_denovo_allele": str,
        "children_with_denovo_allele_strict": str,
    }

    # read in the (lightly) prefiltered de novo dataframe.
    filtered_denovos = pd.read_csv(
        args.raw_denovos,
        sep="\t",
        dtype=DTYPES,
    )
    TRGTDeNovoSchema.validate(filtered_denovos)

    # read in the file containing transmission evidence for
    # each TR de novo, if applicable
    if args.transmission_evidence != "UNK":
        transmission_evidence = pd.read_csv(
            args.transmission_evidence,
            sep=";",
            dtype=DTYPES,
        )
        TransmissionSchema.validate(transmission_evidence)

        filtered_denovos = filtered_denovos.merge(transmission_evidence)

    # read in phasing information
    phasing = pd.read_csv(
        args.phasing,
        sep="\t",
        dtype=DTYPES,
    )
    PhasedSchema.validate(phasing)

    mutations = filtered_denovos.merge(phasing, how="left").fillna({"phase_summary": "unknown"})
    mutations["likely_denovo_size_parsimony"] = mutations.apply(
        lambda row: utils.add_likely_de_novo_size(row, use_parsimony=True),
        axis=1,
    )
    mutations["likely_denovo_size"] = mutations.apply(
        lambda row: utils.add_likely_de_novo_size(row, use_parsimony=False),
        axis=1,
    )

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
