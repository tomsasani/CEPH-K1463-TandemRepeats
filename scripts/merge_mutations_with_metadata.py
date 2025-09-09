import pandas as pd
import numpy as np
from FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS

import utils

from schema import TRGTDeNovoSchema, TransmissionSchema, PhasedSchema


DTYPES = {
    "sample_id": "string",
    "trid": "string",
    "index": "int",
    "children_with_denovo_allele": "string",
    "grandparents_with_denovo_allele": "string",
}

# read in the (lightly) prefiltered de novo dataframe.
filtered_denovos = pd.read_csv(
    snakemake.input.raw_denovos,
    sep="\t",
    dtype=DTYPES,
)
TRGTDeNovoSchema.validate(filtered_denovos)

# read in the file containing transmission evidence for
# each TR de novo, if applicable
transmission_evidence = pd.read_csv(
    snakemake.input.transmission,
    sep="\t",
    dtype=DTYPES,
    usecols=["sample_id", "trid", "index", "children_with_denovo_allele"],
)
TransmissionSchema.validate(transmission_evidence)

filtered_denovos = filtered_denovos.merge(transmission_evidence)
assert transmission_evidence.shape[0] == filtered_denovos.shape[0]

# read in the file containing transmission evidence for
# each TR de novo, if applicable
grandparental_evidence = pd.read_csv(
    snakemake.input.grandparental,
    sep="\t",
    dtype=DTYPES,
    usecols=["sample_id", "trid", "index", "grandparents_with_denovo_allele"],
)
filtered_denovos = filtered_denovos.merge(grandparental_evidence)

# read in phasing information
phasing = pd.read_csv(
    snakemake.input.phasing,
    sep="\t",
    dtype={"sample_id": "string"},
)
PhasedSchema.validate(phasing)


# NOTE: this only works if the transmission data matches up with new DNM data
assert phasing.shape[0] == filtered_denovos.shape[0]

mutations = filtered_denovos.merge(phasing)

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
    (mutations["motif_size"].isin([-1, 1])) | (mutations["likely_denovo_size"].abs() >= 2)
]
assert mutations.query("motif_size > 1 and likely_denovo_size == 1").shape[0] == 0

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

mutations.to_csv(snakemake.output.out, index=False, sep="\t")
