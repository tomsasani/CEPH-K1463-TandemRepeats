import pandas as pd
from typing import List
from collections import defaultdict

from schema import OrthogonalSchema

def annotate_with_concordance(row: pd.Series) -> str:
    """given orthogonal evidence at a locus, figure out whether
    the orthogonal evidence is concordant with the expected genotype
    in the child.

    Args:
        row (pd.Series): pd.Series object

    Returns:
        str: annotation label
    """

    # gather diffs in all members of the trio
    mom_diffs, dad_diffs, kid_diffs = [], [], []
    for col in ("mom_evidence", "dad_evidence", "kid_evidence"):
        for diff_count in row[col].split("|"):
            diff, count = list(map(int, diff_count.split(":")))
            if col == "mom_evidence":
                mom_diffs.extend([diff] * count)
            elif col == "dad_evidence":
                dad_diffs.extend([diff] * count)
            else:
                kid_diffs.extend([diff] * count)

    # figure out the difference between the expected
    # sizes of the de novo and non-de novo alleles
    denovo_al = int(row["exp_allele_diff_denovo"])
    non_denovo_al = int(row["exp_allele_diff_non_denovo"])
    exp_diff = denovo_al - non_denovo_al
    # how many alleles do we expect to see in the kid? i.e., is the kid HOM or HET?
    n_alleles_exp = 1 if exp_diff == 0 else 2

    # figure out the observed number of distinct alleles in the kid
    n_alleles = len(set(kid_diffs))
    if n_alleles_exp > 1 and n_alleles == 1:
        return "not_ibs"

    # validate by asking if the de novo allele is observed in the
    # child and absent from both parents
    out_dict = defaultdict(int)
    out_dict_off_by_one = defaultdict(int)

    for diffs, label in zip(
        (kid_diffs, mom_diffs, dad_diffs),
        ["kid", "mom", "dad"],
    ):
        dn_overlap, dn_overlap_off = 0, 0
        # allow for off-by-one errors when asking if kid has de novo.
        # don't allow for off-by-one errors when asking if parents have de novo.
        for d in diffs:
            if denovo_al == d: dn_overlap += 1
            if d - 1 <= denovo_al <= d + 1: dn_overlap_off += 1
        out_dict[label] = dn_overlap
        out_dict_off_by_one[label] = dn_overlap_off

    return (
        "pass"
        if (
            out_dict_off_by_one["kid"] >= 1
            and out_dict_off_by_one["mom"] == 0
            and out_dict_off_by_one["dad"] == 0
        )
        else "fail"
    )


mutations = pd.read_csv(
    snakemake.input.mutations,
    sep="\t",
    dtype={"sample_id": str},
)
ortho_evidence = pd.read_csv(
    snakemake.input.orthogonal_evidence,
    dtype={"sample_id": str},
    sep="\t",
)
OrthogonalSchema.validate(ortho_evidence)

res: List[pd.DataFrame] = []

for i, row in ortho_evidence.iterrows():
    row_dict = row.to_dict()
    parental_overlap = annotate_with_concordance(row)
    row_dict.update({"validation_status": parental_overlap})
    res.append(row_dict)

ortho_validation = pd.DataFrame(res)

# merge mutations with orthogonal validation
mutations = mutations.merge(ortho_validation, how="left").fillna(
    {
        "validation_status": "no_data",
    }
)

mutations.to_csv(snakemake.output.out, sep="\t", index=False)



