import pandas as pd
import glob
import matplotlib.pyplot as plt
import utils
from FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS

def assess_manual_validation(row: pd.Series):

    label = []

    if abs(row["likely_denovo_size"]) >= 50:
        if row["trid"] in FAIL_SVS:
            label.append("SV_FAIL")
        else:
            label.append("SV_PASS")
    else:
        label.append("not_SV")

    if row["motif_size"] >= 6:
        if row["trid"] in FAIL_VNTRS:
            label.append("VNTR_FAIL")
        else:
            label.append("VNTR_PASS")
    else:
        label.append("not_VNTR")
    return "|".join(label)

plt.rc("font", size=13)

t2t = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.CHM13v2.element.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    t2t.append(df)
t2t = pd.concat(t2t).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
t2t["assembly"] = "T2T"

hg38 = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.GRCh38.element.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    hg38.append(df)
hg38 = pd.concat(hg38).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})
hg38["assembly"] = "hg38"

mutations = pd.concat([hg38, t2t])

mutations = mutations[mutations["paternal_id"] == 2209]
mutations = mutations[mutations["likely_denovo_size"] != 0]

sample_ids = mutations["sample_id"].unique()

mutations["variant_type"] = "de novo"
mutations["inheritance"] = mutations.apply(lambda row: (
        "paternal"
        if row["phase_summary"].split(":")[0] == "dad"
        else (
            "maternal"
            if row["phase_summary"].split(":")[0] == "mom"
            else "cannot_determine"
        )
    ), axis=1
)

mutations["is_transmitted"] = mutations.apply(
    lambda row: (
        "YES"
        if row["children_with_denovo_allele"] != "unknown"
        else (
            "NO"
            if (
                row["children_with_denovo_allele"] == "unknown"
                and row["sample_id"] in ("2189", "2216")
            )
            else "UNK"
        )
    ), axis=1
)

# element filtering
orthogonal = mutations[
        (mutations["simple_motif_size"] == "STR")
        & (mutations["validation_status"] != "no_data")
    ]
# make sure we're only looking at sites with max AL <= 120
orthogonal["max_al"] = orthogonal["child_AL"].apply(lambda a: max(map(int, a.split(","))))
orthogonal = orthogonal[orthogonal["max_al"] <= 120]

fail_trids = orthogonal[orthogonal["validation_status"] != "pass"]["trid"].unique()

mutations["element_validation_status"] = mutations.apply(
    lambda row: (
        "not_evaluated"
        if row["validation_status"] == "no_data"
        else "fail" if row["trid"] in fail_trids else "pass"
    ),
    axis=1,
)

mutations["manual_validation_status"] = mutations.apply(lambda row: assess_manual_validation(row), axis=1)

mutations["child_to_reference_diff"] = mutations.apply(lambda row: row["denovo_al"] - (row["end"] - row["start"]), axis=1)
mutations["child_to_parent_diff"] = mutations["likely_denovo_size"]

mutations["min_motif_size_in_locus"] = mutations["motif_size"]
mutations["max_motif_size_in_locus"] = mutations["max_motiflen"]

for assembly, assembly_df in mutations.groupby("assembly"):

    true_assembly = "CHM13v2" if assembly == "T2T" else "GRCh38" if assembly == "hg38" else None

    recurrent = pd.read_csv(f"tr_validation/data/recurrent_trids.{true_assembly}.tsv", sep="\t")
    recurrent["contains_G3"] = recurrent["samples_with_denovo"].apply(lambda samples: any([s in sample_ids for s in samples.split(",")]))
    recurrent = recurrent[recurrent["contains_G3"] == True]
    recurrent_trids = recurrent["trid"].unique()
    assembly_df["recurrent_status"] = assembly_df["trid"].apply(lambda trid: "recurrent" if trid in recurrent_trids else "not_recurrent")

    assembly_df["pass_all_filters"] = assembly_df.apply(lambda row: 
        "YES" if (row["element_validation_status"] in ["not_evaluated", "pass"])
        and (not "FAIL" in row["manual_validation_status"])
        and (row["recurrent_status"] != "recurrent") else "NO", axis=1)

    assembly_df[
        [
            "alt_sample_id",
            "trid",
            "#chrom",
            "start",
            "end",
            "motifs",
            "min_motif_size_in_locus",
            "max_motif_size_in_locus",
            "simple_motif_size",
            "likely_denovo_size",
            "variant_type",
            "inheritance",
            "children_with_denovo_allele",
            "is_transmitted",
            "element_validation_status",
            "child_to_reference_diff",
            "child_to_parent_diff",
            "manual_validation_status",
            "recurrent_status",
            "pass_all_filters",
        ]
    ].rename(
        columns={
            "alt_sample_id": "sample_id",
            "#chrom": "chrom",
            "simple_motif_size": "simplified_motif_structure",
            
        }
    ).to_csv(
        f"{assembly}.michelle.tsv", sep="\t", index=False
    )