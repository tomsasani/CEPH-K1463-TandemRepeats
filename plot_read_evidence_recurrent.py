
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv

def assign_pedigree_id(row: pd.Series):
    comb_id = []
    if row["sample_id_with_evidence"] in ("2280", "2281") or row["paternal_id_y"] == "2281":
        comb_id.append("G2A")
    if row["sample_id_with_evidence"] in ("2214", "2213") or row["paternal_id_y"] == "2214":
        comb_id.append("G2B")
    if row["sample_id_with_evidence"] in ("2209", "2188") or row["paternal_id_y"] == "2209":
        comb_id.append("G3")
    if row["sample_id_with_evidence"] in ("2216", "200080") or row["paternal_id_y"] == "200080":
        comb_id.append("G4A")
    if row["sample_id_with_evidence"] in ("2189", "200100") or row["paternal_id_y"] == "2189":
        comb_id.append("G4B")

    return ",".join(comb_id)


ASSEMBLY = "CHM13v2"
TECH = "hifi"

# get mutations
mutations = []
for fh in glob.glob(f"csv/orthogonal_support/all/*.{ASSEMBLY}.{TECH}.v4.0.tsv"):
    sample_id = fh.split("/")[-1].split(".")[0]
    df = pd.read_csv(fh, sep="\t")#, dtype={"sample_id": str, "paternal_id": str, "maternal_id": str})
    df["sample_id_with_evidence"] = sample_id
    mutations.append(df)
mutations = pd.concat(mutations).dropna(subset=["kid_evidence"])

recurrents = pd.read_csv("csv/recurrent/CHM13v2.v4.0.recurrent.tsv", sep="\t")["trid"].to_list()
mutations = mutations[mutations["trid"].isin(recurrents)]

mutations["reference_al"] = mutations["end"] - mutations["start"]

mutations["denovo_al"] = mutations.apply(lambda row: list(map(int, row["child_AL"].split(",")))[row["index"]], axis=1)


ped = pd.read_csv("tr_validation/data/file_mapping.csv", dtype={"sample_id": str, "paternal_id": str, "maternal_id": str})
mutations = mutations.merge(ped, left_on="sample_id_with_evidence", right_on="sample_id", how="left")

# need evidence in all members of pedigree
count_per_trid = mutations.drop_duplicates(["sample_id_with_evidence", "trid"]).groupby("trid").size().to_dict()

good_trids = {k for k,v in count_per_trid.items() if v == mutations["sample_id_with_evidence"].nunique()}

mutations["pedigree_id"] = mutations.apply(lambda row: assign_pedigree_id(row), axis=1)

mutations["is_complex"] = mutations["motifs"].apply(lambda m: "," in m)

mutations["generation"] = mutations["sample_id_with_evidence"].apply(
    lambda s: (
        "2"
        if s in ("2209", "2188")
        else (
            "4"
            if (s.startswith("200") and s not in ("200080", "200100"))
            else "1" if s in ("2281", "2280", "2213", "2214") else "3"
        )
    )
)

mutations["parent_status"] = mutations["sample_id_with_evidence"].apply(lambda s: (
            "G4A"
            if s in ("2216", "200080") else "G4B" if s in ("2189", "200100")
            else (
                "G3"
                if s in ("2209", "2188")
                else "G2A" if s in ("2280", "2281") else "G2B" if s in ("2214", "2213") else "UNK"
            )
        ))


for i, (trid, trid_df) in enumerate(mutations.groupby("trid")):
    n_samples = trid_df["sample_id_with_evidence"].nunique()
    print (n_samples)
    if n_samples != 28: continue

    samples_with_denovo = list(map(str, trid_df["sample_id_x"].unique()))
    print (samples_with_denovo)
    if not all([s.startswith("200") for s in samples_with_denovo]): continue
    res = []

    # gather diffs in all members of the trio
    for sample_id, sample_df in trid_df.groupby("sample_id_with_evidence"):
        sample_df = sample_df.drop_duplicates("trid")
        has_denovo = str(sample_id) in samples_with_denovo

        assert sample_df.shape[0] == 1

        # get pedigree IDs of this individual (parents have multiple)
        pedigree_ids = sample_df["pedigree_id"].unique()[0]
        generation = sample_df["generation"].unique()[0]
        parent_status = sample_df["parent_status"].unique()[0]
        for pedigree_id in pedigree_ids.split(","):

            # is this individual a parent within the pedigree ID
            is_parent = parent_status == pedigree_id
            is_child = generation == pedigree_id[1]

            for diff_count in sample_df["kid_evidence"].values[0].split("|"):
                diff, count = list(map(int, diff_count.split(":")))
                for _ in range(count):
                    if is_parent:
                        res.append(
                            {
                                "sample_id": sample_id,
                                "diff": diff,
                                "status": "parent",
                                "generation": pedigree_id,
                                "has_denovo": has_denovo,
                            }
                        )
                    if is_child:
                        res.append(
                            {
                                "sample_id": sample_id,
                                "diff": diff,
                                "status": "child",
                                "generation": pedigree_id,
                                "has_denovo": has_denovo,
                            }
                        )

    res_df = pd.DataFrame(res)

    #res_df = res_df[res_df["generation"].isin(["G3"])]

    f, axarr = plt.subplots(res_df["generation"].nunique(), figsize=(8, 18), sharex=True)
    for gen_i, (gen, gen_df) in enumerate(res_df.groupby("generation")):
        sns.stripplot(data=gen_df.sort_values("status", ascending=True), y="sample_id", x="diff", hue="has_denovo", alpha=0.5, ax=axarr[gen_i])
        axarr[gen_i].set_ylabel("Sample ID")
        axarr[gen_i].set_xlabel(f"Allele length (w/r/t {ASSEMBLY} genome)")
        sns.despine(ax=axarr[gen_i])
        axarr[gen_i].set_title(gen)
    f.tight_layout()
    f.savefig(f'fig/recurrent/{ASSEMBLY}.{trid}.png', dpi=200)
    plt.close()
