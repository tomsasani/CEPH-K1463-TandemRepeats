import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob


ASSEMBLY = "CHM13v2"
ORTHOGONAL_TECH = "hifi"
MIN_SIZE = 1

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna({"children_with_denovo_allele": "unknown", "phase_summary": "unknown"})

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == 2209]

# FILTERING
mutations = mutations[np.abs(mutations["likely_denovo_size"]) >= MIN_SIZE]
mutations = mutations[mutations["simple_motif_size"] == "VNTR"]

mutations = mutations[mutations["validation_status"] != "no_data"]

mutations["reference_al"] = mutations["end"] - mutations["start"]

# calculate child-to-reference diff
mutations["child_to_reference_diff"] = (
    mutations["denovo_al"] - mutations["reference_al"]
)
mutations["child_to_parent_diff"] = mutations["likely_denovo_size"]


mutations = mutations.sort_values(["alt_sample_id", "trid"])

for (sample, trid, genotype), trid_df in mutations.groupby(
    ["alt_sample_id", "trid", "genotype"]
):
    validation = trid_df["validation_status"].unique()[0]
    motif_size = trid_df["motif_size"].values[0]
    poi = trid_df["phase_summary"].unique()[0]
    size = trid_df["likely_denovo_size"].unique()[0]


    if poi != "unknown":
        poi = poi.split(":")[0] if int(poi.split(":")[1]) > 5 else "unknown"

    chrom, start, end, _ = trid.split("_")

    is_male_sex_chrom = "," not in trid_df["child_AL"].unique()[0]


    trid_df["denovo_al"] = trid_df.apply(
        lambda row: row["child_AL"].split(",")[row["index"]], axis=1
    )
    if not is_male_sex_chrom:
        trid_df["non_denovo_al"] = trid_df.apply(
            lambda row: row["child_AL"].split(",")[1 - row["index"]], axis=1
        )

    assert trid_df.shape[0] == 1

    mom_diffs, dad_diffs, kid_diffs = [], [], []
    for col in ("mom_evidence", "dad_evidence", "kid_evidence"):
        for diff_count in trid_df[col].values[0].split("|"):
            diff, count = list(map(int, diff_count.split(":")))
            if col == "mom_evidence":
                mom_diffs.extend([diff] * count)
            elif col == "dad_evidence":
                dad_diffs.extend([diff] * count)
            else:
                kid_diffs.extend([diff] * count)

    exp_denovo_diff = trid_df["exp_allele_diff_denovo"].unique()[0]
    exp_non_denovo_diff = trid_df["exp_allele_diff_non_denovo"].unique()[0]

    exp_denovo = trid_df["denovo_al"].unique()[0]
    if not is_male_sex_chrom:
        exp_non_denovo = trid_df["non_denovo_al"].unique()[0]

    f, ax = plt.subplots()

    plot_df = []
    for diff_list, label in zip(
        (kid_diffs, mom_diffs, dad_diffs), ("kid", "mom", "dad")
    ):
        is_poi = None
        if label == "kid":
            is_poi = "kid"
        else:
            if poi == "unknown":
                is_poi = False
            else:
                is_poi = poi == label
        for diff in diff_list:
            plot_df.append({"sample": label, "is_poi": is_poi, "diff": diff})
    plot_df = pd.DataFrame(plot_df)
    sns.stripplot(
        data=plot_df,
        x="sample",
        y="diff",
        ax=ax,
        hue="is_poi",
        palette={True: "green", False: "red", "kid": "blue"},
    ) 
    ax.axhline(
        exp_denovo_diff,
        ls=":",
        c="dodgerblue",
        label="Expected allele length (de novo)",
    )
    if not is_male_sex_chrom:
        ax.axhline(
            exp_non_denovo_diff,
            ls=":",
            c="firebrick",
            label="Expected allele length (non de novo)",
        )
    ax.set_ylabel(f"Allele length inferred using {ORTHOGONAL_TECH} reads")
    ax.set_xlabel("Sample")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Kid", "Mom", "Dad"])
    ax.set_title(f"{trid}\nMotif size = {motif_size}\nDNM size = {size}")
    # ax.legend()
    sns.despine(ax=ax)
    f.tight_layout()
    f.savefig(
        f"tr_validation/fig/denovo/{ORTHOGONAL_TECH}/{ASSEMBLY}/{sample}.{trid}.png",
        dpi=200,
    )
    plt.close()
