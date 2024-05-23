import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import utils
from sklearn.mixture import GaussianMixture
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance


ASSEMBLY = "GRCh38"
# ASSEMBLY = "CHM13v2"

# get mutations
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.element.tsv"):
    df = pd.read_csv(
        fh,
        sep="\t",
        dtype={"sample_id": str},
    )
    df = df[df["motif"].str.len() <= 6]
    mutations.append(df)

mutations = pd.concat(mutations)
mutations = mutations[mutations["likely_denovo_size"] != 0]
mutations = mutations[mutations["validation_status"] != "no_data"]


validation_outfh = open("tr_validation/validation.bed", "w")


for (trid, genotype), trid_df in mutations.groupby(["trid", "genotype"]):

    validation = trid_df["validation_status"].unique()[0]
    motif_size = trid_df["motif"].str.len().values[0]
    sample = trid_df["sample_id"].unique()[0]

    chrom, start, end, _ = trid.split("_")
    print ("\t".join([chrom, start, end, validation]), file=validation_outfh)

    if validation in ("fail", "not_ibs"):
        print (sample, trid, motif_size)

    trid_df["denovo_al"] = trid_df.apply(lambda row: row["child_AL"].split(",")[row["index"]], axis=1)
    trid_df["non_denovo_al"] = trid_df.apply(lambda row: row["child_AL"].split(",")[1 - row["index"]], axis=1)

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
    exp_non_denovo = trid_df["non_denovo_al"].unique()[0]

    f, ax = plt.subplots()

    sns.stripplot([kid_diffs, mom_diffs, dad_diffs], palette="colorblind", ax=ax, alpha=0.75)
    ax.axhline(exp_denovo_diff, ls=":", c="dodgerblue", label="Expected allele length (de novo)")
    ax.axhline(exp_non_denovo_diff, ls=":", c="firebrick", label="Expected allele length (non de novo)")
    ax.set_ylabel("Allele length inferred using Element reads")
    ax.set_xlabel("Sample")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Kid", "Mom", "Dad"])
    ax.set_title(f"{trid}\n{validation} - {motif_size}")
    ax.legend()
    sns.despine(ax=ax)
    f.savefig(f'tr_validation/fig/denovo/{sample}.{trid}.{validation}.png', dpi=200)
    