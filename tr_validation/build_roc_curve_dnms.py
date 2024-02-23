import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable
import numba
from schema import DeNovoSchema
import cyvcf2
from collections import Counter


from plot_read_support import classify_motif, calculate_min_dist
from utils import get_read_diff
from assess_concordance import als_in_stdev

from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

plt.rc("font", size=10)


def calculate_expected_diffs(
    vcf: cyvcf2.VCF,
    row: pd.Series,
):
    exp_diff_ref, exp_diff_alt = 0, 0

    try:
        ref_al, alt_al = list(map(int, row["child_motif_counts"].split(",")))
    except ValueError:
        return ",".join(list(map(str, [exp_diff_ref, exp_diff_alt])))

    for var in vcf(row["region"]):
        var_trid = var.INFO.get("TRID")
        if row["trid"] != var_trid:
            continue
        # get information about a variant from a VCF
        exp_diff_alt = alt_al - len(var.REF)
        exp_diff_ref = ref_al - len(var.REF)

    return ",".join(list(map(str, [exp_diff_ref, exp_diff_alt])))


def trid_to_motif(
    vcf: cyvcf2.VCF,
    row: pd.Series,
):
    motif = "N"
    for var in vcf(row["region"]):
        var_trid = var.INFO.get("TRID")
        if row["trid"] != var_trid:
            continue
        # get information about a variant from a VCF
        motif = var.INFO.get("MOTIFS")

    return motif


def trid_to_region(t):
    try:
        chrom, start, end, method = t.split("_")
        start, end = int(start) - 1, int(end)
        return f"{chrom}:{start}-{end}"
    except ValueError:
        return "UNK"


def format_validation_df(df, tech: str = "element"):
    df = df[df["sample"] == "kid"]
    df_totals = (
        df.groupby("region")
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": f"total_{tech}_reads"})
    )
    df = df.merge(df_totals, on="region")

    df["Motif type"] = df["motif"].apply(lambda m: classify_motif(m))

    # for each read, calculate the minimum distance to either the REF or ALT allele
    # (i.e., the "closest" allele)
    df["is_ref"] = df.apply(lambda row: calculate_min_dist(row), axis=1)
    df["Expected allele length"] = df.apply(
        lambda row: row["exp_allele_diff_ref"]
        if row["is_ref"] == "REF"
        else row["exp_allele_diff_alt"],
        axis=1,
    )
    df[f"{tech}_error"] = df["diff"] - df["Expected allele length"]
    df[f"frac_{tech}_reads"] = df["Read support"] / df[f"total_{tech}_reads"]
    return df


def get_ab(par: str):
    per_allele_reads = list(map(int, par.split(",")))
    if len(per_allele_reads) == 1:
        return 1
    else:
        return per_allele_reads[-1] / sum(per_allele_reads)


def get_dp(par: str):
    per_allele_reads = list(map(int, par.split(",")))
    return sum(per_allele_reads)


def main():
    # read in raw DNM file
    mutations = pd.read_csv(
        "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/2189_SF_GRCh38_50bp_merge_trgtdn.T95.csv.gz",
        sep="\t",
    )

    KID_STR_VCF = cyvcf2.VCF(
        "/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/data/hiphase/2189_SF_GRCh38_50bp_merge.sorted.phased.vcf.gz",
        gts012=True,
    )

    # sample_ids = ["2189"]
    # # filter to the samples of interest
    # mutations = mutations[mutations["sample_id"].isin(sample_ids)]

    # if we're looking at de novos, ensure that we filter
    # the DataFrame to exclude the non-denovo candidates
    # remove denovos where allele size is unchanged
    mutations = mutations[mutations["denovo_coverage"] >= 1]
    mutations = mutations[~mutations["trid"].str.contains("X|Y|Un|random", regex=True)]

    print (mutations.shape)

    # remove pathogenics
    mutations["region"] = mutations["trid"].apply(lambda t: trid_to_region(t))
    mutations = mutations[mutations["region"] != "UNK"]

    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )

    # merge with read support evidence
    support = pd.read_csv("tr_validation/csv/2189.element.dnm.read_support.csv")[
        [
            "trid",
            "region",
            "exp_allele_diff_ref",
            "exp_allele_diff_alt",
            "sample",
            "diff",
            "Read support",
        ]
    ]

    support_totals = (
        support.groupby("region")
        .agg({"Read support": "sum"})
        .reset_index()
        .rename(columns={"Read support": "Total read support"})
    )
    support = support.merge(support_totals, on="region", how="left")
    # require at least 10 total spanning reads
    support = support[(support["Total read support"] >= 10) & (support["Total read support"] <= 100)]

    mutations = mutations.merge(support, on=["trid", "region"])

    mutations = mutations[mutations["sample"] == "kid"]

    mutations["exp_diff"] = mutations["exp_allele_diff_ref"] - mutations["exp_allele_diff_alt"]
    mutations = mutations.query("exp_diff != 0")

    res = []

    for trid, kid_df_trid in mutations.groupby("trid"):
        is_concordant = als_in_stdev(kid_df_trid)
        kid_df_trid["is_concordant"] = is_concordant 
        res.append(kid_df_trid)

    mutations = pd.concat(res)

    print (mutations.columns)

    # remove complex STRs
    mutations = mutations[~mutations["child_motif_counts"].str.contains("_")]
    mutations["motif"] = mutations.apply(lambda row: trid_to_motif(KID_STR_VCF, row), axis=1)
    mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

    mutations["expansion_size"] = mutations["child_motif_counts"].apply(lambda mc: abs(int(mc.split(",")[1]) - int(mc.split(",")[0])) if len(mc.split(",")) > 1 else np.nan)
    mutations.dropna(subset=["expansion_size"], inplace=True)

    candidates = mutations["trid"].nunique()
    print(f"{candidates} loci that are candidate DNMs")

    # add a transmission column if it's not in there
    # require that the samples have children
    # mutations["expected_children"] = mutations["sample_id"].apply(lambda s: 6 if s in ("2189", "2216") else 8 if s in ("2209", "2188") else np.nan)
    # mutations = mutations[mutations["n_children"] == mutations["expected_children"]]

    # MIN_TRANSMISSIONS = 1
    # if "n_with_denovo_allele_strict" in mutations.columns:
    #     mutations["Is transmitted?"] = mutations["n_with_denovo_allele"].apply(
    #         lambda n: 1 if n >= MIN_TRANSMISSIONS else 0
    #     )
    # else:
    #     print("No transmission column!")
    #     mutations["Is transmitted?"] = "unknown"

    non_inf = np.where(~np.isinf(mutations["allele_ratio"].values))[0]
    mutations = mutations.iloc[non_inf]

    FEATURES = [
        "child_ratio",
        "allele_ratio",
        "denovo_coverage",
        "allele_coverage",
        # "child_coverage",
        # "mom_dp",
        # "dad_dp",
        "mean_diff_father",
        "mean_diff_mother",
        "father_dropout_prob",
        "mother_dropout_prob",
        # "motif_size",
        "expansion_size",
    ]

    X = mutations[FEATURES]
    y = mutations["is_concordant"].values

    f, ax = plt.subplots()
    sns.heatmap(X.corr(method="spearman"), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_title("Correlation between features")
    ax.set_yticks(np.arange(len(FEATURES)) + 0.5)
    ax.set_yticklabels(FEATURES, rotation=0)
    f.tight_layout()
    f.savefig("corr.png", dpi=200)

    X_train, X_test, y_train, y_test = train_test_split(
        X.values, y, shuffle=True, test_size=0.2
    )

    f, axarr = plt.subplots(
        len(FEATURES),
        2,
        figsize=(12, (4 * len(FEATURES))),
        sharex=False,
        sharey=False,
    )

    for fi, feat in enumerate(FEATURES):
        # get feature values for this feature
        X_feat = X.values[:, fi]

        # infs
        # inf_idxs = np.where(np.isinf(X_feat))
        # non_inf_idxs = np.where(~np.isinf(X_feat))
        # X_feat = X_feat[non_inf_idxs]
        # y_prime = y[non_inf_idxs]

        bins = np.linspace(np.min(X_feat), np.max(X_feat), num=50)
        bins = np.arange(np.min(X_feat), np.max(X_feat), step=1 if np.max(X_feat) > 1 else 0.05)

        for y_ in np.unique(y):
            y_i = np.where(y == y_)[0]
            X_feat_sub = X_feat[y_i]

            histvals, edges = np.histogram(X_feat_sub, bins=bins)
            histfracs = histvals / np.sum(histvals)
            cumsum = np.cumsum(histfracs)

            label = "TP" if y_ == 1 else "TN"

            # axarr[fi, 0].plot(
            #     bins[1:],
            #     cumsum,
            #     label=label,
            #     lw=2,
            # )
            axarr[fi, 0].hist(
                X_feat_sub,
                bins=bins,
                ec="k",
                lw=1,
                label=label,
                density=True,
                alpha=0.25,
            )

            axarr[fi, 0].set_xlabel(f"{feat}")
            axarr[fi, 0].set_ylabel("Number of DNMs")
            axarr[fi, 0].legend()

            axarr[fi, 0].set_title(f"Distribution of {feat} values")

            sns.despine(ax=axarr[fi, 0], top=True, right=True)

        fpr, tpr, thresholds = roc_curve(y, X_feat)
        auc = round(roc_auc_score(y, X_feat), 3)

        # youden's j
        idx = np.argmax(tpr - fpr)
        if auc < 0.5: idx = np.argmax(fpr - tpr)
        threshold = thresholds[idx]

        # if feat == "kid AL diff.":
        # print (feat, thresholds)
        axarr[fi, 1].plot(fpr, tpr, c="firebrick", lw=2)
        axarr[fi, 1].axline(xy1=(0, 0), slope=1, c="k", ls=":")
        axarr[fi, 1].set_title(f"ROC curve using {feat} only\n(AUC = {auc}\nYouden's J at {threshold} for TPR of {round(tpr[idx], 2)} and FPR of {round(fpr[idx], 2)})")
        axarr[fi, 1].set_xlabel("FPR")
        axarr[fi, 1].set_ylabel("TPR")
        sns.despine(ax=axarr[fi, 1], top=True, right=True)

    f.tight_layout()
    f.savefig("hist.png", dpi=300)

    clf = RandomForestClassifier(max_depth=2, class_weight="balanced")
    clf.fit(X_train, y_train)
    importances = clf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    forest_importances = pd.Series(importances, index=FEATURES)

    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using random forest")
    ax.set_ylabel("Mean decrease in impurity")
    ax.set_xlabel("Feature name")
    sns.despine(ax=ax, top=True, right=True)
    fig.tight_layout()
    fig.savefig("importances.png", dpi=200)

    probs = clf.predict_proba(X_test)[:, 1]

    fpr, tpr, thresholds = roc_curve(y_test, probs)
    auc = round(roc_auc_score(y_test, probs), 3)

    f, ax = plt.subplots()
    ax.plot(fpr, tpr)
    ax.axline(xy1=(0, 0), slope=1, c="k", ls=":")
    sns.despine(ax=ax, top=True, right=True)
    ax.set_title(
        f"Binary classification of TRGT calls using\na class-balanced RF (AUC = {auc})"
    )
    ax.set_ylabel("TPR")
    ax.set_xlabel("FPR")
    f.tight_layout()
    f.savefig("roc.png", dpi=200)


if __name__ == "__main__":
    # p = argparse.ArgumentParser()
    # p.add_argument(
    #     "--mutations",
    #     type=str,
    #     help="""Path to CSV file containing estimated allele lengths.""",
    # )
    # p.add_argument(
    #     "--element_validation",
    #     type=str,
    #     help="""Path to output image.""",
    # )
    # p.add_argument(
    #     "--ont_validation",
    #     type=str,
    #     default="illumina",
    #     help="""Technology used to generate allele length estimates.""",
    # )
    # p.add_argument("--sample_id")

    # args = p.parse_args()
    # main(args)
    main()
