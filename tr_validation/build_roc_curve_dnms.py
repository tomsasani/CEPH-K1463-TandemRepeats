import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable, List
import numba
from schema import DeNovoSchema
import cyvcf2
from collections import Counter
import glob


from utils import filter_mutation_dataframe
from annotate_with_orthogonal_evidence import annotate_with_concordance

from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

plt.rc("font", size=12)


def annotate_with_allr(
    row: pd.Series,
    mom_vcf: cyvcf2.VCF,
    dad_vcf: cyvcf2.VCF,
):
    trid = row["trid"]
    chrom, start, end, _ = trid.split("_")
    region = f"{chrom}:{start}-{end}"

    dad_allr, mom_allr = "?", "?"
    for v in mom_vcf(region):
        if v.INFO.get("TRID") != trid: continue
        mom_allr = v.format("ALLR")[0]
    for v in dad_vcf(region):
        if v.INFO.get("TRID") != trid: continue
        dad_allr = v.format("ALLR")[0]
    return ";".join([dad_allr, mom_allr])

def check_for_overlap(row: pd.Series):

    allr = row["allr"]
    denovo_al = row["denovo_al"]
    

    is_overlap = False
    for parent in allr.split(";"):
        for allele_range in parent.split(","):
            a1, a2 = list(map(int, allele_range.split("-")))
            if a1 <= int(denovo_al) <= a2: is_overlap = True
    return int(is_overlap)

def get_dp(par: str):
    per_allele_reads = list(map(int, par.split(",")))
    return sum(per_allele_reads)


def main():

    ped = pd.read_csv(
        "tr_validation/data/file_mapping.csv",
        dtype={"paternal_id": str, "maternal_id": str, "sample_id": str},
    )

    # read in mutations
    ASSEMBLY = "GRCh38"
    mutations = []
    for fh in glob.glob(f"tr_validation/data/denovos/orthogonal_support/*.{ASSEMBLY}.element.read_support.csv"):
        df = pd.read_csv(
            fh,
            dtype={"sample_id": str},
        )
        res = []
        for i, row in df.iterrows():
            row_dict = row.to_dict()
            parental_overlap = annotate_with_concordance(row)
            row_dict.update({"validation_status": parental_overlap})
            res.append(row_dict)

        res = pd.DataFrame(res)
        mutations.append(res)

    mutations = pd.concat(mutations)
    print (mutations)
    
    # if we're looking at de novos, ensure that we filter
    # the DataFrame to exclude the non-denovo candidates
    # remove denovos where allele size is unchanged
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=False,
        remove_duplicates=False,
        remove_gp_ev=False,
        parental_overlap_frac_max=1,
        denovo_coverage_min=2,
    )

    print (mutations)
    # add parental depth columns
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )
    
    mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))
    mutations["has_transmission"] = mutations["sample_id"].apply(lambda s: s in ("2189", "2216"))

    mutations["father_overlap_coverage"] = mutations["father_overlap_coverage"].apply(lambda f: get_dp(f))
    mutations["mother_overlap_coverage"] = mutations["mother_overlap_coverage"].apply(lambda f: get_dp(f))
    mutations["denovo_al"] = mutations.apply(lambda row: int(row["child_AL"].split(",")[row["index"]]), axis=1)
    # mutations["is_transmitted"] = mutations["children_with_denovo_allele"].apply(lambda c: c != "unknown")
    mutations["is_validated"] = mutations["validation_status"].apply(lambda p: 1 if p == "pass" else 0)

    FEATURES = [
        "child_ratio",
        # "allele_ratio",
        "denovo_coverage",
        # "allele_coverage",
        # "child_coverage",
        # "mom_dp",
        # "dad_dp",
        "motif_size",
        "denovo_al",        
        # "parental_overlap_coverage_frac"
    ]

    CLASSIFICATION = "is_validated"

    if CLASSIFICATION == "is_validated":
        mutations = mutations[mutations["validation_status"] != "no_element_data"]

    candidates = mutations["trid"].nunique()
    print(f"{candidates} loci that are candidate DNMs")

    mutations = mutations.drop_duplicates(["trid", "genotype", "sample_id"])

    X = mutations[FEATURES]
    y = mutations[CLASSIFICATION].values

    f, ax = plt.subplots()
    sns.heatmap(X.corr(method="spearman"), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_title("Correlation between features")
    ax.set_yticks(np.arange(len(FEATURES)) + 0.5)
    ax.set_yticklabels(FEATURES, rotation=0)
    f.tight_layout()
    f.savefig("corr.png", dpi=200)

    X = X.values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, test_size=0.2
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
        X_feat = X[:, fi]

        fpr, tpr, thresholds = roc_curve(y, X_feat)
        auc = round(roc_auc_score(y, X_feat), 3)

        # youden's j
        idx = np.argmax(tpr - fpr)
        if auc < 0.5: idx = np.argmax(fpr - tpr)
        threshold = thresholds[idx]

        axarr[fi, 1].plot(fpr, tpr, c="firebrick", lw=2)
        axarr[fi, 1].axline(xy1=(0, 0), slope=1, c="k", ls=":")
        axarr[fi, 1].set_title(f"ROC curve using {feat} only\n(AUC = {auc}, Youden's J at {threshold})")
        axarr[fi, 1].set_xlabel("FPR")
        axarr[fi, 1].set_ylabel("TPR")
        sns.despine(ax=axarr[fi, 1], top=True, right=True)


        bins = np.linspace(np.min(X_feat), np.max(X_feat), num=50)

        for y_ in np.unique(y):
            y_i = np.where(y == y_)[0]
            X_feat_sub = X_feat[y_i]

            histvals, edges = np.histogram(X_feat_sub, bins=bins)
            histfracs = histvals / np.sum(histvals)
            cumsum = np.cumsum(histfracs)

            axarr[fi, 0].plot(
                bins[1:],
                cumsum,
                label=f"{y_} (n = {y_i.shape[0]})",
                lw=2,
                
            )

            # axarr[fi, 0].hist(
            #     X_feat_sub,
            #     bins=bins,
            #     ec="k",
            #     lw=1,
            #     label=f"{y_} (n = {y_i.shape[0]})",
            #     density=True,
            #     alpha=0.25,
            # )

            axarr[fi, 0].set_xlabel(f"{feat}")
            axarr[fi, 0].set_ylabel("Cumulative fraction of DNMs")
            axarr[fi, 0].legend(title=CLASSIFICATION)

            axarr[fi, 0].set_title(f"Distribution of {feat} values")

            sns.despine(ax=axarr[fi, 0], top=True, right=True)

        # subset to good ones

    f.tight_layout()
    f.savefig("G2.hist.png", dpi=300)

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
