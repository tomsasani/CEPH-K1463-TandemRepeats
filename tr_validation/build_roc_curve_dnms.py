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
from assess_concordance import assign_allele

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
    # read in raw DNM files
    STR_VCF_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/hiphase"

    MIN_DN_COVERAGE = 1
    MIN_TOTAL_COVERAGE = 1

    mutations: List[pd.DataFrame] = []

    ped = pd.read_csv(
        "tr_validation/data/file_mapping.csv",
        dtype={"paternal_id": str, "maternal_id": str, "sample_id": str},
    )
    SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
    SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
    SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))


    G4a = ["200081", "200082", "200084", "200085", "200086", "200087"]
    G4b = ["200101", "200102", "200103", "200104", "200105", "200106"]
    G2 = ["2209", "2188"]

    # loop over phased mutation dataframes
    for i, fh in enumerate(glob.glob("tr_validation/csv/phased/*.phased.2gen.tsv")):
        sample = fh.split("/")[-1].split(".")[0]
        #if sample not in G4a + G4b + G2: continue

        mdf = pd.read_csv(
            fh,
            sep="\t",
        )
        mdf = mdf[mdf["denovo_coverage"] >= MIN_DN_COVERAGE]
        mdf["phase"] = mdf["phase_summary"].apply(lambda p: p.split(":")[0])
        mdf["sample_id"] = sample

        # mom, dad = SMP2MOM[sample], SMP2DAD[sample]
        # mom_vcf = cyvcf2.VCF(
        #     f"{STR_VCF_PREF}/{mom}_{SMP2SUFF[mom]}_GRCh38_50bp_merge.sorted.phased.vcf.gz",
        #     gts012=True,
        # )
        # dad_vcf = cyvcf2.VCF(
        #     f"{STR_VCF_PREF}/{dad}_{SMP2SUFF[dad]}_GRCh38_50bp_merge.sorted.phased.vcf.gz",
        #     gts012=True,
        # )

        # mdf["allr"] = mdf.apply(lambda row: annotate_with_allr(row, mom_vcf, dad_vcf), axis=1)
        mutations.append(mdf)

    mutations = pd.concat(mutations)
    
    # if we're looking at de novos, ensure that we filter
    # the DataFrame to exclude the non-denovo candidates
    # remove denovos where allele size is unchanged
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=True,
        remove_duplicates=True,
    )

    print (mutations["allele_ratio"].min())

    mutations["generation"] = mutations["sample_id"].apply(lambda s: 1 if str(s).startswith("200") else 0)


    # add columns and filter
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )
    mutations["maternal_to_paternal_ratio"] = mutations["mom_dp"] / mutations["dad_dp"]

    mutations = mutations[mutations["maternal_to_paternal_ratio"].between(0.75, 1.25)]

    #mutations["denovo_overlap"] = mutations.apply(lambda row: check_for_overlap(row), axis=1)
    mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

    mutations["phase"] = mutations["phase_summary"].apply(lambda p: 1 if p.split(":")[0] == "dad" else 0)

    # mutations = mutations[
    #     (mutations["child_coverage"] >= MIN_TOTAL_COVERAGE)
    #     & (mutations["mom_dp"] >= MIN_TOTAL_COVERAGE)
    #     & (mutations["dad_dp"] >= MIN_TOTAL_COVERAGE)
    # ]

    # print (mutations["allele_ratio"].min())


    # combine with element validation data
    validation_df: List[pd.DataFrame] = []
    for fh in glob.glob("tr_validation/csv/*.element.dnm.read_support.csv"):
        df = pd.read_csv(fh)
        validation_df.append(df)

    if len(validation_df) > 0:
        validation_df = pd.concat(validation_df)

        mutations = mutations.merge(validation_df)

        GROUP_COLS = [
            "trid",
            "sample_id",
            "genotype",
            "exp_allele_diff_denovo",
            "exp_allele_diff_non_denovo",
        ]
        mutations = mutations[mutations["sample"] == "kid"]
        mutation_totals = (
            mutations.groupby(GROUP_COLS)
            .agg({"Read support": "sum"})
            .reset_index()
            .rename(columns={"Read support": "Total read support"})
        )
        mutations = mutations.merge(mutation_totals, on=GROUP_COLS, how="left")
        mutations["Read support frac"] = mutations["Read support"] / mutations["Total read support"]

        # require at least 10 total spanning reads and at least 2 reads supporting a given diff
        mutations = mutations[
            (mutations["Total read support"] >= 10)
            & (mutations["Total read support"] <= 100)
        ]

        # for every observed number of net CIGAR operations in a read,
        # figure out if that net CIGAR is closer to the expected net CIGAR
        # operation for a REF or ALT allele.
        mutations["allele_assignment"] = mutations.apply(lambda row: assign_allele(row), axis=1)
        mutation_total_support = mutations.groupby(["trid", "allele_assignment"]).agg({"Read support": "sum"}).reset_index().rename(columns={"Read support": "Total read support for allele"})
        mutations = mutations.merge(mutation_total_support, how="left").fillna(value=0)

        res: List[pd.DataFrame] = []
        for trid, trid_df in mutations.groupby(["trid", "genotype"]):
            if "denovo" not in trid_df["allele_assignment"].to_list():
                is_concordant = False
            else:
                # assumption is at least 1 read supports the denovo
                is_concordant = True
            trid_df["is_concordant"] = is_concordant
            trid_df = trid_df.drop_duplicates(["trid", "genotype"])
            res.append(trid_df)
        mutations = pd.concat(res)
        
        #mutations = mutations[mutations["allele_assignment"] == "denovo"]

    mutations["father_overlap_coverage"] = mutations["father_overlap_coverage"].apply(lambda f: get_dp(f))
    mutations["mother_overlap_coverage"] = mutations["mother_overlap_coverage"].apply(lambda f: get_dp(f))

    candidates = mutations["trid"].nunique()
    print(f"{candidates} loci that are candidate DNMs")

    # add a transmission column if it's not in there
    # require that the samples have children
    # mutations["expected_children"] = mutations["sample_id"].apply(
    #     lambda s: 6 if s in ("2189", "2216") else 8 if s in ("2209", "2188") else np.nan
    # )
    # mutations = mutations[mutations["n_children"] == mutations["expected_children"]]

    # MIN_TRANSMISSIONS = 1

    # mutations["Is transmitted?"] = mutations["n_with_denovo_allele"].apply(
    #     lambda n: 1 if n >= MIN_TRANSMISSIONS else 0
    # )

    non_inf = np.where(~np.isinf(mutations["allele_ratio"].values))[0]
    mutations = mutations.iloc[non_inf]
    print (mutations["allele_ratio"].min())

    FEATURES = [
        "child_ratio",
        "allele_ratio",
        "denovo_coverage",
        "allele_coverage",
        "child_coverage",
        "mom_dp",
        "dad_dp",
        "motif_size",
        "phase",
        "denovo_al",
        # "is_concordant",
        "phase",
        
    ]

    CLASSIFICATION = "is_concordant"

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

        if feat == "allele_ratio":
            print (np.min(X_feat))

        # y_clean_idxs = np.where(y)[0]
        # y_clean = y[y_clean_idxs]
        # y_clean[y_clean > 1] = 1
        fpr, tpr, thresholds = roc_curve(y, X_feat)
        auc = round(roc_auc_score(y, X_feat), 3)

        # youden's j
        idx = np.argmax(tpr - fpr)
        if auc < 0.5: idx = np.argmax(fpr - tpr)
        threshold = thresholds[idx]

        axarr[fi, 1].plot(fpr, tpr, c="firebrick", lw=2)
        axarr[fi, 1].axline(xy1=(0, 0), slope=1, c="k", ls=":")
        axarr[fi, 1].set_title(f"ROC curve using {feat} only\n(AUC = {auc}, Youden's J at {threshold})")# for TPR of {round(tpr[idx], 2)} and FPR of {round(fpr[idx], 2)})")
        axarr[fi, 1].set_xlabel("FPR")
        axarr[fi, 1].set_ylabel("TPR")
        sns.despine(ax=axarr[fi, 1], top=True, right=True)

        # infs
        # inf_idxs = np.where(np.isinf(X_feat))
        # non_inf_idxs = np.where(~np.isinf(X_feat))
        # X_feat = X_feat[non_inf_idxs]
        # y_prime = y[non_inf_idxs]

        bins = np.linspace(np.min(X_feat), np.max(X_feat), num=50)
        # bins = np.arange(np.min(X_feat), np.max(X_feat), step=1 if np.max(X_feat) > 1 else 0.05)

        for y_ in np.unique(y):
            y_i = np.where(y == y_)[0]
            X_feat_sub = X_feat[y_i]

            histvals, edges = np.histogram(X_feat_sub, bins=bins)
            histfracs = histvals / np.sum(histvals)
            cumsum = np.cumsum(histfracs)

            #label = "Maternal" if y_ == 0 else "Paternal"

            # axarr[fi, 0].plot(
            #     bins[1:],
            #     cumsum,
            #     label=f"{y_} (n = {y_i.shape[0]})",
            #     lw=2,
                
            # )

            axarr[fi, 0].hist(
                X_feat_sub,
                bins=bins,
                ec="k",
                lw=1,
                label=f"{y_} (n = {y_i.shape[0]})",
                density=True,
                alpha=0.25,
            )

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
