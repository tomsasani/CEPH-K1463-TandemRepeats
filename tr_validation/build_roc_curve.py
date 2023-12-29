import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Callable
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import tqdm
import csv
import cyvcf2
import pysam
from validate_dnms import get_read_diff
from collections import Counter
import scipy.stats as ss


def query_vcf(
    vcf: cyvcf2.VCF,
    region: str,
    trid: str,
    ref_al: int,
    alt_al: int,
):
    exp_diff_ref, exp_diff_alt = None, None

    for var in vcf(region):
        var_trid = var.INFO.get("TRID")
        if trid != var_trid:
            continue
        # get information about a variant from a VCF
        exp_diff_alt = alt_al - len(var.REF)
        exp_diff_ref = ref_al - len(var.REF)

    return exp_diff_ref, exp_diff_alt


def apply_op(value: str, operation: Callable):
    try:
        return operation(list(map(float, value.split(","))))
    except AttributeError:
        return np.nan


def find_length(value: str):
    return len(value)


def calculate_ab(sd: str):
    try:
        a_sd, b_sd = list(map(float, sd.split(",")))
        return a_sd / sum([a_sd, b_sd])
    except ValueError:
        return np.nan


plt.rc("font", size=12)


N = 100_000


FEATURES = [
    "kid AB",
    "mom AB",
    "dad AB",
    "kid total SD",
    "mom total SD",
    "dad total SD",
    # "kid min. AL",
    "kid AL diff.",
    "motif length",
    "is homopolymer",
    # "element spanning reads",
    # "element AB diff",
    # "pass element",
]


X = np.zeros((N, len(FEATURES)), dtype=np.float32)
y = np.zeros(N)

HEADER = None

# AP = AL * AP (i.e., # of base pairs that match, rather than fraction)
# just divide by AL to get back to fraction
# haven't been using AM, MS, ALLR
# could use ALLR as metric, how wide the interval is?
# MS is motif span per allele
# AM is methylation
# parent_ht is the set of alleles that produced the specified manhattan distance in `dist`
# try with just motifs of a certain size?

# 1, 2, 3, 4, 5, 6, 7-10, 11-15, 16-20, 21-30, 30+, compound with 2, compound with 3+
# ignore complex motifs?
# come up with a metric that characterizes how much the spanning reads look good?
# convert SD into adjusted allele balance, minimum allele balance, binomial p-value, total depth
# reference AL, diff in ALs
# include kid genotype

# TODO: add in element validation here!

KID, MOM, DAD = "NA12879", "NA12878", "NA12877"

VCF_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/GRCh38/dashnow/fd6baac/"
BAM_PREF = (
    "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/complete/"
)

KID_VCF = cyvcf2.VCF(VCF_PREF + "2216_DM_GRCh38.sorted.vcf.gz", gts012=True)
MOM_VCF = cyvcf2.VCF(VCF_PREF + "2188_DM_GRCh38.sorted.vcf.gz", gts012=True)
DAD_VCF = cyvcf2.VCF(VCF_PREF + "2209_SF_GRCh38.sorted.vcf.gz", gts012=True)

KID_BAM = pysam.AlignmentFile(BAM_PREF + "2216_1_merged_sort.bam", "rb")
MOM_BAM = pysam.AlignmentFile(BAM_PREF + "2188_1_merged_sort.bam", "rb")
DAD_BAM = pysam.AlignmentFile(BAM_PREF + "2209_2_merged_sort.bam", "rb")

SMP2VCF = dict(zip([KID, MOM, DAD], [KID_VCF, MOM_VCF, DAD_VCF]))
SMP2BAM = dict(zip([KID, MOM, DAD], [KID_BAM, MOM_BAM, DAD_BAM]))

# NOTE: probably don't want to double count observations?

STRICT = True

with open(
    "/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/GRCh38/dashnow/fd6baac/all.dists.sorted.txt",
    "r",
) as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for li, l in tqdm.tqdm(enumerate(csvf)):
        if li > N:
            break
        if li == 0:
            HEADER = l
            continue
        v = dict(zip(HEADER, l))

        if v["kid_id"] not in SMP2VCF:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue
        if "pathogenic" in v["#variant"]:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue

        trid = v["#variant"]
        chrom, start, end, _ = trid.split("_")
        region = f"{chrom}:{start}-{end}"

        try:
            ref_al, alt_al = list(map(int, v["kid_AL"].split(",")))
        except ValueError:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue

        sample_vcf = SMP2VCF[v["kid_id"]]

        # use Element data to figure out allele balances in the focal sample
        # first, get the expected read diffs to the reference

        exp_allele_diff_ref, exp_allele_diff_alt = query_vcf(
            sample_vcf,
            region,
            trid,
            ref_al,
            alt_al,
        )

        sample_bam = SMP2BAM[v["kid_id"]]

        # then, query the BAMs to see if the reads support the expectation
        diffs = []
        for read in sample_bam.fetch(chrom, int(start), int(end)):
            diff = get_read_diff(
                read,
                int(start),
                int(end),
                slop=max([abs(exp_allele_diff_ref), abs(exp_allele_diff_alt)]),
            )
            diffs.append(diff)

        diff_counts = Counter(diffs).most_common()

        # calculate element allele balance

        # make sure we only consider sites that *Could* be validated by
        # element
        spanning_reads = sum([v for k, v in diff_counts if k is not None])

        if spanning_reads < 10:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue

        match_ref, match_alt = 0, 0

        match_ref = sum([v for k, v in diff_counts if k == exp_allele_diff_ref])
        match_alt = sum([v for k, v in diff_counts if k == exp_allele_diff_alt])

        ref_ab = match_ref / spanning_reads if spanning_reads > 0 else 0.
        alt_ab = match_alt / spanning_reads if spanning_reads > 0 else 0.

        # calculate "degree of support" for the specified call by
        # performing a binomial test on the allele balance and asking
        # if the p-value is < 0.05
        element_ab_diff = 0

        if exp_allele_diff_alt == exp_allele_diff_ref:
            expected_p, alternative = 1.0, "less"
            element_ab_diff += abs(ref_ab - expected_p)
        else:
            expected_p, alternative = 0.5, "two-sided"
            element_ab_diff += abs(ref_ab - expected_p)
            element_ab_diff += abs(alt_ab - expected_p)

        #print (exp_allele_diff_ref, exp_allele_diff_alt, ref_ab, alt_ab, v["dist"], element_ab_diff)

        # NOTE: confounder!!! if mendelian consistent calls are mostly HOM
        # and violations are mostly HET, comparing allele balance diffs is
        # not appropriate

        p = ss.binom_test(
            match_ref,
            spanning_reads,
            p=expected_p,
            alternative=alternative,
        )

        pass_element_p = p > 0.05

        kid_ab, mom_ab, dad_ab = [
            calculate_ab(v[k]) for k in ("kid_SD", "mom_SD", "dad_SD")
        ]
        kid_al_min = apply_op(v["kid_AL"], np.min)
        kid_al_max = apply_op(v["kid_AL"], np.max)

        # add in AL max-min for both parents as well
        kid_al_diff = abs(kid_al_max - kid_al_min)

        motifs = v["motifs"].split(",")
        if len(motifs) > 1:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue

        motif_lens = [len(m) for m in motifs]
        motif_len = min(motif_lens)
        if motif_len != 1:
            X[li - 1, :] = np.nan
            y[li - 1] = np.nan
            continue

        kid_dp, mom_dp, dad_dp = [
            apply_op(v[k], np.sum) for k in ("kid_SD", "mom_SD", "dad_SD")
        ]

        is_homopolymer = len(set(list(motifs[0]))) == 1

        for fi, (feat, feat_val) in enumerate(
            zip(
                FEATURES,
                [
                    kid_ab,
                    mom_ab,
                    dad_ab,
                    kid_dp,
                    mom_dp,
                    dad_dp,
                    # kid_al_min,
                    kid_al_diff,
                    motif_len,
                    is_homopolymer,
                    # spanning_reads,
                    # element_ab_diff,
                    # pass_element_p,
                ],
            )
        ):
            X[li - 1, fi] = feat_val


        # if strict, require both a dist of 0 and a passing element AB
        label = float(v["dist"])

        class_label = None
        if label == 0 and pass_element_p:
            class_label = 1
        elif label == 0 and (not pass_element_p):
            class_label = None
        else:
            class_label = 0
        
        if class_label is None:
            y[li - 1] = np.nan
        else:
            y[li - 1] = class_label
        


good_idxs = np.unique(np.where(np.all(~np.isnan(X), axis=1))[0])
good_idxs = np.intersect1d(good_idxs, np.where(~np.isnan(y)))

X = X[good_idxs]
y = y[good_idxs].astype(np.int8)

# find TNs
tn_idxs = np.where(y == 0)[0]
tp_idxs = np.where(y == 1)[0]
print(tn_idxs.shape[0], tp_idxs.shape[0])
# tp_idxs = np.random.choice(tp_idxs, size=tn_idxs.shape[0], replace=False)
idxs = np.concatenate((tn_idxs, tp_idxs))
print(tp_idxs.shape[0])
X = X[idxs]
y = y[idxs]

print(X[0])

print(np.sum(y) / y.shape[0])

X_df = pd.DataFrame(X, columns=FEATURES)
f, ax = plt.subplots()
sns.heatmap(X_df.corr(method="spearman"), cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_title("Correlation between features")
ax.set_yticks(np.arange(len(FEATURES)) + 0.5)
ax.set_yticklabels(FEATURES, rotation=0)
f.tight_layout()
f.savefig("corr.png", dpi=200)

X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, test_size=0.2)

print(X_train.shape, X_test.shape)


f, axarr = plt.subplots(
    len(FEATURES), 2, figsize=(12, (4 * len(FEATURES))), sharex=False, sharey=False
)

for fi, feat in enumerate(FEATURES):
    # get feature values for this feature
    X_feat = X[:, fi]

    # # for some features, only plot 99% of the data that's
    # # not a huge outlier
    # if feat in ("kid_AL_min", "kid_AL_diff", "motif_len"):
    #     pctile = np.percentile(X_feat, q=99)
    #     print (pctile, feat)
    #     idxs = np.where(X_feat <= pctile)[0]

    bins = np.linspace(np.min(X_feat), np.max(X_feat), num=100)

    for y_ in np.unique(y):
        y_i = np.where(y == y_)[0]
        X_feat_sub = X_feat[y_i]

        # print (feat, y_, np.min(X_feat_sub), np.max(X_feat_sub))

        # histogram
        histvals, edges = np.histogram(X_feat_sub, bins=bins)
        histfracs = histvals / np.sum(histvals)
        cumsum = np.cumsum(histfracs)

        label = "Consistent" if y_ == 1 else "Mendel. violation"

        if feat not in ("is homopolymer", "pass element"):
            axarr[fi, 0].plot(
                bins[1:],
                cumsum,
                label=label,
                lw=2,
            )
            # if feat in ("kid min. AL", "kid AL diff."):
            #     axarr[fi, 0].set_xscale("log", base=10)
            axarr[fi, 0].set_xlabel(f"{feat}")
            axarr[fi, 0].set_ylabel("Cumulative proportion of calls")
            axarr[fi, 0].legend()

        else:
            axarr[fi, 0].bar(
                [y_], np.sum(X_feat_sub) / X_feat_sub.shape[0], 1, label=label
            )
            axarr[fi, 0].set_xticks([0, 1])
            axarr[fi, 0].set_xticklabels(["Mendelian viol.", "Consistent"])
            axarr[fi, 0].set_ylabel("Proportion of calls")

        # axarr[fi, 0].hist(
        #     X_feat_sub,
        #     histtype="step",
        #     density=True,
        #     bins=bins,
        #     label="Consistent" if y_ == 1 else "Mendel. violation",
        #     lw=2,
        # )
        axarr[fi, 0].set_title(f"Distribution of {feat} values in TRGT GRCh38 calls")

        sns.despine(ax=axarr[fi, 0], top=True, right=True)

    # # ROC curve
    if feat in ("kid AL diff.", "is homopolymer"):
        fpr, tpr, thresholds = roc_curve(1 - y, X_feat)
        auc = round(roc_auc_score(1 - y, X_feat), 3)
    else:
        fpr, tpr, thresholds = roc_curve(y, X_feat)
        auc = round(roc_auc_score(y, X_feat), 3)

    # if feat == "kid AL diff.":
    # print (feat, thresholds)
    axarr[fi, 1].plot(fpr, tpr, c="firebrick", lw=2)
    axarr[fi, 1].axline(xy1=(0, 0), slope=1, c="k", ls=":")
    axarr[fi, 1].set_title(f"ROC curve using {feat} only\n(AUC = {auc})")
    axarr[fi, 1].set_xlabel("FPR")
    axarr[fi, 1].set_ylabel("TPR")
    sns.despine(ax=axarr[fi, 1], top=True, right=True)

f.tight_layout()
f.savefig("hist.png", dpi=200)

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
# res_df = pd.concat(res_df)

# f, ax = plt.subplots()
# sns.histplot(
#     data=res_df,
#     x="min_SD",
#     hue="TP",
#     element="step",
#     stat="percent",
#     common_norm=False,
#     cumulative=True,
#     fill=False,
#     ax=ax,
# )
# sns.despine(ax=ax, top=True, right=True)
# f.tight_layout()
# f.savefig("dist.png", dpi=200)
