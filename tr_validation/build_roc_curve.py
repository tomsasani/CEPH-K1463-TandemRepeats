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


N = 1_000_000


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
        # if li == 1:
        #     for key, val in v.items():
        #         print("\t".join([key, val]))

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
        if motif_len == 0: print (motifs)

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
                ],
            )
        ):
            X[li - 1, fi] = feat_val

        label = float(v["dist"])
        y[li - 1] = int(label == 0)


good_idxs = np.unique(np.where(np.all(~np.isnan(X), axis=1))[0])

X = X[good_idxs]
y = y[good_idxs].astype(np.int8)

# find TNs
tn_idxs = np.where(y == 0)[0]
tp_idxs = np.where(y == 1)[0]
print (tn_idxs.shape[0], tp_idxs.shape[0])
tp_idxs = np.random.choice(tp_idxs, size=tn_idxs.shape[0], replace=False)
idxs = np.concatenate((tn_idxs, tp_idxs))
print(tp_idxs.shape[0])
X = X[idxs]
y = y[idxs]

print (X[0])

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

print (X_train.shape, X_test.shape)


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

        if feat != "is homopolymer":
            axarr[fi, 0].plot(
                bins[1:], cumsum, label=label, lw=2,
            )
            if feat in ("kid min. AL", "kid AL diff."):
                axarr[fi, 0].set_xscale("log", base=10)
            axarr[fi, 0].set_xlabel(f"{feat}")
            axarr[fi, 0].set_ylabel("Cumulative proportion of calls")
            axarr[fi, 0].legend()

        else:
            axarr[fi, 0].bar([y_], np.sum(X_feat_sub) / X_feat_sub.shape[0], 1, label=label)
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

    #if feat == "kid AL diff.":
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
