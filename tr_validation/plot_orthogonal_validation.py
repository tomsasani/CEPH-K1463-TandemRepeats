import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import utils
from sklearn.mixture import GaussianMixture

ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

df = pd.read_csv(f"tr_validation/data/denovos/orthogonal_support/2189.{ASSEMBLY}.element.read_support.csv")
df = utils.filter_mutation_dataframe(
    df,
    remove_duplicates=False,
    remove_gp_ev=False,
    remove_complex=True,
    denovo_coverage_min=2,
    child_ratio_min=0.2,
)

df = df[df["motif"].str.len() == 1]

trid = np.random.choice(df["trid"].unique())
print (trid)
df_sub = df[df["trid"] == trid]

df_sub["denovo_al"] = df_sub.apply(lambda row: row["child_AL"].split(",")[row["index"]], axis=1)
df_sub["non_denovo_al"] = df_sub.apply(lambda row: row["child_AL"].split(",")[1 - row["index"]], axis=1)

print (df_sub[["exp_allele_diff_denovo", "exp_allele_diff_non_denovo", "kid_evidence", "dad_evidence", "mom_evidence"]])
assert df_sub.shape[0] == 1

# gather diffs in all members of the trio
for i,row in df_sub.iterrows():
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


# fit gaussian mixture to child's reads
X = np.array(kid_diffs)
gm = GaussianMixture(n_components=2).fit(X.reshape(-1, 1))
# get the means and covariance sof the two components
mixture_means = gm.means_ 
mixture_cov = gm.covariances_
# calculate standard deviation from covariance matrix
# standard deviation (sigma) is sqrt of the product of the residuals with themselves
# (i.e., sqrt of the mean of the trace of the covariance matrix)
mixture_std = np.array(
    [np.sqrt(np.trace(mixture_cov[i]) / 2) for i in range(2)]
)
# figure out if either parent has a read that overlaps a normal distribution
# centered around the child's means + stdevs

exp_denovo_diff = df_sub["exp_allele_diff_denovo"].unique()[0]
exp_non_denovo_diff = df_sub["exp_allele_diff_non_denovo"].unique()[0]
exp_denovo = df_sub["denovo_al"].unique()[0]
exp_non_denovo = df_sub["non_denovo_al"].unique()[0]

f, ax = plt.subplots()

sns.stripplot([kid_diffs, mom_diffs, dad_diffs], palette="colorblind", alpha=0.5)
ax.axhline(exp_denovo_diff, ls=":", c="dodgerblue", label="Expected allele length (de novo)")
ax.axhline(exp_non_denovo_diff, ls=":", c="firebrick", label="Expected allele length (non de novo)")
ax.set_ylabel("Allele length inferred using Element reads")
ax.set_xlabel("Sample")
ax.set_xticks(range(3))
ax.set_xticklabels(["Kid", "Mom", "Dad"])
ax.legend()
sns.despine(ax=ax)
f.savefig('o.png', dpi=200)