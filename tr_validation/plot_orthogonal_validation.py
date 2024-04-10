import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import utils
from sklearn.mixture import GaussianMixture

df = pd.read_csv("tr_validation/csv/orthogonal_support/2217.ont.read_support.csv").dropna()
df["sample_id"] = "2217"
df = utils.filter_mutation_dataframe(
    df,
    remove_duplicates=False,
    remove_gp_ev=True,
    remove_complex=True,
)

# print (df.groupby(["sample", "trid"]).size().sort_values())

trid = np.random.choice(df["trid"].unique())
#trid = "chr17_80963375_80963824_trsolve"
#trid = "chr9_40311198_40311211_TRF"
trid = "chr9_68323850_68324860_trsolve"
print (trid)
df_sub = df[df["trid"] == trid]
print (df_sub["region"].unique()[0])
print (df_sub)
df_sub["denovo_al"] = df_sub.apply(lambda row: row["child_AL"].split(",")[row["index"]], axis=1)
df_sub["non_denovo_al"] = df_sub.apply(lambda row: row["child_AL"].split(",")[1 - row["index"]], axis=1)

new = []
for _, row in df_sub.iterrows():
    n_reads = row["Read support"]
    for i in range(n_reads):
        new.append(row.to_dict())

new = pd.DataFrame(new)
# fit gaussian mixture to child's reads
X = new[new["sample"] == "kid"]["diff"].values.reshape(-1, 1)
gm = GaussianMixture(n_components=2).fit(X)
mixture_means = gm.means_ 
mixture_cov = gm.covariances_
mixture_std = np.array([np.sqrt(np.trace(mixture_cov[i]) / 2) for i in range(2)])
# figure out if either parent has a read that overlaps a normal distribution
# centered around the child's means + stdevs

exp_denovo_diff = df_sub["exp_allele_diff_denovo"].unique()[0]
exp_non_denovo_diff = df_sub["exp_allele_diff_non_denovo"].unique()[0]
exp_denovo = df_sub["denovo_al"].unique()[0]
exp_non_denovo = df_sub["non_denovo_al"].unique()[0]

f, ax = plt.subplots()
sns.stripplot(data=new, x="sample", y="diff", ax=ax)
ax.axhline(y=exp_denovo_diff, c="firebrick", label=f"De novo (AL = {exp_denovo})")
ax.axhline(y=exp_non_denovo_diff, c="dodgerblue", label=f"Non de novo (AL = {exp_non_denovo})")


# keep track of how many alleles have overlap
allele_overlap = np.zeros((2, 2))
for allele_i in range(2):
    mean, std = mixture_means[allele_i], 2 * mixture_std[allele_i]
    lo, hi = mean - std, mean + std
    ax.fill_between(x=[-0.5, 2.5], y1=lo, y2=hi, color="gainsboro", alpha=0.5)
    ax.axhline(y=mean, c="k", ls=":")
    for parent_i, parent in enumerate(("mom", "dad")):
        p_diffs = df_sub[df_sub["sample"] == parent]["diff"].values
        p_reps = df_sub[df_sub["sample"] == parent]["Read support"].values
        print (allele_i, parent, p_reps)
        X_p = np.repeat(p_diffs, p_reps)
        in_parent = np.sum((X_p >= lo) & (X_p <= hi))
        allele_overlap[allele_i, parent_i] += in_parent / np.sum(p_reps) 

print (allele_overlap)

ax.set_ylabel("CIGAR operations w/r/t reference")
ax.set_xlabel("Sample in trio")
ax.legend(frameon=False, title="Expected CIGAR operations w/r/t reference")
ax.set_title(trid)
sns.despine(ax=ax)
f.savefig("o.png", dpi=200)
