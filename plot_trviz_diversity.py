import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

TRID = "chr1_54510591_54510710_trsolve"
TRID = "chr8_2623352_2623487_trsolve"

COHORT = "ceph"

key = pd.read_csv(f"tr_validation/trviz/{COHORT}/{TRID}.GRCh38.key.tsv", sep="\t", names=["sequence", "encoded_motif", "n_obs"])
seq = pd.read_csv(f"tr_validation/trviz/{COHORT}/{TRID}.GRCh38.encoded.tsv", sep="\t")

seq["allele_length"] = seq["encoded_motifs"].apply(lambda s: len(s))
seq["contains_motif"] = seq["encoded_motifs"].str.contains("a")
seq["sample"] = seq["sample_id"].apply(lambda s: s.split("_")[0])
seq["hap"] = seq["sample_id"].apply(lambda s: int(s.split("_")[-1]))

# for each "encoded motif" in the key...
res = []
for encoded_motif in key["encoded_motif"].unique():
    # figure out how many copies of that encoded motif every sample had
    motif_counts_per_sample = []
    for i, row in seq.iterrows():
        sample_id = row["sample_id"]
        sample_encoded_motifs = row["encoded_motifs"]
        motif_count = sum([m == encoded_motif for m in list(sample_encoded_motifs)])
        res.append(
            {
                "sample_id": sample_id,
                "encoded_motif": encoded_motif,
                "motif_count": motif_count,
            }
        )
res_df = pd.DataFrame(res)

f, ax = plt.subplots()
sns.stripplot(data=res_df, x="encoded_motif", y="motif_count", hue="encoded_motif", ax=ax)
ax.set_xlabel("Motif label")
ax.set_ylabel("Number of motifs on haplotype")
sns.despine(ax=ax)
f.savefig("diversity.png", dpi=200)

f, ax = plt.subplots()
for sample, sample_df in seq.groupby("sample"):
    # if both haplotypes have the motif (or don't have it), skip
    if sample_df.shape[0] != 2: continue
    if len(sample_df["contains_motif"].unique()) == 1: continue
    # sort by hap
    lengths = sample_df.sort_values("hap")["allele_length"].values
    ax.plot([0, 1], lengths, c="gainsboro", alpha=0.5)
    ax.scatter([0, 1], lengths, c="k")
ax.set_xticks([0, 1])
ax.set_xticklabels(["No hyper-mutable motif", "Has hyper-mutable motif"])
ax.set_ylabel("Allele length")
ax.set_xlabel("HPRC haplotypes")
f.tight_layout()
f.savefig("paired.png", dpi=200)

# permute haplotype status (whether it does or does not contains the hyper mutable motif) and emasure diff in allele lengths

allele_lengths = seq["allele_length"].values
motif_status = seq["contains_motif"].values

is_mutable = np.where(motif_status == True)

median_allele_length = np.median(allele_lengths[is_mutable])
dist = []
for t in range(1_000):
    shuffled_status = np.random.permutation(motif_status)
    is_mutable = np.where(shuffled_status == True)
    median_allele_length_shuffled = np.median(allele_lengths[is_mutable])
    dist.append(median_allele_length_shuffled)
dist = np.array(dist)

pval = (np.sum(dist >= median_allele_length) + 1) / 1_000
f, ax = plt.subplots()
ax.hist(dist, label="Permuted", ec="w", lw=1)
ax.axvline(median_allele_length, label=f"Empirical (p = {pval})", ls=":", c="firebrick")
ax.legend(title = "Median allele length of haplotypes\nwith hyper-mutable motif")
f.tight_layout()
f.savefig("t.png", dpi=200)

motif_status = seq["contains_motif"].values

f, ax = plt.subplots()
sns.stripplot(data=seq, x="contains_motif", y="allele_length", ax=ax)
ax.set_xlabel("Contains hyper-mutable motif?")
ax.set_ylabel("Allele length (# of motifs)")
f.savefig("o.png", dpi=200)
