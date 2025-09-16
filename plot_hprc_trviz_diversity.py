import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

plt.rc("font", size=11, )

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

TRID = "chr8_2623352_2623487_trsolve"

COHORT = "hprc"
GENOME = "GRCh38"
PALETTE = "colorblind"

MOTIF = "a"

key = pd.read_csv(f"trviz/{COHORT}/{TRID}.GRCh38.key.tsv", sep="\t", names=["sequence", "encoded_motif", "n_obs"])
seq = pd.read_csv(f"trviz/{COHORT}/{TRID}.GRCh38.encoded.tsv", sep="\t")

seq["allele_length"] = seq["encoded_motifs"].apply(lambda s: len(s))
seq["allele_length"] = seq["original_sequence"].apply(lambda s: len(s))
seq["contains_motif"] = seq["encoded_motifs"].str.contains(MOTIF)
seq["sample"] = seq["sample_id"].apply(lambda s: s.split("_")[0])
seq["hap"] = seq["sample_id"].apply(lambda s: int(s.split("_")[6]))

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

# get the motifs present in the CEPH individuals
key = pd.read_csv(f"trviz/ceph/{TRID}.{GENOME}.key.tsv", sep="\t", names=["motif", "encoded_motif", "count"])
encoded2motif = dict(zip(key["encoded_motif"], key["motif"]))

res_df = res_df[res_df["encoded_motif"].isin(key["encoded_motif"].unique())]

# identify the unique characters in the encoded
# motif arrays, including "-" gaps
uniq_chars = res_df["encoded_motif"].unique()#["a", "b", "c", "d", "e", "f", "g", "h", "i"]
uniq_vals = np.arange(len(uniq_chars))
cmap = sns.color_palette(PALETTE, uniq_vals.shape[0])

f, ax = plt.subplots(figsize=(5, 6))
sns.stripplot(
    data=res_df,
    orient="h",
    x="motif_count",
    y="encoded_motif",
    hue="encoded_motif",
    # marker_kws={"edgecolor": "w", "linewidth": 1},
    palette=PALETTE,
    edgecolor="w",
    linewidth=0.5,
    ax=ax,
)
ax.set_ylabel("Motif")
ax.set_xlabel("Number of motifs on HPRC haplotype")

# make custom legend
# get colors used in palette
custom_legend = []
for char, color in zip(uniq_chars, cmap):
    custom_legend.append(mpatches.Patch(color=color))


# ax.legend(
#     custom_legend,
#     key["motif"].to_list(),
#     prop={"family": "monospace", "size": 10},
#     frameon=True,
#     shadow=True,
#     #loc="lower left",
# )
sns.despine(ax=ax)
f.savefig("hprc_motifs.png", dpi=200)




f, ax = plt.subplots(figsize=(10, 5))

res_df = res_df[res_df["encoded_motif"].isin(["a", "c", "b"])]

max_count = res_df["motif_count"].max()

bins = np.arange(0, max_count + 1)

out = []
for motif, motif_df in res_df.groupby("encoded_motif"):
    histogram, edges = np.histogram(motif_df["motif_count"].values, bins=bins)
    for hi,h in enumerate(histogram):
        out.append({"motif": motif, "n_copies": hi, "count": h})
out_df = pd.DataFrame(out)

out_df = pd.DataFrame(
    {
        "n_copies": out_df["n_copies"].unique(),
        "c": out_df[out_df["motif"] == "c"]["count"].to_list(),
        "b": out_df[out_df["motif"] == "b"]["count"].to_list(),
        "a": out_df[out_df["motif"] == "a"]["count"].to_list(),
    }
)

out_df.to_csv("test.csv", sep=",", index=False)

cmap = sns.color_palette("colorblind", 6)
cmap = [cmap[1], cmap[2], cmap[4]]
cdict = dict(zip(("a", "b", "c"), cmap))
for motif in ("a", "b", "c"):
    motif_counts = res_df[res_df["encoded_motif"] == motif]
    ax.hist(motif_counts["motif_count"].values, bins=range(0, 100), label=motif, color=cdict[motif])
# ax.set_xlim(0, 8)
f.savefig("diversity.pruned.png", dpi=200)


sns.despine(ax=ax)

xp1, xp2 = 0.2, 0.8
diffs = (
    seq.groupby("sample")
    .agg(
        both_have_motif=("contains_motif", lambda x: len(set(x)) == 1),
        diff=("allele_length", lambda x: np.diff(x)[0]),
    )
    .reset_index()
)
diffs = diffs[diffs["both_have_motif"] == False]
max_diff = diffs["diff"].max()
min_diff = diffs["diff"].min()

f, ax = plt.subplots(figsize=(5, 5))
colors = sns.color_palette("colorblind", 2)
for sample, sample_df in seq.groupby("sample"):
    # if both haplotypes have the motif (or don't have it), skip
    if sample_df.shape[0] != 2: continue
    if len(sample_df["contains_motif"].unique()) == 1: continue
    # sort by hap
    lengths = sample_df.sort_values("hap")["allele_length"].values
    jitter = np.random.normal(0, 0.025)
    slope = lengths[1] - lengths[0]
    adj_slope = (slope - min_diff) / (max_diff - min_diff)
    if adj_slope == 0: adj_slope += 0.05
    ax.plot([xp1 + jitter, xp2 + jitter], lengths, c="darkgrey", alpha=adj_slope, zorder=0)
    ax.scatter([xp1 + jitter, xp2 + jitter], lengths,  ec="w", lw=0.5, zorder=1, color=colors)
ax.set_xticks([xp1, xp2])
ax.set_xticklabels(["No hyper-mutable motif", "Has hyper-mutable motif"])
ax.set_ylabel("Allele length (bp)")
ax.set_xlabel(f"Haplotype state of {COHORT.upper()} samples who are\nheterozygous for the hyper-mutable motif")
# ax.set_title("Haplotypes with the hyper-mutable motif are longer\nthan their homologs without hyper-mutable motifs")
sns.despine(ax=ax)
ax.set_xlim(0, 1)
f.tight_layout()
f.savefig("paired.png", dpi=200)

# permute haplotype status (whether it does or does not contains the hyper mutable motif) and emasure diff in allele lengths
allele_lengths = seq["allele_length"].values
motif_status = seq["contains_motif"].values

has_hypermutable_motif = np.where(motif_status == True)[0]
print (has_hypermutable_motif.shape)
median_allele_length = np.mean(allele_lengths[has_hypermutable_motif])
dist = []
for t in range(10_000):
    shuffled_status = np.random.permutation(motif_status)
    has_hypermutable_motif = np.where(shuffled_status == True)[0]
    median_allele_length_shuffled = np.mean(allele_lengths[has_hypermutable_motif])
    dist.append(median_allele_length_shuffled)
dist = np.array(dist)

pval = (np.sum(dist >= median_allele_length) + 1) / 10_000

colors = sns.color_palette("colorblind", 2)

all_dist = np.concatenate([dist, [median_allele_length]])
bins = np.arange(min(all_dist) * 0.9, max(all_dist) * 1.1, 2)
f, ax = plt.subplots(figsize=(6, 6))
ax.hist(dist, label="Permuted", ec="w", color=colors[0], lw=1, bins=bins)
ax.axvline(median_allele_length, label=f"Empirical (p = {pval})", ls=":", lw=2, c=colors[1])
ax.set_xlabel("Tandem repeat allele length")
ax.set_ylabel("Count")
ax.legend(
    title="Average size of TR alleles with HM motif",
    frameon=True,
    shadow=True,
    bbox_to_anchor=(0.425, 0.9)#, fontsize=11,
)
sns.despine(ax=ax)
f.tight_layout()
f.savefig("permuted.allele_length.png", dpi=200)

motif_status = seq["contains_motif"].values

f, ax = plt.subplots()
sns.stripplot(data=seq, x="contains_motif", y="allele_length", ax=ax)
ax.set_xlabel("Contains hyper-mutable motif?")
ax.set_ylabel("Allele length (# of motifs)")
f.savefig("o.png", dpi=200)
