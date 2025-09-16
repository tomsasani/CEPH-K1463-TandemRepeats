import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

df = pd.read_csv(
    "csv/hprc/combined.GRCh38.heterozygosity.tsv",
    sep="\t",
    names=[
        "chrom",
        "start",
        "end",
        "min_motiflen",
        "sample_id",
        "is_het",
        "is_denovo",
    ],
)
df = df[df["chrom"] != "chrX"]
df["motif"] = df["min_motiflen"].apply(lambda m: m if m < 7 else 7)
# get frac polymorphic at each TRID

df = df[df["n_motifs"] == 1]


# get number of HPRC samples with genotypes at each site
n_hprc = df.groupby(["chrom", "start", "end"]).agg(n_hprc=("sample_id", lambda s: len(set(s)))).reset_index()
df = df.merge(n_hprc)
# require all HPRC samples to be genotyped at each site
df = df.query("n_hprc == 100")
grouped = df.groupby(["chrom", "start", "end", "motif", "is_denovo"]).agg(n=("is_het", "sum")).reset_index()
grouped["n"] = grouped["n"] / 100

f, ax = plt.subplots(figsize=(8, 6))
sns.stripplot(data=grouped, x="motif", y="n", hue="is_denovo", dodge=True, ax=ax, alpha=0.5, palette="colorblind")
sns.boxplot(data=grouped, x="motif", y="n", hue="is_denovo", dodge=True, ax=ax, fliersize=0, palette="colorblind")
sns.despine(ax=ax)
ax.set_ylabel("Fraction of HPRC samples who are heterozygous")
ax.set_xlabel("Minimum motif size in TR locus")
ax.legend_.remove()
ax.set_xticks(range(7))
ax.set_xticklabels(list(map(str, range(1, 7))) + ["7+"])
f.savefig("hprc.png", dpi=200)
