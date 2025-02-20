from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TRID = "chr8_2623352_2623487_trsolve"
GENOME = "GRCh38"

# store unique characters
smps = []

hapseq = []
for record in SeqIO.parse(f"tr_validation/trviz/ceph/{TRID}.{GENOME}_alignment_output.fa", "fasta"):
    hapseq.append(list(record.seq))
    smps.append(record.id)

hapseq = pd.DataFrame(hapseq)

orig_cols = hapseq.columns

# get vectorized mapping to integers
uniq_chars = np.unique(hapseq.values)
uniq_vals = np.arange(uniq_chars.shape[0])

cmap = sns.color_palette("colorblind", uniq_chars.shape[0] - 1)

mapping_dict = dict(zip(uniq_chars, uniq_vals))
vectorized_map = np.vectorize(lambda x: mapping_dict[x]) 

hapseq["sample"] = smps
hapseq["sample_id"] = hapseq["sample"].apply(lambda s: s.split("_")[0])

hapseq["generation"] = hapseq["sample_id"].apply(
    lambda s: (
        "G2A"
        if s == "2209"
        else (
            "G2B"
            if s == "2188"
            else (
                "G4A"
                if (s.startswith("2000") and s != "200080")
                else (
                    "G4B"
                    if (s.startswith("2001") and s != "200100")
                    else "G1" if s in ("2281", "2280", "2213", "2214") else "G3"
                )
            )
        )
    )
)


hapseq["parent_status"] = hapseq["sample_id"].apply(lambda s: (
            "G4A"
            if s in ("2216", "200080") else "G4B" if s in ("2189", "200100")
            else (
                "G3"
                if s in ("2209", "2188")
                else "G2A" if s in ("2280", "2281") else "G2B" if s in ("2214", "2213") else "UNK"
            )
        ))

hapseq["has_denovo"] = hapseq["sample"].apply(lambda s: s.split("_")[-1])

# don't plot G1 as a separate generation
gen_to_plot = [g for g in hapseq["generation"].unique() if g != "G1"]

height_ratios = [1]
n_per_gen = hapseq.groupby("generation").size().to_dict()

min_per_gen = min([v for k,v in n_per_gen.items()])

height_ratios = [v for k,v in n_per_gen.items()]
height_ratios = [1, 0.75, 0.75, 3, 2.5, 2.5]

f, axarr = plt.subplots(
    hapseq["generation"].nunique(),
    sharex=False,
    figsize=(8, 18),
    height_ratios=height_ratios,
)

# loop over the generations
for gen_i, (gen, _) in enumerate(n_per_gen.items()):

    # collect the kids in this generation
    kid_df = hapseq[hapseq["generation"] == gen]
    # plot the kids second
    kid_df["order"] = 1
    # collect the parents in this generation
    par_df = hapseq[hapseq["parent_status"] == gen]
    # plot the parents first
    par_df["order"] = 0

    gen_df = pd.concat([par_df, kid_df]).sort_values("order")

    # get a list of haplotypes with DNMs
    has_denovo = gen_df["has_denovo"].values == "denovo"

    gen_hapseq = gen_df[orig_cols].values
    gen_hapseq = vectorized_map(gen_hapseq)

    gen_hapseq = np.where(gen_hapseq == 0, np.nan, gen_hapseq)

    sns.heatmap(gen_hapseq, ax=axarr[gen_i], lw=0.5, ec="w", cbar=False, cmap=cmap)
    axarr[gen_i].set_yticks(np.arange(gen_df.shape[0]) + 0.5)
    axarr[gen_i].set_yticklabels(gen_df["sample_id"].values, rotation=0)
    axarr[gen_i].set_title(gen)

    for i, t in enumerate(axarr[gen_i].yaxis.get_ticklabels()):
        if has_denovo[i] and t.get_text() not in par_df["sample_id"].to_list():
            t.set_color("red")
            t.set_fontstyle("italic")
        if t.get_text() in par_df["sample_id"].to_list():
            t.set_fontweight("bold")


    # add dividing lines between samples and haplotypes
    for i in np.arange(1, gen_df.shape[0], 2):
        axarr[gen_i].axhline(i, lw=1, c="gainsboro", ls=":")
    for i in np.arange(0, gen_df.shape[0], 2):
        axarr[gen_i].axhline(i, lw=2, c="w")
    axarr[gen_i].axhline(4, lw=6, c="w")
    if gen_i != len(height_ratios) - 1:
        axarr[gen_i].set_xticks([])


# make custom legend
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

# get colors used in palette
custom_legend = []
for char, color in zip(uniq_chars[1:], cmap):
    custom_legend.append(mpatches.Patch(color=color))

key = pd.read_csv(f"tr_validation/trviz/ceph/{TRID}.{GENOME}.key.tsv", sep="\t", names=["motif", "encoded_motif", "count"])

encoded2motif = dict(zip(key["encoded_motif"], key["motif"]))

print (key)

axarr[-1].legend(
    custom_legend,
    key["motif"].to_list(),
    prop={"family": "monospace", "size": 10},
    frameon=True,
    shadow=True,
)

axarr[-1].set_xlabel("Motif number")
# f.suptitle(f"{TRID} ({GENOME})")
f.tight_layout()
f.savefig("trviz_haplotypes.png", dpi=200)
