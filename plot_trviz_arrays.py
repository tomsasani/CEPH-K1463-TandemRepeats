from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from collections import defaultdict, Counter
import matplotlib.patches as patches


plt.rc("font", size=18, )

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

TRID = "chr8_2623352_2623487_trsolve"


GENOME = "GRCh38"

COHORT = "ceph"

PALETTE = "colorblind"

# store unique characters
smps = []

hapseq = []
for record in SeqIO.parse(f"trviz/ceph/{TRID}.{GENOME}_alignment_output.fa", "fasta"):
    hapseq.append(list(record.seq))
    smps.append(record.id)

hapseq = pd.DataFrame(hapseq)
print (hapseq)
# figure out how variable each motif across haplotypes
motif_variability = defaultdict(list)
for i, row in hapseq.iterrows():
    motif_counts = Counter(row.to_list()).most_common()
    for motif, count in motif_counts:
        if motif == "-": continue
        motif_variability[motif].append(count)


# remove "static" motifs if desired
static_motifs = []
for motif, cts in motif_variability.items():
    if len(set(cts)) == 1:# hapseq.shape[0]:
        static_motifs.append(motif)

# hapseq = hapseq.replace(to_replace={m: "-" for m in static_motifs})

orig_cols = hapseq.columns

# identify the unique characters in the encoded
# motif arrays, including "-" gaps
uniq_chars = np.unique(hapseq.values)
print (uniq_chars)

# map characters to colors
cmap = sns.color_palette(PALETTE, uniq_chars.shape[0] - 1)
cdict = dict(zip(uniq_chars[1:], cmap))


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

n_per_gen = hapseq.groupby("generation").size().to_dict()


gs_kw = dict(width_ratios=[1., 1], height_ratios=[1, 1., 3])
f, axarr = plt.subplot_mosaic(
    [
        ["G2A", "G4A"],
        ["G2B", "G4A"],
        ["G3", "G4B"],
    ],
    gridspec_kw=gs_kw,
    figsize=(18, 12),
    layout="constrained",
    sharex=True,
)

# f, axarr = plt.subplots(
#     hapseq["generation"].nunique() - 1,
#     sharex=True,
#     sharey=False,
#     figsize=(8, 18),
#     height_ratios=height_ratios,
# )

# loop over the generations
for gen_i, gen in enumerate(["G2A", "G2B", "G3", "G4A", "G4B"]):

    # collect the kids in this generation
    kid_df = hapseq[hapseq["generation"] == gen]


    # ignore spouses outside the family
    kid_df = kid_df[~kid_df["sample_id"].isin(["200080", "200100"])]
    
    # plot the kids second
    kid_df["order"] = 0
    # collect the parents in this generation
    par_df = hapseq[hapseq["parent_status"] == gen]
    # plot the parents first
    par_df["order"] = 1

    gen_df = pd.concat([par_df, kid_df]).sort_values("order")

    order_vals = gen_df["order"].values

    # get a list of haplotypes with DNMs
    has_denovo = gen_df["has_denovo"].values == "denovo"

    # convert characters to integers, including "-",
    # which will be converted to 0s
    gen_hapseq = gen_df[orig_cols].values

    # NOTE: specific for this TRID to ignore static motifs
    # gen_hapseq = gen_hapseq[:, :-6]

    n_haps, n_motifs = gen_hapseq.shape

    # get colors for the specific range of values seen in the hapseq array
    # ax = axarr[gen_i]
    ax = axarr[gen]

    ytick_vals = []

    cur_y = 0
    parent_spaced = False
    for hap_i, hap_seq in enumerate(gen_hapseq):
        # manual hacking to separate parents and individual samples
        # while keeping samples' haplotypes plotted together
        if order_vals[hap_i] == 1 and not parent_spaced:
            cur_y += 1.5
            parent_spaced = True
        if hap_i % 2 == 1:
            cur_y += 0.25
        if hap_i % 2 == 0:
            cur_y += 1
        
        ytick_vals.append(cur_y)

        for motif_i, motif in enumerate(hap_seq):
            if motif == "-": continue
            rect = patches.FancyBboxPatch(
                (n_motifs - motif_i, cur_y),
                0.33,
                0.33,
                linewidth=1,
                edgecolor="black",
                facecolor=cdict[motif],
                clip_on=False,
                boxstyle="round"
            )
            ax.add_patch(rect)
        cur_y += 1
    ax.set_xlim(0, n_motifs + 1.5)
    ax.set_ylim(0, max(ytick_vals) + 1.5)

    
    ax.set_yticks(ytick_vals)
    ax.set_yticklabels(gen_df["sample_id"].values, rotation=0)
    ax.set_title(gen)

    sns.despine(ax=ax)

    for i, t in enumerate(ax.yaxis.get_ticklabels()):
        if has_denovo[i] and t.get_text() not in par_df["sample_id"].to_list():
            t.set_color("red")
            t.set_fontstyle("oblique")
        if t.get_text() in par_df["sample_id"].to_list():
            t.set_fontweight("bold")


# make custom legend
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

print (uniq_chars)
# get colors used in palette
custom_legend = []
for char, color in zip(uniq_chars[1:], cmap[:-1]):
    custom_legend.append(mpatches.Patch(color=color))

key = pd.read_csv(f"trviz/ceph/{TRID}.{GENOME}.key.tsv", sep="\t", names=["motif", "encoded_motif", "count"])

encoded2motif = dict(zip(key["encoded_motif"], key["motif"]))

print (encoded2motif)

axarr["G4B"].legend(
    custom_legend,
    key["motif"].to_list(),
    prop={"family": "monospace", "size": 12},
    frameon=True,
    shadow=True,
    # loc="lower left",
    # bbox_to_anchor=(0.5, 1.05)
)

axarr["G4B"].set_xlabel("Motif number")
# f.suptitle(f"{TRID} ({GENOME})")
f.tight_layout()
f.savefig("trviz_haplotypes.png", dpi=200)
