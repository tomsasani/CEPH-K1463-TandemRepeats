import pandas as pd
import glob
import matplotlib.pyplot as plt

def get_complete_motif(row: pd.Series):
    complete_motif = ""
    motif_counts = list(map(int, row["denovo_MC"].split("_")))
    motifs = row["motifs"].split(",")
        
    assert len(motifs) == len(motif_counts)
    for m, mc in zip(motifs, motif_counts):
        complete_motif += "".join([m] * mc)
    return complete_motif


ASSEMBLY = "CHM13v2"
TECH = "hifi"

PASS = """chr1_54393726_54394070_trsolve
chr8_2376919_2377075_trsolve
chr7_2500010_2500042_trsolve
chr4_79949242_79949442_trsolve
chr4_21696993_21697153_trsolve
chr12_119907035_119907158_trsolve
chr12_114852499_114852706_trsolve
chr7_42892201_42892385_trsolve
chr21_33731357_33731465_trsolve
chr9_36529968_36530006_trsolve
chr7_6540708_6540973_trsolve
chr7_152489617_152489683_trsolve
chr12_95884953_95885246_trsolve
chr7_13334154_13334671_trsolve
chr15_32243116_32243499_trsolve
chr14_95031468_95031513_trsolve
chrX_49532133_49532368_trsolve
chrX_45925532_45926294_trsolve
chr8_91356817_91357059_trsolve
chr7_2490933_2492021_trsolve
chr7_19629376_19630006_trsolve
chr6_51818304_51818650_trsolve
chr2_121116093_121116213_trsolve
chr2_102413803_102413989_trsolve
chr20_38333021_38333121_trsolve
chr1_153340368_153340485_trsolve
chr19_6620632_6620703_trsolve
chr19_23902332_23902399_trsolve
chr13_71023191_71023231_trsolve
chr11_56923470_56923565_trsolve
chr10_409491_410951_trsolve
chr10_2390774_2391407_trsolve""".split()


mutations = pd.read_csv(f"tr_validation/data/recurrent_trids.{ASSEMBLY}.complete.tsv", sep="\t")
recurrent = mutations[mutations["trid"].isin(PASS)]
print (recurrent)

recurrent["unique_denovo_events"] = recurrent["samples_with_denovo"].apply(lambda s: len(s.split(",")))
recurrent["samples_with_denovo"] = recurrent["samples_with_denovo"].apply(lambda s: s.replace(",", "|"))
recurrent["generations_with_denovo"] = recurrent["generations_with_denovo"].apply(lambda s: s.replace(",", "|"))

recurrent = recurrent.query("unique_denovo_events >= 2")

# read in orthogonal evidence and filter to relevant TRIDs
# get mutations
# ortho = []
# for fh in glob.glob(f"tr_validation/data/denovos/orthogonal_support_recurrent/*.{ASSEMBLY}.{TECH}.read_support.csv"):
#     df = pd.read_csv(
#         fh,
#     )
#     sample_id = fh.split("/")[-1].split(".")[0]
#     df["sample_id"] = sample_id
#     ortho.append(df)

# ortho = pd.concat(ortho)

# need evidence in all members of pedigree
# count_per_trid = ortho.drop_duplicates(["sample_id", "trid"]).groupby("trid").size().to_dict()
# good_trids = {k for k,v in count_per_trid.items() if v == ortho["sample_id"].nunique()}

recurrent = recurrent[recurrent["trid"].isin(PASS)]

recurrent["chrom"] = recurrent["trid"].apply(lambda t: t.split("_")[0])
recurrent["start"] = recurrent["trid"].apply(lambda t: t.split("_")[1])
recurrent["end"] = recurrent["trid"].apply(lambda t: t.split("_")[2])

recurrent["Telomere"] = recurrent["Telomere"].apply(lambda t: "Y" if t == 1 else "N")
recurrent["Centromere"] = recurrent["Centromere"].apply(lambda t: "Y" if t == 1 else "N")

recurrent["denovo_MC"] = recurrent.apply(lambda row: row["child_MC"].split(",")[row["index"]], axis=1)

recurrent["complete_motif"] = recurrent.apply(lambda row: get_complete_motif(row), axis=1)

# output motifs
# motifs = recurrent["complete_motif"].to_list()

# outfh = open(f"recurrent.{ASSEMBLY}.motifs.fa", "w")
# for i, m in enumerate(motifs):
#     if len(m) < 8: continue
#     print (f">seq_{i}\n{m}", file=outfh)
# outfh.close()


recurrent_grouped = recurrent.groupby(
    [
        "chrom",
        "start",
        "end",
        "trid",
        "struc",
        "generations_with_denovo",
        "samples_with_denovo",
        "unique_denovo_events",
        "Telomere",
        "Centromere",
        # "nearest_gene",
        # "distance_to_tss",
    ]
).agg(
    denovo_als=("denovo_al", lambda a: "|".join(list(map(str, a)))),
    likely_denovo_sizes=("likely_denovo_size", lambda a: "|".join(list(map(str, a)))),
    denovo_al_minmax=("denovo_al", lambda a: f"{min(a)},{max(a)}"),
    denovo_size_minmax=("likely_denovo_size", lambda a: f"{min(a)},{max(a)}")

).reset_index()

column_dict = {
    "Telomere": "overlaps_telomere",
    "Centromere": "overlaps_centromere",
    "unique_denovo_events": "high_confidence_de_novo_events",
    "generations_with_denovo": "generations_with_de_novo_event",
    "struc": "motif_structure",
    "samples_with_denovo": "samples_with_de_novo_allele",
    "nearest_gene": "nearest_protein_coding_gene",
    "distance_to_tss": "distance_to_tss_nearest_gene",
}

recurrent_grouped.rename(columns=column_dict).sort_values(
    ["high_confidence_de_novo_events", "trid"],
    ascending=False,
).to_csv(f"recurrent.{ASSEMBLY}.summary.tsv", index=False, sep="\t")
