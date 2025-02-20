import pandas as pd
import glob
import utils
import matplotlib.pyplot as plt

df = []
for fh in glob.glob("tr_validation/csv/filtered_and_merged/*.GRCh38.tsv"):
    d = pd.read_csv(fh, sep="\t").fillna({"children_with_denovo_allele": "unknown"})
    df.append(d)
df = pd.concat(df).fillna("unknown")

df = df[df["sample_id"] == 2216]

df = utils.filter_mutation_dataframe(
    df,
    remove_complex=False,
    remove_gp_ev=False,
    parental_overlap_frac_max=1,
    denovo_coverage_min=2,
    child_ratio_min=0.2,
)


df["chrom"] = df["trid"].apply(lambda t: t.split("_")[0])
df["start"] = df["trid"].apply(lambda t: t.split("_")[1])
df["end"] = df["trid"].apply(lambda t: t.split("_")[2])
df["region"] = df.apply(lambda row: row["chrom"] + ":" + row["start"] + "-" + row["end"], axis=1)

# phased = df[df["phase_summary"] != "unknown"].sample(100)
# phased[["region", "trid", "genotype", "motif", "phase_summary"]].to_csv("2216.phased.tsv", sep="\t", index=False)# out = pd.concat([untransmitted, transmitted])
prev = """chr1_101253182_101253964_trsolve
chr1_104976437_104976497_trsolve
chr1_97316513_97316562_trsolve
chr11_82499385_82499467_trsolve
chr12_119919699_119920004_trsolve
chr12_58148951_58149459_trsolve
chr15_51783741_51783793_trsolve
chr18_41675165_41675281_trsolve
chr18_53351286_53351537_trsolve
chr19_23761901_23761935_trsolve
chr2_162936686_162937236_trsolve
chr3_176605284_176605309_trsolve
chr4_140362375_140362387_trsolve
chr4_167199689_167199738_trsolve
chr5_106331013_106331040_trsolve
chr7_2386586_2386615_trsolve
chr7_42734099_42734214_trsolve
chr7_76996908_76999055_trsolve
chr8_130639150_130639181_trsolve
chr8_60429971_60430012_trsolve""".split()

df["prev"] = df["trid"].apply(lambda t: t in prev)
df = df[df["prev"] == False]

drop_cols = """denovo_al_in_parents	mom_dp	dad_dp	dad_ev	mom_ev	parental_overlap_coverage_total	parental_coverage_total	parental_overlap_coverage_frac	denovo_size	gp_ev	denovo_hap_id	str_chrom	denovo_al	non_denovo_al	denovo_al_diff	non_denovo_al_diff	kid_str_ps	kid_str_gt	kid_str_hap_a	kid_str_hap_b	str_midpoint	str_parent_of_origin	n_upstream  n_downstream    n_children	n_children_consistent	children_consistent""".split()
drop_cols.extend(["chrom", "start", "end"])
df = df.sort_values(["chrom", "start"]).drop(columns=drop_cols)
df.to_csv("out.tsv", sep="\t", index=False)
# print (df[["trid", "denovo_coverage", "child_coverage", "phase_summary"]])
# print (df.columns)
