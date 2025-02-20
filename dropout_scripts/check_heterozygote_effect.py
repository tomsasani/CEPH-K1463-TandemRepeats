import pandas as pd
import glob
import numpy as np
import tqdm
import csv
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns


def parent_is_heterozygote(row: pd.Series):

    phase = row["phase_summary"].split(":")[0]

    father_al = row["father_AL"].split(",")
    mother_al = row["mother_AL"].split(",")

    is_het = True
    if phase == "dad" and len(set(father_al)) == 1:
        is_het = False
    elif phase == "mom" and len(set(mother_al)) == 1:
        is_het = False

    return is_het

def get_parent_heterozygosity(row: pd.Series):

    father_al = row["father_AL"].split(",")
    mother_al = row["mother_AL"].split(",")

    return (len(set(father_al)) == 1) + (len(set(mother_al)) == 1)

ASSEMBLY = "GRCh38"

orig = []
for fh in glob.glob(f"tr_validation/downsampled/csv/2187.{ASSEMBLY}.*.phased.2gen.tsv"):
    trio_depths = fh.split("/")[-1].split(".")[2:5]
    if len(set(trio_depths)) > 1: continue
    print (fh)
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str, "paternal_id": str})
    df["downsampled_to"] = ".".join(trio_depths)
    orig.append(df)

orig = pd.concat(orig)
orig["poi_is_het"] = orig.apply(lambda row: parent_is_heterozygote(row), axis=1)
orig["num_parents_het"] = orig.apply(lambda row: get_parent_heterozygosity(row), axis=1)
MAX_MOTIF_SIZE = 6

res = []

for downsample, downsample_df in orig.groupby("downsampled_to"):
    denovo_trids = downsample_df["trid"].unique()
    dnm_num_parents_het = Counter(downsample_df["num_parents_het"]).most_common()
    # get complete list of sites in this individual
    fh = f"tr_validation/downsampled/csv/2187.{ASSEMBLY}.{downsample}.denovo.csv"
    num_parents_het = []
    with open(fh, "r") as infh:
        csvf = csv.reader(infh, delimiter="\t")
        header = None
        for i, l in tqdm.tqdm(enumerate(csvf)):
            # skip every other
            if i == 0:
                header = l
                continue
            if i % 2 == 0: continue
            d = dict(zip(header, l))
            trid = d["trid"]

            # remove if in the de novos
            if trid in denovo_trids: continue

            father_al = d["father_AL"].split(",")
            mother_al = d["mother_AL"].split(",")
            num_parents_het.append((len(set(father_al)) == 1) + (len(set(mother_al)) == 1))

    total = len(num_parents_het)
    for n, c in Counter(num_parents_het).most_common():
        res.append(
            {
                "denovo_status": "N",
                "downsample": downsample,
                "num_parents_het": n,
                "count": c,
                "frac": c / len(num_parents_het),
            }
        )
    for n, c in dnm_num_parents_het:
        res.append(
            {
                "denovo_status": "Y",
                "downsample": downsample,
                "num_parents_het": n,
                "count": c,
                "frac": c / len(downsample_df["num_parents_het"]),
            }
        )

res_df = pd.DataFrame(res)

g = sns.FacetGrid(data=res_df, row="downsample")
g.map(sns.barplot, "denovo_status", "frac", "num_parents_het")
g.add_legend()
g.savefig("bar.png", dpi=200)
