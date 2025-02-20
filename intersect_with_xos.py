import pandas as pd
import glob

crossover_fh = "/scratch/ucgd/lustre-labs/quinlan/u1006375/old/ceph/recombination/2019-12-16-xos"

all_xos = []

samples = ["2216", "2211", "2212", "2298", "2215", "2217", "2189", "2187"]

for fh in glob.glob(f"{crossover_fh}/*.bed"):
    suff = fh.split("/")[-1]
    sample, chrom = suff.split(".")[:2]
    if sample not in samples: continue
    xos = pd.read_csv(fh, sep="\t")
    all_xos.append(xos)

all_xos = pd.concat(all_xos)
all_xos["chrom"] = all_xos["chrom"].apply(lambda c: f"chr{c}")

all_xos[["chrom", "start", "end"]].to_csv("xos.tsv", sep="\t", index=False)
print (all_xos)