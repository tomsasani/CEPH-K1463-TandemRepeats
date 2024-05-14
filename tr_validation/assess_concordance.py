import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable, List
import numba
import tqdm
import matplotlib.patches as patches
from utils import filter_mutation_dataframe
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance


plt.rc("font", size=10)

ASSEMBLY = "GRCh38"
TECH = "ont"

def bootstrap(arr: np.ndarray, n: int = 100):
    nrow, ncol = arr.shape

    out = np.zeros((n, ncol))

    for i in range(n):
        resampled_idxs = np.random.choice(nrow, size=nrow, replace=True)
        out[i, :] = np.mean(arr[resampled_idxs, :], axis=0)
    ci_lo = np.percentile(out, q=1, axis=0)
    ci_hi = np.percentile(out, q=99, axis=0)
    return ci_lo, ci_hi



def main():

    mutations: List[pd.DataFrame] = []
    for fh in glob.glob(f"tr_validation/data/all_loci/orthogonal_support/*.read_support.csv"):
        tech = fh.split(".")[-3]
        assembly = fh.split(".")[-4]
        df = pd.read_csv(fh)
        df["tech"] = tech
        df["assembly"] = assembly
        mutations.append(df)
        
    mutations = pd.concat(mutations)
    mutations = mutations[mutations["sample_id"] == 2189]

    n_techs = mutations["tech"].nunique()
    n_samps = mutations["sample_id"].nunique()
    n_trids = mutations["trid"].nunique()

    tech2idx = dict(zip(mutations["tech"].unique(), range(n_techs)))
    smp2idx = dict(zip(mutations["sample_id"].unique(), range(n_samps)))
    trid2idx = dict(zip(mutations["trid"].unique(), range(n_trids)))

    # create array to store "off by" information
    MAX_OFFSET = 10
    off = np.zeros((n_techs, n_samps, n_trids, (MAX_OFFSET * 2) + 1))
    seen = np.zeros((n_techs, n_samps, n_trids, (MAX_OFFSET * 2) + 1))


    # figure out how "off" the reads are for this technology
    res = []
    for i, row in tqdm.tqdm(mutations.iterrows()):
        kid_ev = row["kid_evidence"].split("|")
        exp_diff = int(row["exp_allele_diff_denovo"])
        total_ev = sum([int(ev.split(":")[1]) for ev in kid_ev])
        for ev in kid_ev:
            diff, support = list(map(int, ev.split(":")))
            off_by = diff - exp_diff + MAX_OFFSET
            #if diff - exp_diff == 0 and row["tech"] == "element": print (support / total_ev)
            if abs(diff - exp_diff) > MAX_OFFSET: continue
            off[tech2idx[row["tech"]], smp2idx[row["sample_id"]], trid2idx[row["trid"]], off_by] = support / total_ev
            seen[tech2idx[row["tech"]], smp2idx[row["sample_id"]], trid2idx[row["trid"]], off_by] = 1

    f, axarr = plt.subplots(1, n_techs, figsize=(n_techs * 5, 5), sharey=True, )
    
    xvals = np.arange((MAX_OFFSET * 2) + 1)

    for tech, tech_i in tech2idx.items():
        for smp, smp_i in smp2idx.items():
            off_arr = off[tech_i, smp_i, :, :]
            seen_arr = seen[tech_i, smp_i, :, :]
            seen_idxs = np.where(np.sum(seen_arr, axis=1))[0]
            off_arr_tech = off_arr[seen_idxs, :]
            print (tech, np.mean(off_arr_tech[:, MAX_OFFSET]), np.std(off_arr_tech[:, MAX_OFFSET]))
            ci_lo, ci_hi = bootstrap(off_arr_tech)

            axarr[tech_i].plot(xvals, np.mean(off_arr_tech, axis=0))
            axarr[tech_i].fill_between(xvals, ci_lo, ci_hi, alpha=0.5, color="gainsboro")

        axarr[tech_i].set_title(tech.capitalize() if tech in ("illumina", "element") else tech.upper())
        axarr[tech_i].set_xticks(xvals[::5])
        axarr[tech_i].set_xticklabels(np.arange(-MAX_OFFSET, MAX_OFFSET + 1)[::5])
        axarr[tech_i].set_xlabel("Difference between measured and\nexpected allele length (bp)")
        axarr[tech_i].set_ylabel("Mean fraction of reads supporting\nmeasured allele length")
        sns.despine(ax=axarr[tech_i])
    f.suptitle("Average read support at homozygous homopolymer sites (n = 1,000)")
    f.tight_layout()
    
    f.savefig("orthogonal_stutter.png", dpi=200)


if __name__ == "__main__":
    
    main()
