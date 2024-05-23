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

ASSEMBLY = "CHM13v2"
SAMPLE = "2216"

def bootstrap(arr: np.ndarray, n: int = 1_000, pctile: float = 95):
    nrow, ncol = arr.shape

    out = np.zeros((n, ncol))

    for i in range(n):
        resampled_idxs = np.random.choice(nrow, size=nrow, replace=True)
        out[i, :] = np.mean(arr[resampled_idxs, :], axis=0)
    ci_lo = np.percentile(out, q=(100 - pctile) / 2., axis=0)
    ci_hi = np.percentile(out, q=pctile + ((100 - pctile) / 2.), axis=0)
    return ci_lo, ci_hi


def main():

    mutations: List[pd.DataFrame] = []
    for fh in glob.glob(f"tr_validation/data/all_loci/orthogonal_support/{SAMPLE}.{ASSEMBLY}.*.sample.read_support.csv"):
        tech = fh.split(".")[-4]

        df = pd.read_csv(fh)
        df = df[df["motif"].str.len() == 1]

        df["tech"] = tech
        mutations.append(df)

    mutations = pd.concat(mutations)

    n_techs = mutations["tech"].nunique()
    n_trids = mutations["trid"].nunique()

    tech2idx = dict(zip(mutations["tech"].unique(), range(n_techs)))
    trid2idx = dict(zip(mutations["trid"].unique(), range(n_trids)))

    # create array to store "off by" information
    MAX_OFFSET = 10
    off = np.zeros((n_techs, n_trids, (MAX_OFFSET * 2) + 1))
    seen = np.zeros((n_techs, n_trids, (MAX_OFFSET * 2) + 1))

    # figure out how "off" the reads are for this technology
    res = []
    for i, row in tqdm.tqdm(mutations.iterrows()):
        kid_ev = row["kid_evidence"].split("|")
        exp_diff = int(row["exp_allele_diff_denovo"])
        total_ev = sum([int(ev.split(":")[1]) for ev in kid_ev])
        for ev in kid_ev:
            diff, support = list(map(int, ev.split(":")))
            off_by = diff - exp_diff + MAX_OFFSET
            # if diff - exp_diff == 0 and row["tech"] == "element": print (support / total_ev)
            if abs(diff - exp_diff) > MAX_OFFSET: continue
            off[tech2idx[row["tech"]], trid2idx[row["trid"]], off_by] = support / total_ev
            seen[tech2idx[row["tech"]], trid2idx[row["trid"]], off_by] = 1

    f, ax = plt.subplots(figsize=(6, 5))

    xvals = np.arange((MAX_OFFSET * 2) + 1)

    palette = {"element": "goldenrod", "illumina": "cornflowerblue", "ont": "firebrick", "hifi": "darkgrey"}

    for tech, tech_i in tech2idx.items():
        off_arr = off[tech_i, :]
        seen_arr = seen[tech_i, :]
        seen_idxs = np.where(np.sum(seen_arr, axis=1))[0]
        print (tech, seen_idxs.shape)
        off_arr_tech = off_arr[seen_idxs, :]
        print (tech, np.mean(off_arr_tech[:, MAX_OFFSET]), np.std(off_arr_tech[:, MAX_OFFSET]))
        ci_lo, ci_hi = bootstrap(off_arr_tech)

        ax.plot(
            xvals,
            np.mean(off_arr_tech, axis=0),
            lw=2,
            label=(
                tech.capitalize()
                if tech in ("illumina", "element")
                else tech.upper() if tech == "ont" else "HiFi"
            ),
            color=palette[tech],
        )
        ax.fill_between(xvals, ci_lo, ci_hi, alpha=0.5, color="gainsboro")

        ax.set_xticks(xvals[::5])
        ax.set_xticklabels(np.arange(-MAX_OFFSET, MAX_OFFSET + 1)[::5])
        ax.set_xlabel("Difference between measured and\nexpected allele length (bp)")
        ax.set_ylabel("Mean fraction of reads supporting\nmeasured allele length")
        sns.despine(ax=ax)
    f.suptitle(f"Average read support at homozygous\nhomopolymer sites (n = 1,000) in {ASSEMBLY}")
    ax.legend()
    f.tight_layout()

    f.savefig("orthogonal_stutter.png", dpi=200)


if __name__ == "__main__":
    
    main()
