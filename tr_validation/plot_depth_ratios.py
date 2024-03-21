import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable, List
import numba
from schema import DeNovoSchema
import cyvcf2
from collections import Counter
import glob
import matplotlib.patches as mpatches
import scipy.stats as ss
from statsmodels.stats.multitest import multipletests

from build_roc_curve_dnms import get_dp, annotate_with_allr, check_for_overlap

from utils import filter_mutation_dataframe
from assess_concordance import assign_allele

from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

plt.rc("font", size=12)

def get_chisquare(row):

    dad_support = list(map(int, row["per_allele_reads_father"].split(",")))
    mom_support = list(map(int, row["per_allele_reads_mother"].split(",")))

    if len(dad_support) < 2 or len(mom_support) < 2: return np.nan

    else:
        res = ss.chi2_contingency([dad_support, mom_support])
        return res.pvalue

def get_binomial_p(row):

    dad_support = list(map(int, row["per_allele_reads_father"].split(",")))
    mom_support = list(map(int, row["per_allele_reads_mother"].split(",")))

    res = ss.binomtest(sum(mom_support), sum(dad_support) + sum(mom_support), alternative="less")
    return res.pvalue
def main():
    # read in raw DNM files

    mutations = pd.read_csv(
        "tr_validation/csv/phased/combined/phased.2gen.tsv",
        sep="\t",
        dtype={"sample_id": str},
    )

    # if we're looking at de novos, ensure that we filter
    # the DataFrame to exclude the non-denovo candidates
    # remove denovos where allele size is unchanged
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=True,
        remove_duplicates=True,
    )

    MIN_ALLELE_COVERAGE = 1
    MIN_DENOVO_COVERAGE = 2
    MIN_ALLELE_RATIO = 0.
    MIN_CHILD_RATIO = 0.
    MAX_CHILD_RATIO = 1

    mutations = mutations.query(f"allele_coverage >= {MIN_ALLELE_COVERAGE} and \
                denovo_coverage >= {MIN_DENOVO_COVERAGE} and \
                allele_ratio >= {MIN_ALLELE_RATIO} and \
                child_ratio >= {MIN_CHILD_RATIO} and \
                child_ratio <= {MAX_CHILD_RATIO} and \
                father_dropout == 'N' and \
                mother_dropout == 'N'")

    mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4" if str(s).startswith("200") else "G2" if s in ("2209", "2188") else "G3")

    # add columns and filter
    mutations["mom_dp"] = mutations["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
    mutations["dad_dp"] = mutations["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )

    mutations["Minimum allele support in dad"] = mutations["per_allele_reads_father"].apply(lambda a: min(list(map(int, a.split(",")))))
    mutations["Minimum allele support in mom"] = mutations["per_allele_reads_mother"].apply(lambda a: min(list(map(int, a.split(",")))))
    mutations["Mom allele balance"] = mutations["per_allele_reads_mother"].apply(lambda a: list(map(int, a.split(",")))[0] / sum(list(map(int, a.split(",")))))


    mutations["chisquare"] = mutations.apply(lambda row: get_chisquare(row), axis=1)
    binomial_col = "Binomial p-value (mom DP vs. total DP, alternative = less)"
    mutations[binomial_col] = mutations.apply(lambda row: get_binomial_p(row), axis=1)
    # mutations[binomial_col] = multipletests(mutations[binomial_col], method="fdr_bh")[1]

    mutations.dropna(subset=["chisquare"], inplace=True)

    mutations["maternal_to_paternal_ratio"] = mutations["mom_dp"] / mutations["dad_dp"]

    val = binomial_col

    print (mutations[(mutations[val].between(0.75, 1.25)) & (mutations["sample_id"] == "200081")].head()[["trid", "genotype", "denovo_coverage", "phase_summary", "non_denovo_al_diff", "denovo_al_diff"]])

    mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0])

    f, axarr = plt.subplots(3, figsize=(10, 12))

    gen2idx = {"G2": 0, "G3": 1, "G4": 2}

    weird_ratios = ["200081", "200082", "200084", "200085", "200086", "200087"] + ["2209"]

    mutations["weird"] = mutations["sample_id"].apply(lambda s: "Y" if s in weird_ratios else "N")

    minval, maxval = np.min(mutations[val]), np.max(mutations[val])
    bins = np.linspace(minval, maxval, num=50)

    for (sample, gen), sample_df in mutations.groupby(["sample_id", "generation"]):
        vals = sample_df[val]

        hist, edges = np.histogram(vals, bins=bins)
        hist_fracs = hist / np.sum(hist)
        hist_cumsum = np.cumsum(hist_fracs)

        # axarr[gen2idx[gen]].plot(
        #     bins[:-1],
        #     hist_cumsum,
        #     color="firebrick" if sample in weird_ratios else "dodgerblue",
        # )
        axarr[gen2idx[gen]].hist(
            vals,
            bins=bins,
            color="firebrick" if sample in weird_ratios else "dodgerblue",
            alpha=0.2
        )

    legend_patches = [
        mpatches.Patch(color="firebrick", label="Y"),
        mpatches.Patch(color="dodgerblue", label="N"),
    ]

    for gen, idx in gen2idx.items():
        axarr[idx].set_title(f"Samples in {gen}")
        axarr[idx].set_ylabel("Number of DNMs")
        axarr[idx].set_xlabel(val)
        axarr[idx].legend(handles=legend_patches, title="Exhibits maternal DNM phasing bias?", frameon=False)
        axarr[idx].axvline(0.05, ls=":", c="k")
        #axarr[idx].set_xlim(0, 5)
        sns.despine(ax=axarr[idx])
    f.tight_layout()
    f.savefig('ratio.png', dpi=200)

    # mutations["denovo_overlap"] = mutations.apply(lambda row: check_for_overlap(row), axis=1)
    mutations["motif_size"] = mutations["motif"].apply(lambda m: len(m))

    print (mutations)


if __name__ == "__main__":
    # p = argparse.ArgumentParser()
    # p.add_argument(
    #     "--mutations",
    #     type=str,
    #     help="""Path to CSV file containing estimated allele lengths.""",
    # )
    # p.add_argument(
    #     "--element_validation",
    #     type=str,
    #     help="""Path to output image.""",
    # )
    # p.add_argument(
    #     "--ont_validation",
    #     type=str,
    #     default="illumina",
    #     help="""Technology used to generate allele length estimates.""",
    # )
    # p.add_argument("--sample_id")

    # args = p.parse_args()
    # main(args)
    main()
