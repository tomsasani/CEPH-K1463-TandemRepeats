from cyvcf2 import VCF
import cyvcf2
import pandas as pd
import csv
import tqdm
from typing import List
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import gzip
from schema import DeNovoSchema
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 
from build_roc_curve_dnms import get_dp


def make_tree(fh: str):
    """
    """

    tree = defaultdict(IntervalTree)
    is_zipped = fh.endswith(".gz")

    print("BUILDING TREE")

    with gzip.open(fh, "rt") if is_zipped else open(fh, "rt") as infh:
        csvf = csv.reader(infh, delimiter="\t")
        header = None
        for i,l in enumerate(csvf):
            if header is None:
                header = l
                continue
            d = dict(zip(header, l))
            chrom, start, end = d["chrom"], d["start"], d["end"]
            interval = Interval(int(start), int(end))
            tree[chrom].insert_interval(interval)

    return tree


def in_phase_block(row: pd.Series, phase_tree):
    try:
        chrom, start, end, _ = row["trid"].split("_")
    except ValueError:
        return "UNK"

    res = phase_tree[chrom].find(int(start) - 1, int(end))
    if len(res) == 0:
        return "UNK"
    else:
        return chrom + ":" + "-".join(list(map(str, [res[0].start, res[0].end])))


def main():
    KID_STR_VCF = VCF("tr_validation/data/hiphase/2189_SF_GRCh38_50bp_merge.sorted.phased.vcf.gz", gts012=True)

    DNM_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/"
    KID_DNMS = f"{DNM_PREF}/2189_SF_GRCh38_50bp_merge_trgtdn.csv.gz"

    KID_PHASE_BLOCKS = "tr_validation/data/hiphase/NA12886.GRCh38.hiphase.blocks.tsv"

    CHROMS = [f"chr{c}" for c in range(1, 23)]

    informative_sites = pd.read_csv("inf_sites.csv", dtype={"haplotype_A_parent": str, "haplotype_B_parent": str})

    mutations = pd.read_csv(KID_DNMS, sep="\t")

    mutations = mutations[mutations["denovo_coverage"] >= 1]
    # only look at autosomes
    mutations = mutations[~mutations["trid"].str.contains("X|Y|Un|random", regex=True)]
    mutations = mutations[~mutations["child_motif_counts"].str.contains("_")]

    #mutations = mutations[(mutations["child_dropout"] == "N") & (mutations["mother_dropout"] == "N") & (mutations["father_dropout"] == "N")]

    mutations["mom_coverage"] = mutations["per_allele_reads_mother"].apply(
        lambda a: get_dp(a)
    )
    mutations["dad_coverage"] = mutations["per_allele_reads_father"].apply(
        lambda a: get_dp(a)
    )

    phase_tree = make_tree(KID_PHASE_BLOCKS)

    print (f"Total of {mutations.shape[0]} DNMs")
    mutations["phase_block"] = mutations.apply(lambda row: in_phase_block(row, phase_tree), axis=1)
    mutations = mutations[mutations["phase_block"] != "UNK"]

    n_dnms_in_phase_blocks = mutations[mutations["phase_block"] != "UNK"]["trid"].nunique()
    print (f"Total of {n_dnms_in_phase_blocks} DNMs within phase blocks.")

    dnm_phases = []

    WRONG_TRID, NOT_PHASED, NON_AUTOSOMAL = 0, 0, 0

    # loop over mutations in the dataframe
    for _, row in tqdm.tqdm(mutations.iterrows()):
        # extract chrom, start and end
        trid = row["trid"]
        try:
            chrom, start, end, _ = trid.split("_")
            start, end = int(start) - 1, int(end)
        except ValueError:
            continue

        if chrom not in CHROMS:
            NON_AUTOSOMAL += 1
            continue

        region = f"{chrom}:{start}-{end}"

        # figure out genotype that is the de novo
        denovo_gt = row["genotype"]
        phase_block = row["phase_block"]
        phase_chrom = phase_block.split(":")[0]
        phase_start, phase_end = list(map(int, phase_block.split(":")[1].split("-")))

        # loop over VCF, allowing for slop to ensure we hit the STR
        for var in KID_STR_VCF(region):
            # ensure the VCF TRID matches the TR TRID
            var_trid = var.INFO.get("TRID")
            if var_trid != trid: 
                WRONG_TRID += 1
                continue

            # get phased genotype
            hap_a, hap_b, is_phased = var.genotypes[0]
            if not is_phased: 
                NOT_PHASED += 1
                continue

            # figure out haplotype ID on which DNM resides
            if hap_a == denovo_gt:
                denovo_hap_id = "A"
            elif hap_b == denovo_gt:
                denovo_hap_id = "B"
            else:
                continue

            row_dict = row.to_dict()
            row_dict.update(
                {
                    "trid": trid,
                    "denovo_hap_id": denovo_hap_id,
                    "chrom": phase_chrom,
                    "start": phase_start,
                    "end": phase_end,
                }
            )
            dnm_phases.append(row_dict)

    dnm_phases = pd.DataFrame(dnm_phases)

    n_phased_denovos = dnm_phases["trid"].nunique()
    print (f"Total of {n_phased_denovos} phased DNMs")
    print (f"{WRONG_TRID} with non-matching TRIDs, {NOT_PHASED} that weren't phased, {NON_AUTOSOMAL} not on autosomes.")

    merged_dnms_inf = dnm_phases.merge(informative_sites, on=["chrom", "start", "end"], how="left")

    merged_dnms_inf["phase"] = merged_dnms_inf.apply(lambda row: row["haplotype_{}_parent".format(row["denovo_hap_id"])], axis=1)
    print (merged_dnms_inf.groupby("phase").size())

    merged_dnms_inf.drop_duplicates(["trid"]).drop(
        columns=[
            "chrom",
            "start",
            "end",
            "haplotype_A_parent",
            "haplotype_B_parent",
            "phase_block",
            "pos",
            "alt_parent",
            "ref_parent",
            "gt",
            "haplotype_A_allele",
            "haplotype_B_allele",
            "denovo_hap_id",
        ]
    ).fillna({"phase": "NA"}).to_csv("phased.csv", index=False)

    print (merged_dnms_inf.dtypes)



    # add midpoint of STR to dataframe
    merged_dnms_inf["str_midpoint"] = merged_dnms_inf["trid"].apply(lambda t: np.mean(list(map(int, t.split("_")[1:-1]))))

    ###
    # require informative sites on either side?
    ###

    N_INF = 10

    res = []
    for trid, trid_df in merged_dnms_inf.groupby("trid"):
        # sort by informative site
        trid_df = trid_df.sort_values("pos")
        trid_df["diff_to_str"] = trid_df["pos"] - trid_df["str_midpoint"]
        # pick N sites up and downstream of the STR
        upstream_sites = trid_df.query("diff_to_str < 0")
        downstream_sites = trid_df.query("diff_to_str > 0")
        if upstream_sites.shape[0] < N_INF: continue
        if downstream_sites.shape[0] < N_INF: continue
        # get nearest N
        nearest_up = upstream_sites.tail(N_INF)
        nearest_dn = downstream_sites.head(N_INF)
        res.append(pd.concat([nearest_up, nearest_dn]))
        #break
    merged_dnms_inf = pd.concat(res)

    merged = merged_dnms_inf.groupby(["trid", "phase"]).size().reset_index().rename(columns={0: "count"})
    merged_totals = merged.groupby("trid").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})
    merged = merged.merge(merged_totals, on="trid")
    merged["frac"] = merged["count"] / merged["total"]
    merged = merged.sort_values(["trid", "frac"], ascending=False).drop_duplicates(["trid"], keep="first")

    # filter on number of INF sites that support
    merged = merged[merged["frac"] == 1.]
    # merged = merged.drop_duplicates(["trid"])
    # at each informative site, examine the phase of the de novo allele with the parental alleles.
    print (merged.groupby("phase").size())

    # at various filtering values, what is the ratio of paternal to maternal mutations?
    # coverage_thresh = np.arange(0, 30, 5)
    # denovo_thresh = np.arange(0, 10, 1)
    # res = []
    # for cov in coverage_thresh:
    #     for den in denovo_thresh:
    #         mutations_filt = merged_dnms_inf[(merged_dnms_inf["child_coverage"] >= cov) & (merged_dnms_inf["mom_coverage"] >= cov) & (merged_dnms_inf["dad_coverage"] >= cov)]
    #         mutations_filt = mutations_filt[mutations_filt["denovo_coverage"] >= den]

    #         mutations_filt["is_paternal"] = mutations_filt["phase"] == 2209
    #         #phased_counts = mutations_filt.groupby("phase").size().reset_index().rename(columns={0: "count"})
    #         #if phased_counts.shape[0] == 0: continue
            
    #         #phased_counts["fraction"] = phased_counts["count"] / mutations_filt.shape[0]
    #         #if "2209" not in phased_counts["phase"].to_list():
    #         #     paternal = 0.
    #         # else:
    #         #     paternal = phased_counts[phased_counts["phase"] == "2209"]["fraction"].values[0]
    #         phased_counts = mutations_filt[["is_paternal"]]
    #         phased_counts["cov"] = cov
    #         phased_counts["den"] = den
    #         res.append(phased_counts)
    #         # res.append({"paternal_frac": paternal, "Min. DP across trio": cov, "Min. DN coverage in child": den})
    # res_df = pd.concat(res)

    # print (res_df)

    # f, ax = plt.subplots()
    # sns.pointplot(x="den", y="is_paternal")
    # sns.scatterplot(data=res_df.query("paternal_frac > 0"), x="Min. DN coverage in child", y="paternal_frac", hue="Min. DP across trio", ax=ax)
    # ax.set_ylabel("Fraction of putative DNMs\nphased to paternal haplotype")
    # sns.despine(ax=ax, top=True, right=True)
    # f.tight_layout()
    # f.savefig("scatterphase.png", dpi=200)
    #mutations = mutations[mutations["allele_ratio"] >= 0.5]

    f, ax = plt.subplots()
    sns.histplot(
        data=merged,
        x="frac",
        hue="phase",
        stat="proportion",
        #cumulative=True,
        element="step",
        fill=False,
        ax=ax,
    )
    ax.set_xlabel("Fraction of informative sites\nthat support phase assignment")
    ax.set_ylabel("Fraction of DNMs")
    sns.despine(ax=ax, top=True, right=True)
    f.tight_layout()
    f.savefig('phase.png', dpi=200)

    n_phased_denovos = merged_dnms_inf["trid"].nunique()
    print (f"Total of {n_phased_denovos} confidently phased DNMs")

    trids = merged_dnms_inf["trid"].unique()
    focal_trid = trids[np.random.randint(len(trids))]

    merged_focal = merged_dnms_inf[merged_dnms_inf["trid"] == focal_trid]

    xvals = merged_focal["pos"].values
    y1 = merged_focal["haplotype_A_parent"] == "2209"
    y2 = merged_focal["haplotype_B_parent"] == "2209"

    chrom, start, end, _ = focal_trid.split("_")
    focal_hap_id = merged_focal["denovo_hap_id"].unique()[0]

    f, (ax1, ax2) = plt.subplots(2, figsize=(9, 5), sharex=True)

    ax1.scatter(xvals, y1, alpha=0.25)
    ax2.scatter(xvals, y2, alpha=0.25)

    if focal_hap_id == "A":
        ax1.scatter((int(end) + int(start)) / 2, 0.5, marker="s", c="k", label="DNM")
    else:
        ax2.scatter((int(end) + int(start)) / 2, 0.5, marker="s", c="k", label="DNM")

    for ax in (ax1, ax2):
        sns.despine(ax=ax, top=True, right=True)
        ax.set_yticks([-0.05, 1.05])
        ax.set_yticklabels(["maternal", "paternal"])
        ax.set_ylabel("Inferred phase")
        ax.legend(frameon=False)
    ax2.set_xlabel("Position (bp)")
    ax1.set_title("Haplotype A")
    ax2.set_title("Haplotype B")

    f.tight_layout()
    f.savefig('phase_assignment.png', dpi=200)


if __name__ == "__main__":
    main()
