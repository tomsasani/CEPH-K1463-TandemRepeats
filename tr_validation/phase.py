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

    KID_DNMS = "tr_validation/data/df_transmission_C0_T0.parquet.gzip"
    DNM_PREF = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/493ef25/trgt-denovo/"
    KID_DNMS = f"{DNM_PREF}/2189_SF_GRCh38_50bp_merge_trgtdn.csv.gz"

    KID_PHASE_BLOCKS = "tr_validation/data/hiphase/NA12886.GRCh38.hiphase.blocks.tsv"

    CHROMS = [f"chr{c}" for c in range(1, 23)]
    #CHROMS = ["chr1"]

    informative_sites = pd.read_csv("inf_sites.csv")

    mutations = pd.read_csv(KID_DNMS, sep="\t")

    # if we're looking at de novos, ensure that we filter
    # the DataFrame to exclude the non-denovo candidates
    mutations = mutations[mutations["denovo_status"].isin(["Y:+", "Y:-"])]
    print (mutations.shape[0])
    mutations = mutations[(mutations["child_coverage"] >= 10)]# & (mutations["child_ratio"] >= 0.2) & (mutations["child_ratio"] <= 0.8)]
    mutations = mutations[(mutations["father_dropout"] == "N") & (mutations["mother_dropout"] == "N")]
    f, ax = plt.subplots()
    ax.hist(mutations["allele_ratio"])
    f.savefig("ratio.png")
    phase_tree = make_tree(KID_PHASE_BLOCKS)

    # filter to STRs in phase blocks. seems to remove around 100 / 2500 DNMs.
    print (f"Total of {mutations.shape[0]} DNMs")
    mutations["phase_block"] = mutations.apply(lambda row: in_phase_block(row, phase_tree), axis=1)
    mutations = mutations[mutations["phase_block"] != "UNK"]

    n_dnms_in_phase_blocks = mutations["trid"].nunique()
    print (f"Total of {n_dnms_in_phase_blocks} DNMs within phase blocks.")

    print (mutations.columns)

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
        #denovo_idx = row["index"]
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

            #print (hap_a, hap_b)
            #denovo_gt = hap_a if denovo_idx == 0 else hap_b 
            # figure out haplotype ID on which DNM resides
            if hap_a == denovo_gt:
                denovo_hap_id = "A"
            elif hap_b == denovo_gt:
                denovo_hap_id = "B"
            else:
                print (denovo_gt, hap_a, hap_b, is_phased)# += 1
                continue
            dnm_phases.append(
                {
                    "trid": trid,
                    "denovo_hap_id": denovo_hap_id,
                    "chrom": phase_chrom,
                    "start": phase_start,
                    "end": phase_end,
                }
            )

        

    dnm_phases = pd.DataFrame(dnm_phases)

    n_phased_denovos = dnm_phases["trid"].nunique()
    print (f"Total of {n_phased_denovos} phased DNMs")
    print (f"{WRONG_TRID} with non-matching TRIDs, {NOT_PHASED} that weren't phased, {NON_AUTOSOMAL} not on autosomes.")

    merged = dnm_phases.merge(informative_sites, on=["chrom", "start", "end"])

    # add midpoint of STR to dataframe
    merged["str_midpoint"] = merged["trid"].apply(lambda t: np.mean(list(map(int, t.split("_")[1:-1]))))

    N_INF = 5

    res = []
    for trid, trid_df in merged.groupby("trid"):
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
    merged_filtered = pd.concat(res)

    n_phased_denovos = merged["trid"].nunique()
    print (f"Total of {n_phased_denovos} confidently phased DNMs")

    trids = merged["trid"].unique()
    focal_trid = trids[np.random.randint(len(trids))]

    merged_focal = merged[merged["trid"] == focal_trid]

    merged_focal["phase"] = merged_focal.apply(lambda row: row["haplotype_{}_parent".format(row["denovo_hap_id"])], axis=1)
    #print (merged_focal)
    
    xvals = merged_focal["pos"].values
    y1 = merged_focal["haplotype_A_parent"] == 2209
    y2 = merged_focal["haplotype_B_parent"] == 2209

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

    merged_filtered["phase"] = merged_filtered.apply(lambda row: row["haplotype_{}_parent".format(row["denovo_hap_id"])], axis=1)

    # merge with original mutation file to see how well the `allele_origin` column corroborates
    # merged_orig = merged.merge(mutations, on="trid")

    # merged_orig["matches"] = merged_orig.apply(lambda row: row["allele_origin"].split(":")[0] == row["phase"], axis=1)
    # print (merged_orig.drop_duplicates("trid")[["allele_origin", "phase"]])
    # print (merged_orig.drop_duplicates("trid").groupby("matches").size())


    merged_filtered = merged_filtered.groupby(["trid", "phase"]).size().reset_index().rename(columns={0: "count"})
    merged_totals = merged_filtered.groupby("trid").agg({"count": "sum"}).reset_index().rename(columns={"count": "total"})
    merged_filtered = merged_filtered.merge(merged_totals, on="trid")
    merged_filtered["frac"] = merged_filtered["count"] / merged_filtered["total"]
    merged_filtered = merged_filtered.sort_values(["trid", "frac"], ascending=False).drop_duplicates(["trid"], keep="first")
    # merged = merged.drop_duplicates(["trid"])
    # at each informative site, examine the phase of the de novo allele with the parental alleles.
    print (merged_filtered.groupby("phase").size())

    f, ax = plt.subplots()
    sns.histplot(
        data=merged_filtered,
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

if __name__ == "__main__":
    main()
