from cyvcf2 import VCF
import pandas as pd
import csv
import tqdm
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict, Counter
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 
import argparse
from utils import filter_mutation_dataframe

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


def main(args):

    KID_STR_VCF = VCF(args.str_vcf, gts012=True)

    CHROMS = [f"chr{c}" for c in range(1, 23)]

    informative_sites = pd.read_csv(
        args.informative_sites,
        dtype={"haplotype_A_parent": str, "haplotype_B_parent": str},
    )

    mutations = pd.read_csv(args.mutations, sep="\t")
    mutations = filter_mutation_dataframe(
        mutations,
        remove_complex=True,
        remove_duplicates=True,
    )

    phase_tree = make_tree(args.phase_blocks)

    print (f"Total of {mutations.shape[0]} DNMs")
    mutations["phase_block"] = mutations.apply(lambda row: in_phase_block(row, phase_tree), axis=1)
    mutations = mutations[mutations["phase_block"] != "UNK"]

    n_dnms_in_phase_blocks = mutations[mutations["phase_block"] != "UNK"]["trid"].nunique()
    print (f"Total of {n_dnms_in_phase_blocks} DNMs within phase blocks.")

    dnm_phases = []

    WRONG_TRID, NOT_PHASED, NON_AUTOSOMAL = 0, 0, 0

    # index let's us figure out which of the two alleles in the VCF have the allele

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
            is_phased = False
            try:
                hap_a, hap_b, is_phased = var.genotypes[0]
            except ValueError: 
                print ("### WARNING ###", var.genotypes[0])
                continue
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
                    "denovo_hap_id": denovo_hap_id,
                    "phase_block_chrom": phase_chrom,
                    "phase_block_start": phase_start,
                    "phase_block_end": phase_end,
                }
            )
            dnm_phases.append(row_dict)

    dnm_phases = pd.DataFrame(dnm_phases)

    n_phased_denovos = dnm_phases["trid"].nunique()
    print (f"Total of {n_phased_denovos} phased DNMs")
    print (f"{WRONG_TRID} with non-matching TRIDs, {NOT_PHASED} that weren't phased, {NON_AUTOSOMAL} not on autosomes.")

    merged_dnms_inf = dnm_phases.merge(
        informative_sites,
        on=["phase_block_chrom", "phase_block_start", "phase_block_end"],
    )

    merged_dnms_inf.dropna(inplace=True)

    merged_dnms_inf["phase"] = merged_dnms_inf.apply(lambda row: row["haplotype_{}_origin".format(row["denovo_hap_id"])], axis=1)

    # for a given informative site, we can now ask if we were able
    # to determine the actual *haplotype* the parent-of-origin donated
    # to the child
    merged_dnms_inf["phase_haplotype_origin"] = merged_dnms_inf.apply(lambda row: row[f"{row['phase']}_haplotype_origin"], axis=1,)

    # add midpoint of STR to dataframe
    merged_dnms_inf["str_midpoint"] = merged_dnms_inf["trid"].apply(lambda t: np.mean(list(map(int, t.split("_")[1:-1]))))
    merged_dnms_inf["diff_to_str"] = merged_dnms_inf["pos"] - merged_dnms_inf["str_midpoint"]

    # limit to DNMs with at least one upstream and one downstream site. require both to support
    # the same phase, otherwise ambiguous.

    merged_dnms_inf["is_upstream"] = merged_dnms_inf["diff_to_str"] < 0
    merged_dnms_inf["is_downstream"] = merged_dnms_inf["diff_to_str"] > 0

    TRID_COLS = ["trid", "index"]#, "phase_haplotype_origin"]

    updown = (
        merged_dnms_inf.groupby(TRID_COLS)
        .agg({"is_upstream": "sum", "is_downstream": "sum"})
        .reset_index()
        .rename(
            columns={
                "is_upstream": "n_upstream",
                "is_downstream": "n_downstream",
                "phase_haplotype_origin": "phase_haplotype_origin_summary",
            }
        )
    )

    merged_dnms_inf = merged_dnms_inf.merge(updown, on=TRID_COLS)

    # require a certain number of informative sites on either side
    MIN_INF = 1
    merged_dnms_inf = merged_dnms_inf.query(f"n_upstream >= {MIN_INF} and n_downstream >= {MIN_INF}")

    # for each DNM, create a new dataframe that only includes the informative sites that match the phase
    # at the focal (DNM) site.
    closest_inf_sites = []
    for trid, trid_df in merged_dnms_inf.groupby(TRID_COLS):
        # sort the informative sites by absolute distance to the STR
        trid_df["abs_diff_to_str"] = np.abs(trid_df["diff_to_str"])
        trid_df_sorted = trid_df.sort_values("abs_diff_to_str", ascending=True).reset_index()
        sorted_phases = trid_df_sorted["phase"].values

        # figure out how many sites are consistent closest to the STR. we can do this
        # simply by figuring out the first index where they are *inconsistent.*
        inconsistent_phases = np.where(sorted_phases[1:] != sorted_phases[:-1])[0]
        # if none are inconsistent, all are consistent!
        if inconsistent_phases.shape[0] == 0:
            final_consistent = -1
        else:
            final_consistent = inconsistent_phases[0]

        trid_df_sorted_consistent = trid_df_sorted.iloc[:final_consistent]

        consistent_diffs = trid_df_sorted_consistent["diff_to_str"]

        trid_df["n_consistent_upstream"] = np.sum(consistent_diffs < 0)
        trid_df["n_consistent_downstream"] = np.sum(consistent_diffs > 0)

        # now, using the informative sites with consistent phase, figureo out consistency of
        # the parental haplotype of origin.
        trid_df_sorted_consistent = trid_df_sorted_consistent[trid_df_sorted_consistent["phase_haplotype_origin"] != "?"]
        assert np.all(np.diff(trid_df_sorted_consistent["abs_diff_to_str"]) >= 0)

        sorted_haplotype_origins = trid_df_sorted_consistent["phase_haplotype_origin"].values

        inconsistent_haplotype_origins = np.where(sorted_haplotype_origins[1:] != sorted_haplotype_origins[:-1])[0]
        if inconsistent_haplotype_origins.shape[0] == 0:
            final_consistent = -1
        else:
            final_consistent = inconsistent_haplotype_origins[0]

        consistent_diffs = trid_df_sorted_consistent["diff_to_str"].values[:final_consistent]
        consistent_haplotype_origins = trid_df_sorted_consistent["phase_haplotype_origin"].values[:final_consistent]

        trid_df["n_consistent_upstream_origin"] = np.sum(consistent_diffs < 0)
        trid_df["n_consistent_downstream_origin"] = np.sum(consistent_diffs > 0)

        agg_hap_origins = Counter(consistent_haplotype_origins).most_common()
        agg_hap_origins = "|".join([":".join(list(map(str, el))) for el in agg_hap_origins])
        trid_df["phase_haplotype_origin_summary"] = agg_hap_origins

        closest_inf_sites.append(trid_df)
    merged_dnms_inf = pd.concat(closest_inf_sites)

    merged_dnms_inf["sample_id"] = args.focal
    merged_dnms_inf["true_phase"] = merged_dnms_inf["phase"]#.apply(lambda p: phase_dict[str(p)])

    merged_dnms_inf.drop_duplicates(TRID_COLS).drop(
        columns=[
            "phase_block_chrom",
            "phase_block_start",
            "phase_block_end",
            "haplotype_A_origin",
            "haplotype_B_origin",
            "dad_haplotype_origin",
            "mom_haplotype_origin",
            "phase_block",
            "chrom",
            "pos",
            "dad_gt",
            "mom_gt",
            "kid_gt",
            "phase_haplotype_origin",
            "str_midpoint",
            "diff_to_str",
            "is_upstream",
            "is_downstream",
            "abs_diff_to_str",
            "denovo_hap_id",
        ]
    ).fillna({"phase": "NA"}).to_csv(args.out, index=False)

    trids = merged_dnms_inf.query("n_upstream >= 10 and n_consistent_upstream <= 5 and n_downstream >= 10 and n_consistent_downstream <= 5")["trid"].unique()
    trids = merged_dnms_inf.query("n_upstream >= 20 and n_downstream >= 20")["trid"].unique()

    # print (len(merged_dnms_inf["trid"].unique()))
    # print (len(merged_dnms_inf[(merged_dnms_inf["n_consistent_upstream"] < merged_dnms_inf["n_upstream"]) | (merged_dnms_inf["n_consistent_downstream"] < merged_dnms_inf["n_upstream"])]["trid"].unique()))
    if trids.shape[0] > 0:
        focal_trid = trids[np.random.randint(len(trids))]
        merged_focal = merged_dnms_inf[merged_dnms_inf["trid"] == focal_trid]

        chrom, start, end, _ = focal_trid.split("_")
        focal_hap_id = merged_focal["denovo_hap_id"].unique()[0]
        xvals = merged_focal["pos"].values

        y1 = merged_focal[f"haplotype_A_origin"] == "dad"
        y2 = merged_focal[f"haplotype_B_origin"] == "dad"

        f, ax = plt.subplots(figsize=(9, 4))

        for y_ in np.unique(y1):
            y_idxs = np.where(y1 == y_)[0]

            phase = "dad" if y_ else "mom"

            xs = xvals[y_idxs]
            ys = [
                (
                    -0.25
                    if haplotype_origin == "A"
                    else 0 if haplotype_origin == "?" else 0.25
                )
                for haplotype_origin in merged_focal[
                    f"{phase}_haplotype_origin"
                ].values[y_idxs]
            ]

            ax.scatter(
                xs,
                ys,
                c="firebrick" if y_ else "dodgerblue",
                label="paternal" if y_ else "maternal",
                alpha=0.25,
            )

        for y_ in np.unique(y2):
            y_idxs = np.where(y2 == y_)[0]

            phase = "dad" if y_ else "mom"

            xs = xvals[y_idxs]
            ys = [
                (
                    0.75
                    if haplotype_origin == "A"
                    else 1 if haplotype_origin == "?" else 1.25
                )
                for haplotype_origin in merged_focal[
                    f"{phase}_haplotype_origin"
                ].values[y_idxs]
            ]

            ax.scatter(
                xs,
                ys,
                c="firebrick" if y_ else "dodgerblue",
                #label="paternal" if y_ else "maternal",
                alpha=0.25,
            )

        if focal_hap_id == "A":
            ax.scatter((int(end) + int(start)) / 2, 0, marker="s", c="k", label="DNM")
        else:
            ax.scatter((int(end) + int(start)) / 2, 1, marker="s", c="k", label="DNM")

        # for ax in (ax1, ax2):
        sns.despine(ax=ax, top=True, right=True)
        ax.set_yticks([-0.25, 0, 0.25, 0.75, 1, 1.25])
        ax.set_yticklabels(
            [
                "A (A)",
                "A (?)",
                "A (B)",
                "B (A)",
                "B (?)",
                "B (B)",
            ]
        )
        ax.set_ylabel("Haplotype ID")
        ax.legend(frameon=False)
        ax.set_xlabel("Position (bp)")
        ax.set_title(focal_trid)

        f.tight_layout()
        f.savefig(f'{args.focal}.phase_assignment.png', dpi=200)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--mutations")
    p.add_argument("--phase_blocks")
    p.add_argument("--informative_sites")
    p.add_argument("--ped")
    p.add_argument("--focal")
    p.add_argument("--out")
    p.add_argument("--str_vcf")
    p.add_argument("--dad_id")
    p.add_argument("--mom_id")
    args = p.parse_args()
    main(args)
