import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import utils
from sklearn.mixture import GaussianMixture
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance


ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

TECH = "hifi"

# get mutations
mutations = []
for fh in glob.glob(f"tr_validation/data/denovos/orthogonal_support_recurrent/*.{ASSEMBLY}.{TECH}.read_support.csv"):
    print (fh)
    df = pd.read_csv(
        fh,
    )
    sample_id = fh.split("/")[-1].split(".")[0]
    df["sample_id"] = sample_id
    mutations.append(df)

mutations = pd.concat(mutations)

# need evidence in all members of pedigree
count_per_trid = mutations.drop_duplicates(["sample_id", "trid"]).groupby("trid").size().to_dict()
good_trids = {k for k,v in count_per_trid.items() if v == mutations["sample_id"].nunique()}

mutations = mutations[mutations["trid"].isin(good_trids)]
mutations = mutations.drop_duplicates(["sample_id", "trid", "index"])

print (mutations.groupby("trid").size())

for trid, trid_df in mutations.groupby("trid"):

    n_samples = trid_df["sample_id"].nunique()
    

    res = []

    # gather diffs in all members of the trio
    for sample_i, (sample, sample_df) in enumerate(trid_df.groupby("sample_id")):

        # if this sample has any de novo evidence, color it
        #is_denovo = np.any(sample_df["denovo_coverage"].values > 0)
        generation = 4 if sample.startswith("200") else 2 if sample in ("2188", "2209") else 1 if sample in ("2280", "2281", "2214", "2213") else 3

        # if is_denovo: 
        #     sample_df = sample_df[sample_df["denovo_coverage"] > 0].drop_duplicates("trid")
        # else:
        sample_df = sample_df.drop_duplicates("trid")

        assert sample_df.shape[0] == 1

        sample_diffs = []
        
        for diff_count in sample_df["kid_evidence"].values[0].split("|"):
            diff, count = list(map(int, diff_count.split(":")))
            for _ in range(count):
                res.append({"sample_id": sample, "diff": diff, "generation": generation})#, "has_denovo_coverage": is_denovo})

        

    res_df = pd.DataFrame(res)
    f, ax = plt.subplots()
    sns.stripplot(data=res_df.sort_values("generation"), y="sample_id", x="diff", ax=ax, alpha=0.5)
    ax.set_ylabel("Sample ID")
    ax.set_xlabel(f"Allele length (w/r/t {ASSEMBLY} genome)")
    sns.despine(ax=ax)
    ax.set_title(f"Distribution of inferred allele lengths in G2/3\n at {trid}\nusing HiFi reads")
    f.tight_layout()
    f.savefig(f'tr_validation/fig/recurrent/{ASSEMBLY}.{trid}.png', dpi=200)