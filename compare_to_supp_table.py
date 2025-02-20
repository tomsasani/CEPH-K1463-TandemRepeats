import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss


def calculate_poisson_ci(ct: int, obs: int, alpha: float = 0.95):

    l, u = (1 - alpha) / 2, 1 - ((1 - alpha) / 2)

    lower_bound = ss.chi2.ppf(l, ct * 2) / 2
    upper_bound = ss.chi2.ppf(u, (ct + 1) * 2) / 2

    # lower_bound = ss.poisson.ppf(l, ct)
    # upper_bound = ss.poisson.ppf(u, ct)
    
    # Convert bounds to rates per observation
    rate_lower = lower_bound / obs
    rate_upper = upper_bound / obs
    
    return (rate_lower, rate_upper)

pd.set_option("display.precision", 8)
plt.rc("font", size=16)

# define some global variables we'll access for filtering
FILTER_TRANSMITTED = False
PER_HAPLOTYPE = True
FILTER_RECURRENT = True
FILTER_ELEMENT = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

ORTHOGONAL_TECH = "element"

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"tr_validation/csv/orthogonal_support/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).fillna(
    {
        "children_with_denovo_allele": "unknown",
        "children_with_denovo_allele_strict": "unknown",
        "phase_summary": "unknown",
    }
)

# ensure we're looking at G3 DNMs
mutations = mutations[mutations["paternal_id"] == 2209]

# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"tr_validation/csv/rates/*.{ASSEMBLY}.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

# if we want to calculate mutation rates per-haplotype, rather
# than per-genome, we need to multiply the denominator by 2 to account
# for there being 2 copies of every locus in a diploid genome.
if PER_HAPLOTYPE:
    denoms["denominator"] *= 2

# if desired, remove recurrent sites that are recurrent in G3
if FILTER_RECURRENT:
    recurrent = pd.read_csv(f"tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv", sep="\t")
    recurrent["contains_G3"] = recurrent["samples_with_denovo"].apply(lambda samples: any([s in sample_ids for s in samples.split(",")]))
    recurrent = recurrent[recurrent["contains_G3"] == True]
    recurrent_trids = recurrent["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

print (mutations[(mutations["alt_sample_id"] == "NA12882") & (mutations["trid"].isin(["chr13_13805193_13805236_trsolve", "chr17_83874422_83874511_trsolve"]))])


# if desired, filter DNMs that didn't pass our Element validation checks.
# element validation data should already be in this dataframe as an extra annotation.
if FILTER_ELEMENT:
    # filter to STRs that had orthogonal data
    orthogonal = mutations[
        (mutations["simple_motif_size"] == "STR")
        & (mutations["validation_status"] != "no_data")
    ]
    # make sure we're only looking at sites with max AL <= 120
    orthogonal["max_al"] = orthogonal["child_AL"].apply(lambda a: max(map(int, a.split(","))))
    orthogonal = orthogonal[orthogonal["max_al"] <= 120]

    fail_trids = orthogonal[orthogonal["validation_status"] != "pass"]["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(fail_trids)]
print (mutations[(mutations["alt_sample_id"] == "NA12882") & (mutations["trid"].isin(["chr13_13805193_13805236_trsolve", "chr17_83874422_83874511_trsolve"]))])

# if desired, remove untransmitted DNMs in the samples for whom we can assess that
if FILTER_TRANSMITTED:

    has_transmission = mutations[mutations["sample_id"].isin(["2189", "2216"])]
    is_transmitted = has_transmission[has_transmission["children_with_denovo_allele"] != "unknown"]

    mutations = mutations[
        (~mutations["sample_id"].isin(["2189", "2216"]))
        | (mutations["children_with_denovo_allele"] != "unknown")
    ]

mutations["phase"] = mutations["phase_summary"].apply(lambda p: p.split(":")[0] if p != "unknown" else p)

# output filtered mutations to CSV
mutations.to_csv(f"tr_validation/csv/mutations.{ASSEMBLY}.filtered.tsv", index=False, sep="\t")


supp = pd.read_csv("supptable10.tsv", sep="\t", skiprows=1)
supp = supp[supp["recurrent_status"] == "not_recurrent"]
supp = supp[supp["element_validation_status"].isin(["pass", "not_evaluated"])]
supp = supp[~supp["manual_validation_status"].str.contains("FAIL")]
print (supp.shape, mutations.shape)


a_trids = set(mutations["trid"].unique())
b_trids = set(supp["trid"].unique())

print (b_trids.difference(a_trids))
print (a_trids.difference(b_trids))
a_diff = list(a_trids.difference(b_trids))
b_diff = list(b_trids.difference(a_trids))

print (mutations[mutations["trid"].isin(a_diff)][["child_coverage", "motif_size", "likely_denovo_size", "validation_status"]])
print (supp[supp["trid"].isin(b_diff)][["sample_id", "trid", "min_motif_size_in_locus", "likely_denovo_size"]])
# print (supp[supp["trid"].isin(a_diff)][["element_validation_status"]])