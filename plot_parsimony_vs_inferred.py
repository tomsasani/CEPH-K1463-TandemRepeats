import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import scipy.stats as ss
import tqdm


plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def check_phase(p):
    if p == "unknown": return p
    else:
        support = float(p.split(":")[1])
        if support < 0.75:
            return "unknown"
        else:
            return p.split(":")[0]


ASSEMBLY = "CHM13v2"

mutations = pd.read_csv("CHM13v2.filtered.tsv", sep="\t", dtype={"paternal_id": str})

mutations["phase"] = mutations["phase_consensus"].apply(lambda p: check_phase(p))
phased = mutations[mutations["phase"] != "unknown"]
hap_phased = phased[~phased["haplotype_in_parent_consensus"].str.contains("unknown")]

matching = hap_phased[hap_phased["likely_denovo_size"] == hap_phased["likely_denovo_size_parsimony"]]
matching["Both methods agree?"] = "Yes" + f" (n = {matching.shape[0]})"
different = hap_phased[hap_phased["likely_denovo_size"] != hap_phased["likely_denovo_size_parsimony"]]
different["Both methods agree?"] = "No" f" (n = {different.shape[0]})"

both = pd.concat([different, matching])

f, ax = plt.subplots(figsize=(8, 6))

sns.scatterplot(
    data=both,
    x="likely_denovo_size",
    y="likely_denovo_size_parsimony",
    hue="Both methods agree?",
    ax=ax,
    palette="colorblind",
    alpha=0.5
)
ax.axline((0, 0), slope=1, ls=":", c="gainsboro")
ax.set_xlabel("DNM size (inferred from parental haplotype)")
ax.set_ylabel("DNM size (inferred by allele length parsimony)")
ax.set_yscale("symlog")
ax.set_xscale("symlog")
sns.despine(ax=ax)
f.tight_layout()
f.savefig("parsimony.png", dpi=200)
