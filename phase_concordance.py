import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import utils

ASSEMBLY = "GRCh38"
ASSEMBLY = "CHM13v2"

gen2 = []
for fh in glob.glob(f"tr_validation/csv/phased/*.{ASSEMBLY}.phased.2gen.tsv"):
    df = pd.read_csv(fh, sep="\t")
    gen2.append(df)
gen2 = pd.concat(gen2)

gen3 = []
for fh in glob.glob(f"tr_validation/csv/phased/*.{ASSEMBLY}.phased.3gen.tsv"):
    df = pd.read_csv(fh, sep="\t")
    gen3.append(df)
gen3 = pd.concat(gen3)

gen3 = gen3[gen3["most_common_freq"] >= 0.75]
gen3 = gen3[gen3["candidate_postzygotic"] == False]


print (gen3.shape)

print (gen3.groupby(["sample_id",  "phase"]).size())



combined = gen2.merge(gen3)
# combined = combined[combined["sample_id"].isin([2189, 2216, 2189])]

combined["2gen_phase"] = combined["phase_summary"].apply(lambda p: p.split(":")[0])
combined["3gen_phase"] = combined["phase_summary"].apply(lambda p: p.split(":")[0])
combined = combined[(combined["most_common_freq"] >= 0.75) & ((combined["n_upstream"] >= 1) & (combined["n_downstream"] >= 1))]

print (combined.groupby(["sample_id", "phase", "2gen_phase"]).size())
