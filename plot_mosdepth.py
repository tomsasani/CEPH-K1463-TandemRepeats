import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns

res = []
for fh in glob.glob("tr_validation/data/mosdepth/*.summary.txt"):
    sample = fh.split("/")[-1].split(".")[0]
    generation = (
        "G4A"
        if sample.startswith("20008")
        else (
            "G4B"
            if sample.startswith("2001")
            else (
                "G2"
                if sample in ("2209", "2188")
                else "G1" if sample in ("2281", "2280", "2214", "2213") else "G3"
            )
        )
    )
    #if sample in ("200080", "200100"): generation = "G3_spouse"
    df = pd.read_csv(fh, sep="\t")
    total_length = df["length"].sum()
    total_bases = df["bases"].sum()
    mean_depth = total_bases / total_length
    res.append({"sample": sample, "depth": mean_depth, "generation": generation})

res_df = pd.DataFrame(res).sort_values("generation")
print (res_df)

f, ax = plt.subplots()
sns.barplot(data=res_df, x="sample", y="depth", hue="generation", ax=ax)
plt.xticks(rotation=90)
sns.despine(ax=ax)
ax.set_xlabel("Sample ID")
ax.set_ylabel("Mean genome-wide HiFi depth")
f.tight_layout()

pref = "/scratch/ucgd/lustre-labs/quinlan/u1006375/sasani-lab-notebook/img/dropout"


f.savefig("mosdepth.png", dpi=200)
f.savefig(f"{pref}/mosdepth.png", dpi=200)