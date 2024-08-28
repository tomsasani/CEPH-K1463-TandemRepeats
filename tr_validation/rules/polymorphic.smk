rule find_polymorphic_trids:
    input:
        repeat_gts = "/scratch/ucgd/lustre-work/quinlan/u0055382/src/CEPH-K1463-TandemRepeats/distance/T2T/50bp_merge/493ef25/xy_added/repeat_GTs.tsv"
    output:
        polymorphic_trids = "tr_validation/data/polymorphic.txt"
    run:
        import pandas as pd
        df = pd.read_csv(input.repeat_gts, sep="\t")
        df.columns = ["trid", "missing", "0", "0/0", "0/1", "1", "1/1", "1/2"]
        missing = df["missing"].values
        hom_ref = df["0"].values + df["0/0"].values
        hom_alt = df["1"].values + df["1/1"].values

        df["hom_ref"] = hom_ref
        df["hom_alt"] = hom_alt
        df["is_polymorphic"] = df.apply(lambda row: False if ((row["hom_ref"] == 28 - row["missing"]) or (row["hom_alt"] == 28 - row["missing"])) else True, axis=1)
        df[df["is_polymorphic"] == True][["trid"]].to_csv(output.polymorphic_trids, index=False)
