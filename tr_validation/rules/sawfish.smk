rule prep_sawfish_calls:
    input:
        sawfish="tr_validation/data/{ASSEMBLY}_mergedDenovoSVs.tsv"
    output: 
        fh = "tr_validation/data/sawfish/{ASSEMBLY}_sawfish.{SAMPLE}.tsv"
    params:
        alt_sample_id=lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    run:
        import pandas as pd

        df = pd.read_csv(input.sawfish, sep="\t")
        df = df[df["source"] == "sawfish"]
        df["SVLEN"] = df["SVLEN"].apply(lambda svlen: max(svlen.split(";")))
        df[["start", "end", "SVLEN"]] = df[["start", "end", "SVLEN"]].astype(int)
        df = df[df["sample"] == params.alt_sample_id]

        df["adj_start"] = df["start"] - (df["SVLEN"] // 2)
        df["adj_end"] = df["end"] + (df["SVLEN"] // 2)

        df["trid"] = df.apply(lambda row: f'{row["seqnames"]}_{row["adj_start"]}_{row["adj_end"]}_sawfish', axis=1)
        df["child_AL"] = "0,0"
        df["index"] = 0
        
        df.to_csv(output.fh, sep="\t", index=False)


rule validate_sawfish_calls:
    input:
        kid_sawfish="tr_validation/data/sawfish/{ASSEMBLY}_sawfish.{SAMPLE}.tsv",
        py_script="tr_validation/validate_with_bams.py",
    output: "tr_validation/data/sawfish/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
    params:
        kid_bam=lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam=lambda wildcards: get_bam_fh(SMP2MOM[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2MOM[wildcards.SAMPLE] != "missing" else None,
        dad_bam=lambda wildcards: get_bam_fh(SMP2DAD[wildcards.SAMPLE], TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False) if SMP2DAD[wildcards.SAMPLE] != "missing" else None,
    shell:
        """
        python {input.py_script} --mutations {input.kid_sawfish} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam}
        """
