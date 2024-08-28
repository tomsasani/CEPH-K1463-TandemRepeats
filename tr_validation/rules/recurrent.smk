rule get_recurrent_trids:
    input:
        fhs = expand("tr_validation/csv/filtered_and_merged/{SAMPLE}.{{ASSEMBLY}}.tsv", SAMPLE=CHILDREN),
        ped = "tr_validation/data/file_mapping.csv",
    output:
        recurrent_simple = "tr_validation/data/recurrent_trids.{ASSEMBLY}.simple.tsv",
        recurrent_full = "tr_validation/data/recurrent_trids.{ASSEMBLY}.complete.tsv"
    params: 
        annotations = lambda wildcards: ASSEMBLY2CATLOG[wildcards.ASSEMBLY],
    run:
        import pandas as pd

        dtype_dict = {"sample_id": str, "paternal_id": str, "maternal_id": str}

        dfs = []
        for fh in input.fhs:
            df = pd.read_csv(fh, sep="\t", dtype=dtype_dict)
            dfs.append(df)
        dfs = pd.concat(dfs)
        ped = pd.read_csv(input.ped, dtype=dtype_dict)
        dfs = dfs.merge(ped, on=list(dtype_dict.keys()))
        dfs["pedigree_id"] = dfs.apply(lambda row: assign_pedigree_id(row), axis=1)
        dfs["generation"] = dfs["sample_id"].apply(
                lambda s: (
                    "2"
                    if s in ("2209", "2188")
                    else (
                        "4"
                        if (s.startswith("200") and s not in ("200080", "200100"))
                        else "1" if s in ("2281", "2280", "2213", "2214") else "3"
                    )
                )
            )
        # remove loci at which the likely de novo size is 0
        dfs = dfs[dfs["likely_denovo_size"] != 0]

        # figure out the loci at which >1 sample has a denovo -- require at least one
        # to be in G3
        recurrence = dfs.groupby("trid").agg(recurrence=("generation", lambda g: (len(set(g)) > 1 or any([v > 1 for k, v in Counter(g).most_common()])) and any([k == "3" for k in g]))).reset_index()
        multiple_trids = recurrence[recurrence["recurrence"] == True]["trid"].unique()

        dfs = dfs[dfs["trid"].isin(multiple_trids)]

        # figure out the samples with de novo coverage
        res = []
        for trid, trid_df in dfs.groupby("trid"):
            samples_with_denovo = trid_df.query('denovo_coverage >= 2')["sample_id"].unique()
            generations_with_denovo = trid_df.query('denovo_coverage >= 2')["generation"].unique()
            trid_df["samples_with_denovo"] = ",".join(samples_with_denovo)
            trid_df["generations_with_denovo"] = ",".join(generations_with_denovo)

            res.append(trid_df)

        res = pd.concat(res)

        res["unique_denovo_events"] = res["samples_with_denovo"].apply(lambda s: len(s.split(",")))

        
        res.to_csv(output.recurrent_full, index=False, sep="\t")


        res[["trid", "index", "child_AL", "samples_with_denovo", "generations_with_denovo", "unique_denovo_events"]].drop_duplicates("trid").to_csv(output.recurrent_simple, index=False, sep="\t")


# query orthogonal BAMs for evidence at each of the raw denovo loci
rule validate_recurrent_with_orthogonal_tech:
    input:
        py_script="tr_validation/validate_with_bams.py",
        trids = "tr_validation/data/recurrent_trids.{ASSEMBLY}.simple.tsv",
        trids_complete = "tr_validation/data/recurrent_trids.{ASSEMBLY}.complete.tsv"

    output: "tr_validation/data/denovos/orthogonal_support_recurrent/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv",
    params:
        kid_bam = lambda wildcards: get_bam_fh(wildcards.SAMPLE, TECH2PREF[wildcards.TECH], TECH2SUFF[wildcards.TECH], SMP2ALT, use_alt=True if wildcards.TECH in ("ont", "hifi") else False),
        mom_bam = None,
        dad_bam = None,
    shell:
        """
        python {input.py_script} --mutations {input.trids} \
                                 --kid_bam {params.kid_bam} \
                                 --sample_id {wildcards.SAMPLE} \
                                 --out {output} \
                                 --tech {wildcards.TECH} \
                                 -mom_bam {params.mom_bam} \
                                 -dad_bam {params.dad_bam} \
                        
        """
