import pandas as pd
# from snakemake.script import snakemake
import utils


DTYPES = {
    "trid": "string",
    "genotype": "int",
    "index": "int",
    "denovo_coverage": "int",
    "child_coverage": "int",
    "child_ratio": "float",
    "child_AL": "string",
    "mother_AL": "string",
    "father_AL": "string",
    "per_allele_reads_mother": "string",
    "per_allele_reads_father": "string",
    "per_allele_reads_child": "string",
    "father_overlap_coverage": "string",
    "mother_overlap_coverage": "string",
}


# merge annotations with kid mutation dataframe
mutations = pd.read_csv(
    snakemake.input.kid_mutation_df,
    sep="\t",
    usecols=DTYPES.keys(),
    dtype=DTYPES,
)

# if we're filtering for candidate DNMs, use strict criteria. otherwise,
# just filter on depth.
_strict = snakemake.params.filtering_mode == "strict"

# remove duplicates
if not _strict:
    mutations = mutations.drop_duplicates(["trid"], keep="first")

mutations["sample_id"] = snakemake.wildcards.SAMPLE

file_mapping = pd.read_csv(snakemake.input.ped, dtype={"sample_id": "string"})
mutations = mutations.merge(file_mapping)

mutations = utils.filter_mutation_dataframe(
    mutations,
    remove_inherited=_strict,
    parental_overlap_frac_max=0.05 if _strict else 1,
    denovo_coverage_min=2 if _strict else 0,
    child_ratio_min=0.2 if _strict else 0,
    depth_min=10,
)

annotations = pd.read_csv(snakemake.input.annotations, sep="\t")

merged_mutations = mutations.merge(
    annotations,
    left_on="trid",
    right_on="TRid",
)

merged_mutations["motif_size"] = merged_mutations["min_motiflen"]
merged_mutations["simple_motif_size"] = merged_mutations.apply(
    lambda row: utils.determine_simplified_motif_size(row),
    axis=1,
)

merged_mutations.to_csv(snakemake.output.fh, sep="\t", index=False)
