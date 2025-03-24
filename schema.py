from pandera import Column, Check, DataFrameSchema


BaseSchema = DataFrameSchema(
    {
        "sample_id": Column(str, required=True),
        "trid": Column(str, required=True),
        "index": Column(int, Check.isin([0, 1]), required=True),
    }
)

TRGTDeNovoSchema = BaseSchema.add_columns({
    "genotype": Column(int, Check.isin([0, 1, 2]), required=True),
    "denovo_coverage": Column(int, required=True),
    "child_ratio": Column(float, required=True),
    "per_allele_reads_father": Column(str, required=True, coerce=True),
    "per_allele_reads_mother": Column(str, required=True),
    "child_coverage": Column(int, required=True),
    "father_overlap_coverage": Column(str, required=True, coerce=True),
    "mother_overlap_coverage": Column(str, required=True),
})


TransmissionSchema = BaseSchema.add_columns(
    {
        "n_children": Column(int, required=True),
        "n_children_consistent": Column(int, required=True),
        "children_consistent": Column(str, required=True, nullable=True),
        "children_also_denovo": Column(str, required=True, nullable=True),
        "children_with_denovo_allele": Column(str, required=True, nullable=True),
        "children_with_denovo_allele_strict": Column(str, required=True, nullable=True),
        "missing_gt_in_partner": Column(int, required=True),
    }
)


OrthogonalSchema = TRGTDeNovoSchema.add_columns({
    "kid_evidence": Column(str, required=True),
    "mom_evidence": Column(str, required=True),
    "dad_evidence": Column(str, required=True),
    "exp_allele_diff_denovo": Column(int, required=True),
    "exp_allele_diff_non_denovo": Column(int, required=True),
})


InformativeSiteSchema = TRGTDeNovoSchema.add_columns({
    "inf_chrom": Column(str, required=True),
    "inf_pos": Column(int, required=True, coerce=True),
    "dad_inf_gt": Column(int, required=True),
    "mom_inf_gt": Column(int, required=True),
    "str_parent_of_origin": Column(str, required=True),
    "abs_diff_to_str": Column(float, required=True, nullable=True)
})


HaplotypedSchema = InformativeSiteSchema.add_columns({
    "haplotype_in_parent": Column(str, required=True),
    "allele_length_in_parent": Column(float, required=True),
})


PhasedSchema = TRGTDeNovoSchema.add_columns({
    "phase_summary": Column(str, required=True),
    "allele_length_summary": Column(float, required=True),
})
