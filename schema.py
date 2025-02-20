from pandera import Column, Check, DataFrameSchema

DeNovoSchema = DataFrameSchema({
    "sample_id": Column(str, required=True),
    "denovo_status": Column(str, Check.isin(["Y:-", "Y:=", "Y:+", "Y:?", "X"]), required=True),
    "trid": Column(str, required=True),
    "n_with_denovo_allele_strict": Column(int, coerce=True)
})

InheritedSchema = DataFrameSchema({
    "trid": Column(str, required=True),
    "sample_id": Column(str, required=True),
    "maternal_id": Column(str, required=True),
    "paternal_id": Column(str, required=True),
})