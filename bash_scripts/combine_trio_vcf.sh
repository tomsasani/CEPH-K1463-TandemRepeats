#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
    ${snakemake_input[sif]} \
    /usr/local/bin/glnexus_cli \
    --dir ${snakemake_params[gl_nexus_prefix]} \
    --config DeepVariant_unfiltered \
    --mem-gbytes 32 \
    --threads ${snakemake[threads]} \
    ${snakemake_input[kid_snp_gvcf]} \
    ${snakemake_input[mom_snp_gvcf]} \
    ${snakemake_input[dad_snp_gvcf]} > ${snakemake_output}