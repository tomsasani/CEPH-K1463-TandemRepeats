#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deepvariant \
                --model_type PACBIO \
                --num_shards ${snakemake[threads]} \
                --reads ${snakemake_input[bam]} \
                --output_vcf ${snakemake_output[vcf]} \
                --sample_name ${snakemake_params[alt_sample_id]} \
                --output_gvcf ${snakemake_output[gvcf]} \
                --ref ${snakemake_input[ref]} \
                --regions ${snakemake_wildcards[CHROM]} \
                --make_examples_extra_args "select_variant_types='snps',min_mapping_quality=1"