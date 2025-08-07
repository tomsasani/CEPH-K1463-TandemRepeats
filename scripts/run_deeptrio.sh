#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deeptrio \
                --model_type PACBIO \
                --num_shards ${snakemake[threads]} \
                --reads_child ${snakemake_input[kid_bam]} \
                --reads_parent1 ${snakemake_input[dad_bam]} \
                --reads_parent2 ${snakemake_input[mom_bam]} \
                --output_vcf_child ${snakemake_output[kid_vcf]} \
                --output_vcf_parent1 ${snakemake_output[dad_vcf]} \
                --output_vcf_parent2 ${snakemake_output[mom_vcf]} \
                --sample_name_child ${snakemake_params[kid_name]} \
                --sample_name_parent1 ${snakemake_params[dad_name]} \
                --sample_name_parent2 ${snakemake_params[mom_name]} \
                --output_gvcf_child ${snakemake_output[kid_gvcf]} \
                --output_gvcf_parent1 ${snakemake_output[dad_gvcf]} \
                --output_gvcf_parent2 ${snakemake_output[mom_gvcf]} \
                --ref ${snakemake_input[ref]} \
                --regions ${snakemake_wildcards[CHROM]} \
                --make_examples_extra_args "select_variant_types='snps',min_mapping_quality=1"