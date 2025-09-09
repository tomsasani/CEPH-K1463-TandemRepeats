#!/bin/bash
set -e

/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/hiphase-v1.4.4-x86_64-unknown-linux-gnu/hiphase \
    --threads ${snakemake[threads]} \
    --bam ${snakemake_input[bam]} \
    --vcf ${snakemake_input[snv_vcf]} \
    --vcf ${snakemake_input[str_vcf]} \
    --reference ${snakemake_input[reference]} \
    --output-vcf ${snakemake_output[snv_vcf]} \
    --output-vcf ${snakemake_output[str_vcf]}
