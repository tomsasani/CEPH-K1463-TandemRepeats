#!/bin/bash
set -e

module load bcftools
bcftools view \
    ${snakemake_input[vcf]} \
    | grep "fileformat\|FILTER\|INFO\|FORMAT\|trgt\|bcftools\|${snakemake_wildcards[CHROM]},\|^${snakemake_wildcards[CHROM]}" \
    | bcftools sort -Ob -o ${snakemake_output[vcf]}
