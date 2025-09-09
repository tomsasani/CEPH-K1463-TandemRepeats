#!/bin/bash
set -e

module load bcftools

bcftools sort -Oz -o ${snakemake_output[vcf]} ${snakemake_input[vcf]}