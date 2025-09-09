#!/bin/bash
set -e

module load bcftools

bcftools index --tbi ${snakemake_input[vcf]}