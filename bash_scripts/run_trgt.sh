#!/bin/bash
set -e

${snakemake_input[trgt_binary]} genotype \
    --threads ${snakemake[threads]} \
    --genome ${snakemake_input[reference]} \
    --repeats ${snakemake_input[repeats]} \
    --reads ${snakemake_input[bam]} \
    ${snakemake_params[karyotype_cmd]} \
    --output-prefix ${snakemake_params[output_prefix]} \
    -v
