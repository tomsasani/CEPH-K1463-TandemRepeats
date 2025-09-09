#!/bin/bash
set -e

module load samtools

samtools sort -o ${snakemake_output[bam]} ${snakemake_input[bam]}
