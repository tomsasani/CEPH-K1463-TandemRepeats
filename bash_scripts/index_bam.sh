#!/bin/bash
set -e

module load samtools

samtools index ${snakemake_input[bam]}
