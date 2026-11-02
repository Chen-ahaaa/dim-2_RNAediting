#!/usr/bin/env bash
# Copy to config.sh and adapt these paths to your environment.
# Do not commit config.sh if it contains private filesystem paths.

# Inputs
fq_dir="/path/to/paired_end_fastq"
genome_dir="/path/to/bismark_genome_with_lambda_spikein"
region_bed="/path/to/dup3kinds.bed"
bismark_dir="/path/to/bismark_v0.22.1"

# Outputs and temporary files
methy_dir="/path/to/output/wgbs"
tmp_dir="/path/to/tmp"

# Analysis parameters (retained from the original script)
lambda_chr="Control_lambdaDNA"
lambda_min_cov=1
min_cov=5
max_jobs=5
align_threads=4
extractor_threads=2
