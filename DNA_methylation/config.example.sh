#!/usr/bin/env bash
# Adapt these paths to your environment

# Inputs
fq_dir="/path/to/paired_end_fastq"                # *_1.fq.gz/*_2.fq.gz
genome_dir="/path/to/bismark_genome"              # with spikein
region_bed="/path/to/Genome_classification.bed"   # non-overlapping regions
bismark_dir="/path/to/bismark"                    # v0.22.1

# Outputs and temporary files
methy_dir="/path/to/output/wgbs"
tmp_dir="/path/to/tmp"

# Analysis parameters
lambda_chr="Control_lambdaDNA"                    # spikein
lambda_min_cov=1
min_cov=5
max_jobs=5
align_threads=4
extractor_threads=2
