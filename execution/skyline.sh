#!/usr/bin/env bash
#SBATCH --job-name=viralrecon
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=4:00:00
#SBATCH --output=viralrecon.%j.log
set -euo pipefail

VIRALRECON=/data/openomics/viralrecon/viralrecon
WORKDIR=/data/CHANGE_ME/my_project
FASTQ_DIR=$WORKDIR/fastq
REF=$WORKDIR/target_reference
OUT=$WORKDIR/viralrecon_out

module load snakemake/7.22.0-ufanewz

if [ ! -d "$FASTQ_DIR" ]; then
    echo "FASTQ_DIR not found: $FASTQ_DIR" >&2
    exit 1
fi

$VIRALRECON build --platform SKYLINE --virus SARS --accession NC_045512.2 --output "$REF"

$VIRALRECON run --platform SKYLINE \
    --input "$FASTQ_DIR"/*.fastq.gz \
    --output "$OUT" \
    --genome "$REF/genome.json"
