#!/usr/bin/env bash
#SBATCH --job-name=viralrecon
#SBATCH --partition=norm
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=4:00:00
#SBATCH --output=viralrecon.%j.log
set -euo pipefail

VIRALRECON=/data/RTB_GRS/internal/pipeline/viralrecon/viralrecon
WORKDIR=/data/CHANGE_ME/my_project
FASTQ_DIR=$WORKDIR/fastq
REF=$WORKDIR/target_reference
OUT=$WORKDIR/viralrecon_out

module load singularity python/3.10
export http_proxy=http://dtn20-e0:3128
export https_proxy=http://dtn20-e0:3128

if [ ! -d "$FASTQ_DIR" ]; then
    echo "FASTQ_DIR not found: $FASTQ_DIR" >&2
    exit 1
fi

$VIRALRECON build --platform BIOWULF --virus SARS --accession NC_045512.2 --output "$REF"

$VIRALRECON run --platform BIOWULF \
    --input "$FASTQ_DIR"/*.fastq.gz \
    --output "$OUT" \
    --genome "$REF/genome.json"
