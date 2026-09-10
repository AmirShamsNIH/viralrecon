#!/bin/bash
set -euo pipefail

# ── Edit these ────────────────────────────────────────────────────────────────
VIRALRECON=/data/RTB_GRS/internal/Amir/viralrecon
GENOME_DIR=/data/shamsaddinisha/Test_Space/GRS_virmap/target_reference/data
# ─────────────────────────────────────────────────────────────────────────────

# No module load: viralrecon build runs every tool from a pinned
# Singularity image listed in config/containers.json.

$VIRALRECON/viralrecon build --virus MARBURG  --accession KM261523.1  --output $GENOME_DIR
$VIRALRECON/viralrecon build --virus SARS     --accession NC_045512.2 --output $GENOME_DIR
# $VIRALRECON/viralrecon build --virus BDBV  --accession FJ217161.1  --output $GENOME_DIR
# $VIRALRECON/viralrecon build --virus EBOLA --accession NC_002549.1 --output $GENOME_DIR

echo "Done. genome.json: $GENOME_DIR/genome.json"
