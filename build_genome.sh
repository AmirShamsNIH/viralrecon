#!/bin/bash
set -euo pipefail

VIRALRECON=/data/RTB_GRS/internal/pipeline/viralrecon
GENOME_DIR=/data/shamsaddinisha/Test_Space/GRS_virmap/target_reference

$VIRALRECON/viralrecon build --virus MARBURG  --accession KM261523.1  --output $GENOME_DIR
$VIRALRECON/viralrecon build --virus SARS     --accession NC_045512.2 --output $GENOME_DIR

echo "Done. genome.json: $GENOME_DIR/genome.json"
