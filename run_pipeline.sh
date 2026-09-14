#!/bin/bash
set -euo pipefail

# ── Edit these ────────────────────────────────────────────────────────────────
VIRALRECON=/data/RTB_GRS/internal/pipeline/viralrecon
INPUT_DIR=/data/shamsaddinisha/Test_Space/GRS_virmap/virmap_fastq
OUTPUT_DIR=/data/shamsaddinisha/Test_Space/GRS_virmap/viralrecon_out
GENOME_JSON=/data/shamsaddinisha/Test_Space/GRS_virmap/target_reference/data/genome.json
TARGETS="MARBURG_KM261523.1 SARS_NC_045512.2"   # leave blank for all targets
# ─────────────────────────────────────────────────────────────────────────────

FASTQS=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" | sort)
if [ -z "$FASTQS" ]; then
    echo "ERROR: no *.fastq.gz files found in $INPUT_DIR" >&2
    exit 1
fi

$VIRALRECON/viralrecon run \
    --input  $FASTQS \
    --output $OUTPUT_DIR \
    --genome $GENOME_JSON \
    --targets $TARGETS

echo "Submitted. Monitor:"
echo "  tail -f $OUTPUT_DIR/logfiles/master.err"
echo "  squeue -u \$USER"
