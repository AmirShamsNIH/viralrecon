# ############################################################################
# build_environment.smk — wire pre-built viral reference databases
#
# All indexing (bowtie2-build, samtools faidx/dict, snpEff build) is done
# OUTSIDE the pipeline by 'viralrecon build'.  This rule simply symlinks the
# pre-built files from the reference directory into ref_db/{target}/ so the
# rest of the pipeline can find them under a consistent path.
#
# Pre-built directory layout (produced by viralrecon build):
#   {genome_dir}/
#     {target}.fa          reference FASTA
#     {target}.fa.fai      samtools faidx index
#     {target}.dict        sequence dictionary
#     {target}.1.bt2       Bowtie2 index (+ .2 .3 .4 .rev.1 .rev.2)
#     genes.gff            GFF3 annotation
#     sequences.fa         snpEff copy of FASTA
#     genes.gff (copy)     snpEff copy of annotation
#     snpEff.config        snpEff database config
#
# Run 'viralrecon build --virus X --accession Y --output /data/refs'
# to produce the above layout before launching 'viralrecon run'.
# ############################################################################

import os
from os.path import join
from scripts.common import allocated


rule custom_virmapDB:
    """
    Wire pre-built reference into ref_db/{target}/ via symlinks.

    Inputs are the pre-built index files produced by 'viralrecon build'.
    Snakemake validates they exist before this rule runs — if any are missing,
    the pipeline fails with a clear missing-input error pointing to build.

    @Input:   Pre-built FASTA, fai, dict, bt2, snpEff.config from genome dir
    @Output:  Symlinks in ref_db/{target}/ pointing to each pre-built file
    """
    input:
        fasta   = lambda wc: config["references"]["target"][PLATFORM][wc.target]["fasta"],
        fai     = lambda wc: config["references"]["target"][PLATFORM][wc.target]["fasta"] + ".fai",
        seqdict = lambda wc: config["references"]["target"][PLATFORM][wc.target]["fasta"].replace(".fa", ".dict"),
        bt2     = lambda wc: config["references"]["target"][PLATFORM][wc.target]["fasta"].replace(".fa", ".1.bt2"),
        snpeff  = lambda wc: os.path.join(
            os.path.dirname(config["references"]["target"][PLATFORM][wc.target]["fasta"]),
            "snpEff.config",
        ),
    output:
        fasta   = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
        fai     = join(WORKPATH, "ref_db", "{target}", "{target}.fa.fai"),
        seqdict = join(WORKPATH, "ref_db", "{target}", "{target}.dict"),
        bt2     = join(WORKPATH, "ref_db", "{target}", "{target}.1.bt2"),
        snpeff  = join(WORKPATH, "ref_db", "{target}", "snpEff.config"),
    params:
        rname  = "virmapDB",
        srcdir = lambda wc: os.path.dirname(
            config["references"]["target"][PLATFORM][wc.target]["fasta"]
        ),
        dstdir = join(WORKPATH, "ref_db", "{target}"),
    log:
        join(WORKPATH, "logfiles", "build_environment", "{target}.log"),
    resources:
        partition = allocated("partition", "custom_virmapDB", cluster),
        mem       = allocated("mem",       "custom_virmapDB", cluster),
        time      = allocated("time",      "custom_virmapDB", cluster),
    threads: 1
    shell: """
set -euo pipefail
TARGET="{wildcards.target}"
SRC="{params.srcdir}"
DST="{params.dstdir}"
mkdir -p "$DST"

# Core reference files
ln -sf "$SRC/$TARGET.fa"      "$DST/$TARGET.fa"
ln -sf "$SRC/$TARGET.fa.fai"  "$DST/$TARGET.fa.fai"
ln -sf "$SRC/$TARGET.dict"    "$DST/$TARGET.dict"
ln -sf "$SRC/snpEff.config"   "$DST/snpEff.config"

# Bowtie2 index  (.bt2 = standard,  .bt2l = large index)
for bt2f in "$SRC/$TARGET"*.bt2 "$SRC/$TARGET"*.bt2l; do
    [ -f "$bt2f" ] && ln -sf "$bt2f" "$DST/$(basename "$bt2f")" || true
done

# snpEff data files (sequences.fa, genes.gff / genes.gtf)
for f in "$SRC/sequences.fa" "$SRC/genes.gff" "$SRC/genes.gtf"; do
    [ -f "$f" ] && ln -sf "$f" "$DST/$(basename "$f")" || true
done

# snpEff pre-built binary databases (required for annotation)
for f in "$SRC/snpEffectPredictor.bin" "$SRC/sequence.bin"; do
    [ -f "$f" ] && ln -sf "$f" "$DST/$(basename "$f")" || true
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Linked: $DST → $SRC" >> "{log}"
"""
