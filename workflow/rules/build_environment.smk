# build_environment.smk: symlink a reference prebuilt by `viralrecon build` into
# ref_db/{target}/ so every rule finds it at one path.

import os
from os.path import join
from scripts.common import allocated


rule custom_virmapDB:
    """Symlink the prebuilt reference files into ref_db/{target}/. Missing inputs fail
    before the rule runs, pointing back to `viralrecon build`."""
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
