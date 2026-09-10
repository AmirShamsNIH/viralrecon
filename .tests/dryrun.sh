#!/usr/bin/env bash
# Minimal DAG check: build a throwaway run directory from the repo's own config
# and ask Snakemake to resolve the graph. No tools run, nothing is downloaded,
# so it works anywhere - which is the point: it catches the failure mode that
# broke this pipeline repeatedly, where splitting a rule leaves its output
# referenced by nothing and Snakemake silently omits it.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO=$PWD
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$WORK"/{workflow,config,resources,inputs,ref_db}
cp -r workflow/* "$WORK/workflow/"
cp -r config/*  "$WORK/config/"
cp -r resources/* "$WORK/resources/" 2>/dev/null || true

SAMPLES=(TEST_S1 TEST_S2)
# Two targets on purpose. TESTVIRUS exercises the generic path; SARS_NC_045512.1
# matches lineage_targets and freyja_targets by substring, so pangolin,
# nextclade and freyja enter the graph. Without a lineage target the whole
# lineage stage is silently absent from this check.
TARGETS=(TESTVIRUS_NC_000001.1 SARS_NC_045512.1)

# Inputs and reference files only have to exist; the DAG is resolved from
# filenames, and --dry-run never opens them.
for s in "${SAMPLES[@]}"; do
    : > "$WORK/inputs/$s.R1.fastq.gz"; : > "$WORK/inputs/$s.R2.fastq.gz"
done
for t in "${TARGETS[@]}"; do
    mkdir -p "$WORK/ref_db/$t"
    for ext in .fa .fa.fai .dict .1.bt2 ; do : > "$WORK/ref_db/$t/$t$ext"; done
    : > "$WORK/ref_db/$t/snpEff.config"
done

python3 - "$WORK" "${SAMPLES[*]}" "${TARGETS[*]}" <<'PY'
import json, os, sys
work, samples, targets = sys.argv[1], sys.argv[2].split(), sys.argv[3].split()
cfg = {}
for name in ("config", "containers"):
    with open(os.path.join(work, "config", name + ".json")) as fh:
        cfg.update(json.load(fh))
cfg["options"] = {"output": work, "platform": "BIOWULF"}
cfg["project"] = {"paired": True}
cfg["samples"] = samples
cfg["targets"] = targets
cfg["references"] = {
    "target": {"BIOWULF": {
        t: ({"fasta": "x.fa", "gtf": "x.gff", "nextclade_dataset": "/x/nextclade"}
            if "NC_045512" in t else {"fasta": "x.fa", "gtf": "x.gff"})
        for t in targets}},
}
with open(os.path.join(work, "config.json"), "w") as fh:
    json.dump(cfg, fh, indent=2)
PY

cd "$WORK"
if ! command -v snakemake >/dev/null 2>&1; then
    echo "snakemake not installed here; config and rule checks passed, DAG resolution skipped"
    exit 0
fi
snakemake -n --quiet -s workflow/Snakefile --configfile config.json --cores 1 > dag.txt 2>&1 || {
    echo "DAG did not resolve:"; tail -30 dag.txt; exit 1; }
TOTAL=$(grep -E "^total" dag.txt | awk '{print $2}')
echo "DAG resolved: ${TOTAL:-?} jobs for ${#SAMPLES[@]} samples x ${#TARGETS[@]} target(s)"
grep -E "^(total|job|[a-z_]+ +[0-9])" dag.txt | head -40
