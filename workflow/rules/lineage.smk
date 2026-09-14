# lineage.smk: Nextclade clade assignment for targets whose reference carries a dataset.
# Targets without one are skipped rather than described in another virus's terms.

from os.path import join
from scripts.common import allocated


# ── nextclade_clade ───────────────────────────────────────────────────────────

rule nextclade_clade:
    """Nextclade clade assignment, QC and mutation calling against the dataset fetched
    at build time, since compute nodes cannot reliably download it."""
    input:
        fa      = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.consensus.fa"),
    output:
        tsv = join(WORKPATH, "{sample}", "lineage", "{target}",
                   "{sample}.{target}.nextclade.tsv"),
        # MultiQC's nextclade module needs the semicolon CSV; the TSV is what the
        # summary collector and people read.
        csv = join(WORKPATH, "{sample}", "lineage", "{target}",
                   "{sample}.{target}.nextclade.csv"),
    params:
        rname      = "nextclade_clade",
        dataset_dir = lambda wc: nextclade_dataset_for(wc.target),
        outdir      = join(WORKPATH, "{sample}", "lineage", "{target}"),
        extra       = config["parameters"]["lineage"]["nextclade"],
    log:
        join(WORKPATH, "logfiles", "lineage",
             "{sample}.{target}.nextclade_clade.log"),
    resources:
        partition = allocated("partition", "nextclade_clade", cluster),
        mem       = allocated("mem",       "nextclade_clade", cluster),
        time      = allocated("time",      "nextclade_clade", cluster),
    threads:
        int(allocated("threads", "nextclade_clade", cluster))
    container:
        config["images"]["nextclade"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

# A failed clade call says nothing about the assembly: stub the outputs and let
# collect_final_report record the failure in qc_status.
if ! nextclade run {params.extra} \
        --input-dataset "{params.dataset_dir}" \
        --output-tsv "{output.tsv}" \
        --output-csv "{output.csv}" \
        "{input.fa}" >> "{log}" 2>&1; then
    echo "nextclade failed for {wildcards.sample} / {wildcards.target}; see {log}" | tee -a "{log}"
    printf 'seqName\tclade\tqc.overallStatus\n' > "{output.tsv}"
    printf '%s\t\tnextclade_error\n' "{wildcards.sample}.{wildcards.target}" >> "{output.tsv}"
    printf 'seqName;clade;qc.overallStatus\n' > "{output.csv}"
    printf '%s;;nextclade_error\n' "{wildcards.sample}.{wildcards.target}" >> "{output.csv}"
fi
"""

