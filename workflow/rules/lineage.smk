# ############################################################################
# lineage.smk — viral clade assignment
#
# Rules (all run in pinned Singularity images — no module load, no conda)
# -----------------------------------------------
#   nextclade_clade    – Nextclade clade assignment + QC
#
# The image is pinned by version in config/containers.json: nextclade's
# behaviour depends on the dataset it is given, and a floating tag could change
# how a dataset is interpreted between runs.
# Requires snakemake --use-singularity (set by src/run.sh in the kickoff).
#
# Nextclade is the only lineage caller. It works for any virus whose reference
# carries a dataset, which is why it is the one that survived: pangolin and
# freyja could each describe a fixed set of pathogens and nothing else, and a
# pipeline for any virus should not carry a stage most references cannot use.
# A target with no dataset is skipped rather than described against another
# virus's nomenclature.
# ############################################################################

from os.path import join
from scripts.common import allocated


# ── nextclade_clade ───────────────────────────────────────────────────────────

rule nextclade_clade:
    """
    Nextclade clade assignment, QC, and mutation calling.

    The dataset is read from a pre-fetched directory rather than downloaded at
    runtime: compute nodes reach the internet only through a per-session proxy,
    so `nextclade dataset get` inside a job is not reliably reachable. Refresh
    it deliberately with `nextclade dataset get` on a node that has the proxy.

    @Input:  per-sample consensus FASTA from bcftools_consensus
    @Output: nextclade TSV results
    """
    input:
        fa      = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.consensus.fa"),
    output:
        tsv = join(WORKPATH, "{sample}", "lineage", "{target}",
                   "{sample}.{target}.nextclade.tsv"),
        # MultiQC's nextclade module keys on the literal header "seqName;clade;",
        # i.e. the semicolon-delimited CSV. The TSV above is kept because it is
        # what the summary collector and humans read.
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

# A clade call is an annotation on a consensus that already exists; a failure
# here says nothing about the assembly. Stub both outputs rather than abort the
# run, and let collect_final_report turn the stub back into a qc_status reason
# so the failure is visible instead of silently blank.
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

