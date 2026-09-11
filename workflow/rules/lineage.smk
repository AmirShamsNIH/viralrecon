# ############################################################################
# lineage.smk — viral lineage and clade assignment
#
# Rules (all run in pinned Singularity images — no module load, no conda)
# -----------------------------------------------
#   pangolin_lineage   – SARS-CoV-2 lineage calling (Pango nomenclature)
#   nextclade_clade    – Nextclade clade assignment + QC
#
# Images live in config/containers.json and are pinned by version, because
# pangolin bundles its lineage database inside the image: floating tags would
# silently change lineage calls between runs.
# Requires snakemake --use-singularity (set by src/run.sh in the kickoff).
#
# Each caller runs where it can actually describe the target. Nextclade works
# for any virus whose reference carries a dataset, so it is gated on that.
# Pangolin implements Pango nomenclature, which exists for SARS-CoV-2 alone, so
# it is gated on the target's taxid. A target that qualifies for neither simply
# produces no lineage output rather than a call from the wrong virus.
# ############################################################################

from os.path import join
from scripts.common import allocated


# ── pangolin_lineage ──────────────────────────────────────────────────────────

rule pangolin_lineage:
    """
    Assign SARS-CoV-2 Pango lineage to each per-sample consensus FASTA.
    Runs in the pinned pangolin image; the module (pangolin/2.3.6) predates
    the current nomenclature entirely and cannot be used.

    @Input:  per-sample consensus FASTA from bcftools_consensus
    @Output: pangolin lineage CSV
    """
    input:
        fa = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                  "{sample}.{target}.consensus.fa"),
    output:
        csv = join(WORKPATH, "{sample}", "lineage", "{target}",
                   "{sample}.{target}.pangolin_lineage.csv"),
    params:
        rname  = "pangolin_lineage",
        outdir = join(WORKPATH, "{sample}", "lineage", "{target}"),
        extra  = config["parameters"]["lineage"]["pangolin"],
    log:
        join(WORKPATH, "logfiles", "lineage",
             "{sample}.{target}.pangolin_lineage.log"),
    resources:
        partition = allocated("partition", "pangolin_lineage", cluster),
        mem       = allocated("mem",       "pangolin_lineage", cluster),
        time      = allocated("time",      "pangolin_lineage", cluster),
    threads:
        int(allocated("threads", "pangolin_lineage", cluster))
    container:
        config["images"]["pangolin"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

# scorpio, which pangolin calls, uses multiprocessing.Manager(). That forks a
# helper which binds a Unix domain socket inside TMPDIR. The pipeline points
# TMPDIR at the output directory, which is GPFS, and AF_UNIX sockets do not
# work there: the helper cannot bind, the parent reads an empty pipe, and
# scorpio dies with EOFError. Confirmed by running the same consensus through
# the same image twice, changing only TMPDIR.
#
# So TMPDIR is redirected to a node-local path for this rule alone. This is
# not the pipeline's scratch policy being reversed: everything pangolin writes
# still goes to the output directory via --outdir, and what lands here is a
# socket plus a few KB, not intermediates anyone would want to inspect. It is
# /tmp rather than lscratch precisely so no lscratch allocation is needed.
PANGO_TMP=$(mktemp -d /tmp/viralrecon_pangolin.XXXXXX)
trap 'rm -rf "$PANGO_TMP"' EXIT
export TMPDIR="$PANGO_TMP"

# pangolin runs its own internal Snakemake and shells out to scorpio and
# usher. A non-zero exit from any of those is an annotation failure on a
# consensus that already exists - it says nothing about whether the assembly
# or the variant calls are sound. Letting it abort the run discards the whole
# report for every sample, which is the same failure mode already ruled out
# for a weak-mapping target. On failure, write a header-only CSV so the file
# the DAG promised exists, and record the reason for run_summary.tsv.
if ! pangolin {params.extra} \
        --threads {threads} \
        --outfile "{output.csv}" \
        "{input.fa}" >> "{log}" 2>&1; then
    echo "pangolin failed for {wildcards.sample} / {wildcards.target}; see {log}" | tee -a "{log}"
    printf 'taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note\n' > "{output.csv}"
    printf '%s,,,,,,,,,,,,,fail,pangolin_error,\n' "{wildcards.sample}.{wildcards.target}" >> "{output.csv}"
fi
"""


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
# here says nothing about the assembly. Stub both outputs rather than abort
# the run - see the same reasoning in pangolin_lineage.
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

