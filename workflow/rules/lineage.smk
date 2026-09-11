# ############################################################################
# lineage.smk — viral lineage, clade, and variant-frequency analysis
#
# Rules (all run in pinned Singularity images — no module load, no conda)
# -----------------------------------------------
#   pangolin_lineage   – SARS-CoV-2 lineage calling (Pango nomenclature)
#   nextclade_clade    – Nextclade clade assignment + QC
#   freyja_demix       – Freyja variant decomposition (wastewater / mixed)
#
# Images live in config/containers.json and are pinned by version, because
# pangolin and freyja bundle their lineage / barcode databases inside the
# image: floating tags would silently change lineage calls between runs.
# Requires snakemake --use-singularity (set by src/run.sh in the kickoff).
#
# They run only when SARS-CoV-2 targets are present; skipped silently
# for other viruses because the Nextclade/Pango databases are SARS-specific.
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


# ── freyja_demix ─────────────────────────────────────────────────────────────

rule freyja_demix:
    """
    Freyja variant decomposition: estimates lineage abundances from a BAM.
    Primarily useful for wastewater surveillance or mixed clinical samples.

    `freyja variants` wraps samtools mpileup + iVar internally, so the whole
    step runs inside the freyja image with no external samtools/ivar needed.
    The SARS-CoV-2 barcode set is baked into the pinned image.

    @Input:  duplicate-marked BAM + reference FASTA
    @Output: freyja variants TSV, depths file, demix results
    """
    input:
        bam = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam"),
        bai = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam.bai"),
        fa  = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        variants = join(WORKPATH, "{sample}", "lineage", "{target}",
                        "{sample}.{target}.freyja.variants.tsv"),
        depths   = join(WORKPATH, "{sample}", "lineage", "{target}",
                        "{sample}.{target}.freyja.depths"),
        demix    = join(WORKPATH, "{sample}", "lineage", "{target}",
                        "{sample}.{target}.freyja.demix"),
        # MultiQC's freyja module searches "*.tsv" for the "summarized\t["
        # line, so the demix result is duplicated under a .tsv name. The
        # extension is the whole reason: content already matches.
        demix_tsv = join(WORKPATH, "{sample}", "lineage", "{target}",
                         "{sample}.{target}.freyja.demix.tsv"),
        # Bootstrap replicates give a confidence interval on each lineage
        # abundance. demix alone reports a point estimate, which reads as more
        # certain than it is when coverage is uneven or a lineage is rare.
        boot_lin  = join(WORKPATH, "{sample}", "lineage", "{target}",
                         "{sample}.{target}.freyja.boot_lineages.csv"),
        boot_sum  = join(WORKPATH, "{sample}", "lineage", "{target}",
                         "{sample}.{target}.freyja.boot_summarized.csv"),
    params:
        rname  = "freyja_demix",
        outdir = join(WORKPATH, "{sample}", "lineage", "{target}"),
        extra_variants = config["parameters"]["lineage"]["freyja_variants"],
        extra_demix    = config["parameters"]["lineage"]["freyja_demix"],
        extra_boot     = config["parameters"]["lineage"].get("freyja_boot", "--nb 100"),
        skip_boot      = config["parameters"]["lineage"].get("skip_freyja_boot", "false"),
        # Empty for a reference registered before --freyja-pathogen existed,
        # which leaves freyja on its own default of SARS-CoV-2 - the behaviour
        # those references were gated for.
        pathogen       = lambda wc: freyja_pathogen_for(wc.target),
        boot_base      = join(WORKPATH, "{sample}", "lineage", "{target}",
                              "{sample}.{target}.freyja.boot"),
        refname        = "{target}",
    log:
        join(WORKPATH, "logfiles", "lineage",
             "{sample}.{target}.freyja_demix.log"),
    resources:
        partition = allocated("partition", "freyja_demix", cluster),
        mem       = allocated("mem",       "freyja_demix", cluster),
        time      = allocated("time",      "freyja_demix", cluster),
    threads:
        int(allocated("threads", "freyja_demix", cluster))
    container:
        config["images"]["freyja"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

# Mixture deconvolution is an annotation on an existing alignment: a failure
# describes freyja, not the data, and must not discard the run. Each step is
# guarded, and a failure anywhere stubs every declared output and stops - see
# the same reasoning in pangolin_lineage.
FREYJA_OK=1

# Barcodes are per pathogen and keyed to that pathogen's own reference, so the
# pathogen travels with the reference rather than being a run-wide setting.
PATHOGEN_OPT=""
if [ -n "{params.pathogen}" ]; then
    PATHOGEN_OPT="--pathogen {params.pathogen}"
fi

# Step 1: variant calling + depth (samtools/iVar run inside freyja)
freyja variants {params.extra_variants} \
    --ref "{input.fa}" \
    --variants "{output.variants}" \
    --depths "{output.depths}" \
    "{input.bam}" >> "{log}" 2>&1 || FREYJA_OK=0

# freyja appends .tsv to --variants; normalise the name if it did
if [ ! -s "{output.variants}" ] && [ -s "{output.variants}.tsv" ]; then
    mv "{output.variants}.tsv" "{output.variants}"
fi

# Step 2: deconvolute lineage abundances against the bundled barcodes
if [ "$FREYJA_OK" = "1" ]; then
    freyja demix {params.extra_demix} $PATHOGEN_OPT \
        "{output.variants}" "{output.depths}" \
        --output "{output.demix}" >> "{log}" 2>&1 || FREYJA_OK=0
fi

if [ "$FREYJA_OK" = "0" ]; then
    echo "freyja failed for {wildcards.sample} / {wildcards.target}; see {log}" | tee -a "{log}"
    for f in "{output.variants}" "{output.depths}" "{output.demix}" \
             "{output.demix_tsv}" "{output.boot_lin}" "{output.boot_sum}"; do
        [ -e "$f" ] || : > "$f"
    done
    exit 0
fi

cp "{output.demix}" "{output.demix_tsv}"

# ── Step 3: bootstrap confidence intervals on the abundances ───────────────
if [ "{params.skip_boot}" = "true" ]; then
    echo "skip_freyja_boot=true; not bootstrapping" >> "{log}"
    : > "{output.boot_lin}"; : > "{output.boot_sum}"
else
    freyja boot {params.extra_boot} $PATHOGEN_OPT --nt {threads} \\
        --output_base "{params.boot_base}" \\
        "{output.variants}" "{output.depths}" >> "{log}" 2>&1
    # freyja writes <base>_lineages.csv / <base>_summarized.csv, which with this
    # boot_base is already exactly the declared output path. The previous `mv`
    # therefore had identical source and destination, failed with "are the same
    # file", and the `||` fallback truncated the file freyja had just written -
    # so every bootstrap was computed and then destroyed, and the empty result
    # read as freyja having produced nothing. Move only if the paths differ.
    if [ -s "{params.boot_base}_lineages.csv" ]; then
        [ "{params.boot_base}_lineages.csv" -ef "{output.boot_lin}" ] || \\
            mv "{params.boot_base}_lineages.csv" "{output.boot_lin}"
    else
        echo "freyja boot produced no lineages CSV" >> "{log}"; : > "{output.boot_lin}"
    fi
    if [ -s "{params.boot_base}_summarized.csv" ]; then
        [ "{params.boot_base}_summarized.csv" -ef "{output.boot_sum}" ] || \\
            mv "{params.boot_base}_summarized.csv" "{output.boot_sum}"
    else
        echo "freyja boot produced no summarized CSV" >> "{log}"; : > "{output.boot_sum}"
    fi
fi
"""


