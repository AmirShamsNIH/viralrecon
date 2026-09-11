# ############################################################################
# pre_process.smk — read trimming, decontamination, and pre-alignment QC
#
# Rules
# -----
#   bbtools_reformat      – interleave / normalize paired FASTQ
#   fastp_trim            – adapter trimming and quality filtering
#   kraken2_classify      - Kraken2 classification against the standard DB
#   kraken2_decon         - subtractive host/phiX depletion by taxon
#   krona_plot            - interactive Krona chart of the profile
#   fastqc_pre            – per-sample FastQC (post-decon reads)
# ############################################################################

from os.path import join
from scripts.common import allocated


# ── bbtools_reformat ──────────────────────────────────────────────────────────

rule bbtools_reformat:
    """
    Repair + reformat paired FASTQ files with BBtools.

    Step 1 (paired only): repair.sh — restores proper R1/R2 pairing when
    the two files have mismatched read counts (common with clinical samples
    or files that were processed separately upstream).

    Step 2: reformat — normalises quality encoding, filters junk reads,
    enforces min-length / quality caps.

    @Input:  raw R1 / R2 fastq.gz symlinked into inputs/
    @Output: repaired + reformatted R1 / R2 fastq.gz
    """
    input:
        r1 = join(WORKPATH, "inputs", "{sample}.R1.fastq.gz"),
        r2 = join(WORKPATH, "inputs", "{sample}.R2.fastq.gz") if PAIRED else [],
    output:
        r1 = temp(join(WORKPATH, "{sample}", "pre_process",
                  "{sample}.bbtools_reformat.R1.fastq.gz")),
        r2 = temp(join(WORKPATH, "{sample}", "pre_process",
                  "{sample}.bbtools_reformat.R2.fastq.gz")) if PAIRED else [],
    params:
        rname   = "bbtools_reformat",
        paired  = PAIRED,
        extra   = config["parameters"]["pre_process"]["bbtools_reformat"],
        tmpdir  = join(WORKPATH, "{sample}", "pre_process", "repair_tmp"),
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.bbtools_reformat.log"),
    resources:
        partition = allocated("partition", "bbtools_reformat", cluster),
        mem       = allocated("mem",       "bbtools_reformat", cluster),
        time      = allocated("time",      "bbtools_reformat", cluster),
    threads:
        int(allocated("threads", "bbtools_reformat", cluster))
    container:
        config["images"]["bbtools"]
    shell: """
set -euo pipefail

if [ "{params.paired}" = "True" ]; then
    mkdir -p "{params.tmpdir}"
    REPAIRED_R1="{params.tmpdir}/{wildcards.sample}.repaired.R1.fastq.gz"
    REPAIRED_R2="{params.tmpdir}/{wildcards.sample}.repaired.R2.fastq.gz"
    SINGLETONS="{params.tmpdir}/{wildcards.sample}.singletons.fastq.gz"

    # Step 1 — repair mismatched R1/R2 read counts
    echo "[repair] Fixing pair order for {wildcards.sample}" >> "{log}" 2>&1
    bbtools repair \
        in1="{input.r1}" in2="{input.r2}" \
        out1="$REPAIRED_R1" out2="$REPAIRED_R2" \
        outs="$SINGLETONS" \
        repair=t overwrite=t \
        >> "{log}" 2>&1

    # Step 2 — reformat / quality-normalise
    echo "[reformat] Reformatting {wildcards.sample}" >> "{log}" 2>&1
    bbtools reformat {params.extra} \
        in1="$REPAIRED_R1" in2="$REPAIRED_R2" \
        out1="{output.r1}" out2="{output.r2}" \
        >> "{log}" 2>&1

    rm -rf "{params.tmpdir}"
else
    bbtools reformat {params.extra} \
        in="{input.r1}" out="{output.r1}" \
        >> "{log}" 2>&1
fi
"""


# ── fastp_trim ────────────────────────────────────────────────────────────────

rule fastp_trim:
    """
    Adapter trimming and quality filtering with fastp.

    @Input:  reformatted R1 / R2 from bbtools_reformat
    @Output: trimmed R1 / R2, fastp JSON + HTML reports
    """
    input:
        r1 = rules.bbtools_reformat.output.r1,
        r2 = rules.bbtools_reformat.output.r2 if PAIRED else [],
    output:
        r1   = temp(join(WORKPATH, "{sample}", "pre_process",
                    "{sample}.fastp_trim.R1.fastq.gz")),
        r2   = temp(join(WORKPATH, "{sample}", "pre_process",
                    "{sample}.fastp_trim.R2.fastq.gz")) if PAIRED else [],
        json = join(WORKPATH, "{sample}", "pre_process",
                    "{sample}.fastp_trim.json"),
        html = join(WORKPATH, "{sample}", "pre_process",
                    "{sample}.fastp_trim.html"),
    params:
        rname  = "fastp_trim",
        paired = PAIRED,
        extra  = config["parameters"]["pre_process"]["fastp_trim"],
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.fastp_trim.log"),
    resources:
        partition = allocated("partition", "fastp_trim", cluster),
        mem       = allocated("mem",       "fastp_trim", cluster),
        time      = allocated("time",      "fastp_trim", cluster),
    threads:
        int(allocated("threads", "fastp_trim", cluster))
    container:
        config["images"]["fastp"]
    shell: """
set -euo pipefail

if [ "{params.paired}" = "True" ]; then
    fastp {params.extra} --thread {threads} \\
        -i "{input.r1}" -I "{input.r2}" \\
        -o "{output.r1}" -O "{output.r2}" \\
        -j "{output.json}" -h "{output.html}" \\
        >> "{log}" 2>&1
else
    fastp {params.extra} --thread {threads} \\
        -i "{input.r1}" -o "{output.r1}" \\
        -j "{output.json}" -h "{output.html}" \\
        >> "{log}" 2>&1
fi
"""


# ── kraken2_classify / kraken2_decon / krona_plot ─────────────────────────────
#
# These three were a single rule while the pipeline ran on modules. A rule can
# declare only one container, and no single image carries kraken2, python and
# Krona, so the work is split along tool boundaries - which is also the
# Snakemake idiom and how nf-core and RNA-seek organise theirs.

rule kraken2_classify:
    """
    Classify trimmed reads against the standard Kraken2 database.

    On the 150 GB memory request, which looks wildly over-provisioned against
    any single measurement and is not. --memory-mapping means kraken2 mmaps
    hash.k2d rather than reading it into heap, so reported peak RSS is how much
    of that ~97 GiB file became resident, which depends on how many distinct
    hash buckets the reads touch and on the node's page-cache state - not on
    how many reads there are. Measured peaks across four runs were 79, 80, 90,
    92, 93, 93 and 94 GB, with no relation to input size: the single-end test
    on 65 MB of reads used 79 GB, while a 512 MB run used 94 GB.

    The worst case is the whole database resident, about 104 GB, so do not trim
    this request toward an observed peak. Anything near 100 GB will OOM the
    first time a run touches most of the database.

    --memory-mapping is also what lets several samples on one node share those
    pages through the page cache instead of each loading a private copy.

    @Input:  trimmed reads from fastp_trim
    @Output: Kraken2 report, plus per-read classifications (temp: several GB
             at this depth, consumed only by kraken2_decon)
    """
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2 if PAIRED else [],
    output:
        report = join(WORKPATH, "{sample}", "pre_process", "kraken2",
                      "{sample}.kraken2_decon.report.txt"),
        classified = temp(join(WORKPATH, "{sample}", "pre_process", "kraken2",
                               "{sample}.kraken2_classified.tsv")),
    params:
        rname  = "kraken2_classify",
        paired = PAIRED,
        extra  = config["parameters"]["pre_process"]["kraken2_decon"],
        db     = config["parameters"]["pre_process"]["kraken2_decon_db"],
        outdir = join(WORKPATH, "{sample}", "pre_process", "kraken2"),
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.kraken2_classify.log"),
    resources:
        partition = allocated("partition", "kraken2_classify", cluster),
        mem       = allocated("mem",       "kraken2_classify", cluster),
        time      = allocated("time",      "kraken2_classify", cluster),
    threads:
        int(allocated("threads", "kraken2_classify", cluster))
    container:
        config["images"]["kraken2"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

if [ "{params.paired}" = "True" ]; then
    kraken2 {params.extra} --threads {threads} \
        --db "{params.db}" --paired \
        --report "{output.report}" --output "{output.classified}" \
        "{input.r1}" "{input.r2}" >> "{log}" 2>&1
else
    kraken2 {params.extra} --threads {threads} \
        --db "{params.db}" \
        --report "{output.report}" --output "{output.classified}" \
        "{input.r1}" >> "{log}" 2>&1
fi
"""


rule kraken2_decon:
    """
    Deplete host / contaminant reads, and profile the library as received.

    Depletion is subtractive by taxon: reads are removed only when Kraken2
    assigned them inside the subtrees named in kraken2_decon_taxids - by
    default Homo sapiens (9606), phiX (10847) and Bacteria (2). Reads that are
    unclassified, or classified as viral, are always kept - the target virus is
    in the database, so a "keep unclassified" filter would discard exactly the
    reads this pipeline exists to analyse.

    Bacteria are depleted rather than merely reported. They stay visible in the
    composition table either way, since that profiles the library before any
    depletion, so nothing is lost from the QC picture by removing them from the
    reads that go on to alignment.

    @Input:  Kraken2 report + per-read classifications, trimmed reads
    @Output: depleted R1 / R2, depletion summary, pre-depletion composition
    """
    input:
        r1         = rules.fastp_trim.output.r1,
        r2         = rules.fastp_trim.output.r2 if PAIRED else [],
        report     = rules.kraken2_classify.output.report,
        classified = rules.kraken2_classify.output.classified,
        script = join(WORKPATH, "workflow", "scripts", "kraken_deplete.py"),
    output:
        r1      = join(WORKPATH, "{sample}", "pre_process",
                       "{sample}.kraken2_decon.R1.fastq.gz"),
        r2      = join(WORKPATH, "{sample}", "pre_process",
                       "{sample}.kraken2_decon.R2.fastq.gz") if PAIRED else [],
        summary = join(WORKPATH, "{sample}", "pre_process", "kraken2",
                       "{sample}.kraken2_decon.summary.tsv"),
        profile = join(WORKPATH, "{sample}", "pre_process", "kraken2",
                       "{sample}.kraken2_decon.composition.tsv"),
    params:
        rname      = "kraken2_decon",
        paired     = PAIRED,
        taxids     = config["parameters"]["pre_process"]["kraken2_decon_taxids"],
        # Target taxids come from genome.json; the config value is for extra
        # clades a user wants broken out on top of them.
        profile_ids= " ".join(TARGET_TAXIDS + [
                        t for t in config["parameters"]["pre_process"]
                        .get("kraken2_profile_taxids", "").split()
                        if t not in TARGET_TAXIDS]),
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.kraken2_decon.log"),
    resources:
        partition = allocated("partition", "kraken2_decon", cluster),
        mem       = allocated("mem",       "kraken2_decon", cluster),
        time      = allocated("time",      "kraken2_decon", cluster),
    threads:
        int(allocated("threads", "kraken2_decon", cluster))
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail

if [ "{params.paired}" = "True" ]; then
    python3 "{input.script}" \
        --kraken-output "{input.classified}" --report "{input.report}" \
        --deplete-taxids {params.taxids} \
        --in1 "{input.r1}" --in2 "{input.r2}" \
        --out1 "{output.r1}" --out2 "{output.r2}" \
        --summary "{output.summary}" \
        --profile "{output.profile}" --profile-taxids {params.profile_ids} \
        >> "{log}" 2>&1
else
    python3 "{input.script}" \
        --kraken-output "{input.classified}" --report "{input.report}" \
        --deplete-taxids {params.taxids} \
        --in1 "{input.r1}" --out1 "{output.r1}" \
        --summary "{output.summary}" \
        --profile "{output.profile}" --profile-taxids {params.profile_ids} \
        >> "{log}" 2>&1
fi
"""


rule krona_plot:
    """
    Interactive Krona chart of the pre-depletion Kraken2 profile.

    The Krona biocontainer ships a placeholder taxonomy, so the database is
    passed explicitly with -tax from our reference tree. Only taxonomy.tab
    (~123 MB) is needed for a report-based import; the full Biowulf taxonomy
    directory is 36 GB because of all.accession2taxid.sorted, which is not
    used here.

    @Input:  Kraken2 report from kraken2_classify
    @Output: Krona HTML
    """
    input:
        report = rules.kraken2_classify.output.report,
    output:
        krona = join(WORKPATH, "{sample}", "pre_process", "kraken2",
                     "{sample}.kraken2_decon.krona.html"),
    params:
        rname = "krona_plot",
        tax   = config["parameters"]["pre_process"]["krona_taxonomy"],
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.krona_plot.log"),
    resources:
        partition = allocated("partition", "krona_plot", cluster),
        mem       = allocated("mem",       "krona_plot", cluster),
        time      = allocated("time",      "krona_plot", cluster),
    threads:
        int(allocated("threads", "krona_plot", cluster))
    container:
        config["images"]["kronatools"]
    shell: """
set -euo pipefail
ktImportTaxonomy -tax "{params.tax}" -t 5 -m 3 \
    -o "{output.krona}" "{input.report}" >> "{log}" 2>&1
"""


# ── fastqc_pre ────────────────────────────────────────────────────────────────

rule fastqc_pre:
    """
    FastQC quality report on decontaminated reads.

    @Input:  depleted R1 (and R2) from kraken2_decon
    @Output: FastQC HTML + zip reports
    """
    input:
        r1 = rules.kraken2_decon.output.r1,
        r2 = rules.kraken2_decon.output.r2 if PAIRED else [],
    output:
        html_r1 = join(WORKPATH, "{sample}", "pre_process", "fastqc",
                       "{sample}.kraken2_decon.R1_fastqc.html"),
        zip_r1  = join(WORKPATH, "{sample}", "pre_process", "fastqc",
                       "{sample}.kraken2_decon.R1_fastqc.zip"),
    params:
        rname   = "fastqc_pre",
        outdir  = join(WORKPATH, "{sample}", "pre_process", "fastqc"),
        extra   = config["parameters"]["pre_process"]["fastqc"],
    log:
        join(WORKPATH, "logfiles", "pre_process", "{sample}.fastqc_pre.log"),
    resources:
        partition = allocated("partition", "fastqc_pre", cluster),
        mem       = allocated("mem",       "fastqc_pre", cluster),
        time      = allocated("time",      "fastqc_pre", cluster),
    threads:
        int(allocated("threads", "fastqc_pre", cluster))
    container:
        config["images"]["fastqc"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"
fastqc {params.extra} --threads {threads} --outdir "{params.outdir}" \\
    "{input.r1}" {input.r2} >> "{log}" 2>&1
"""


