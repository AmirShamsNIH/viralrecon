# pre_process.smk: read repair, trimming, Kraken2 depletion and pre-alignment QC.

from os.path import join
from scripts.common import allocated


# ── bbtools_reformat ──────────────────────────────────────────────────────────

rule bbtools_reformat:
    """Repair mismatched R1/R2 pairing (paired only) and reformat FASTQ with BBtools,
    normalising quality encoding and enforcing length and quality limits."""
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

    # Step 1: repair mismatched R1/R2 read counts
    echo "[repair] Fixing pair order for {wildcards.sample}" >> "{log}" 2>&1
    bbtools repair \
        in1="{input.r1}" in2="{input.r2}" \
        out1="$REPAIRED_R1" out2="$REPAIRED_R2" \
        outs="$SINGLETONS" \
        repair=t overwrite=t \
        >> "{log}" 2>&1

    # Step 2: reformat / quality-normalise
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
    """Adapter trimming and quality filtering with fastp."""
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
# Split by tool because no single image carries kraken2, python and Krona.

rule kraken2_classify:
    """Classify trimmed reads with Kraken2. Keep the 150 GB request: with memory mapping,
    resident memory can approach the whole ~104 GB database regardless of input size."""
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
    """Deplete reads inside the kraken2_decon_taxids subtrees (human, phiX and bacteria by
    default), keeping unclassified and viral reads, and profile the library first."""
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
    """Interactive Krona chart of the pre-depletion Kraken2 profile, with an explicit
    -tax taxonomy because the Krona image ships only a placeholder."""
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
    """FastQC quality report on depleted reads."""
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


