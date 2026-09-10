# ############################################################################
# alignment.smk — map decontaminated reads to viral targets
#
# Rules
# -----
#   bowtie2_map      – align per-sample to each target, sort, mark-dup, add RG
#   flagstat_align   – samtools flagstat on final BAM
#   idxstats_align   – samtools idxstats on final BAM
# ############################################################################

from os.path import join
from scripts.common import allocated


# ── bowtie2_map ───────────────────────────────────────────────────────────────

rule bowtie2_map:
    """
    Align depleted reads to a viral target with Bowtie2, coordinate-sort, drop
    unmapped/secondary/supplementary records, optionally mark duplicates, and
    add a read group.

    Two things are deliberate here:

    * The flagstat is taken *before* filtering. It is the only place the true
      library mapping rate is visible; once unmapped records are dropped every
      downstream flagstat reads ~100% mapped. MultiQC is told to ignore it (see
      multiqc_ignore) so it does not appear as a second sample.
    * Duplicate marking is off by default. On a deep amplicon/enriched viral
      library nearly every read is a duplicate by position (98%+ observed here)
      and FreeBayes is run with --use-duplicate-reads, so the pass costs a full
      Picard traversal to produce a flag nothing acts on. Set
      skip_markduplicates=false to re-enable it.

    A sample whose mapped-read count falls below min_mapped_reads is flagged in
    the log and via a .qc_warn marker that reaches run_summary.tsv. It is not
    failed: one weak sample should not abort a whole run, and the consensus is
    depth-masked anyway.

    @Input:  depleted R1/R2 and target Bowtie2 index from ref_db/
    @Output: filtered (optionally dup-marked) sorted BAM + index, pre-filter
             flagstat, duplicate metrics
    """
    input:
        r1  = join(WORKPATH, "{sample}", "pre_process",
                   "{sample}.kraken2_decon.R1.fastq.gz"),
        r2  = join(WORKPATH, "{sample}", "pre_process",
                   "{sample}.kraken2_decon.R2.fastq.gz") if PAIRED else [],
        bt2 = join(WORKPATH, "ref_db", "{target}", "{target}.1.bt2"),
    output:
        # Handed to mark_duplicates, which produces the canonical
        # {sample}.{target}.bowtie2_map.bam that everything downstream reads.
        bam     = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                       "{sample}.{target}.aligned.bam")),
        bai     = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                       "{sample}.{target}.aligned.bam.bai")),
        raw_stat  = join(WORKPATH, "{sample}", "alignment", "{target}",
                         "{sample}.{target}.bowtie2_map.raw.flagstat"),
    params:
        rname       = "bowtie2_map",
        paired      = PAIRED,
        extra_bt2   = config["parameters"]["alignment"]["bowtie2_map"],
        extra_picard= config["parameters"]["alignment"]["picard_markdup"],
        filt        = config["parameters"]["alignment"].get("samtools_filter", "-F 0x004"),
        min_mapped  = config["parameters"]["alignment"].get("min_mapped_reads", "1000"),
        idx         = join(WORKPATH, "ref_db", "{target}", "{target}"),
        outdir      = join(WORKPATH, "{sample}", "alignment", "{target}"),
        tmpbam      = join(WORKPATH, "{sample}", "alignment", "{target}",
                           "{sample}.{target}.raw.bam"),
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.bowtie2_map.log"),
    resources:
        partition = allocated("partition", "bowtie2_map", cluster),
        mem       = allocated("mem",       "bowtie2_map", cluster),
        time      = allocated("time",      "bowtie2_map", cluster),
    threads:
        int(allocated("threads", "bowtie2_map", cluster))
    container:
        config["images"]["bowtie2"]
    shell: """
set -euo pipefail

mkdir -p "$(dirname "{output.bam}")"

# Align (with read-group) → sort → mark duplicates
if [ "{params.paired}" = "True" ]; then
    bowtie2 {params.extra_bt2} --threads {threads} \\
        --rg-id "{wildcards.sample}" --rg "LB:{wildcards.sample}" \\
        --rg "PL:ILLUMINA" --rg "PU:unit1" --rg "SM:{wildcards.sample}" \\
        -x "{params.idx}" -1 "{input.r1}" -2 "{input.r2}" \\
        2>> "{log}" \\
    | samtools sort -@ {threads} -o "{params.tmpbam}" - >> "{log}" 2>&1
else
    bowtie2 {params.extra_bt2} --threads {threads} \\
        --rg-id "{wildcards.sample}" --rg "LB:{wildcards.sample}" \\
        --rg "PL:ILLUMINA" --rg "PU:unit1" --rg "SM:{wildcards.sample}" \\
        -x "{params.idx}" -U "{input.r1}" \\
        2>> "{log}" \\
    | samtools sort -@ {threads} -o "{params.tmpbam}" - >> "{log}" 2>&1
fi

samtools index "{params.tmpbam}" >> "{log}" 2>&1

# ── Pre-filter flagstat: the only view of the real mapping rate ────────────
samtools flagstat --threads {threads} "{params.tmpbam}" > "{output.raw_stat}" 2>> "{log}"

# ── Minimum-mapped-reads gate (warn, never abort the whole run) ────────────
MAPPED=$(awk '/ mapped \(/ && !/primary/ {{print $1; exit}}' "{output.raw_stat}")
MAPPED=${{MAPPED:-0}}
if [ "$MAPPED" -lt "{params.min_mapped}" ]; then
    echo "QC WARNING: {wildcards.sample} / {wildcards.target}: $MAPPED mapped reads < min_mapped_reads={params.min_mapped}" \\
        | tee -a "{log}" > "{params.outdir}/{wildcards.sample}.{wildcards.target}.qc_warn"
fi

# ── Drop unmapped / secondary / supplementary before anything downstream ───
samtools view -b {params.filt} -@ {threads} \
    -o "{output.bam}" "{params.tmpbam}" >> "{log}" 2>&1
samtools index "{output.bam}" >> "{log}" 2>&1
rm -f "{params.tmpbam}" "{params.tmpbam}.bai"
"""


# ── mark_duplicates ───────────────────────────────────────────────────────────

rule mark_duplicates:
    """
    Optionally mark PCR duplicates, and produce the canonical alignment BAM.

    This rule always runs, so that everything downstream can depend on one
    filename regardless of configuration; what varies is whether Picard is
    invoked. With skip_markduplicates=true (the default, matching nf-core) the
    filtered BAM is hard-linked through and an empty metrics file is written.

    Marking is off by default because on a deep amplicon or enriched viral
    library nearly every read is a duplicate by position - 98.3% observed on
    this assay - and FreeBayes runs with --use-duplicate-reads, so the pass
    costs a full Picard traversal to produce a flag nothing acts on.

    Picard needs a JVM, which the bowtie2 image does not carry, so this is a
    separate rule rather than a branch inside bowtie2_map.

    @Input:  filtered BAM from bowtie2_map
    @Output: canonical {sample}.{target}.bowtie2_map.bam + index + metrics
    """
    input:
        bam = rules.bowtie2_map.output.bam,
        bai = rules.bowtie2_map.output.bai,
    output:
        bam       = join(WORKPATH, "{sample}", "alignment", "{target}",
                         "{sample}.{target}.bowtie2_map.bam"),
        bai       = join(WORKPATH, "{sample}", "alignment", "{target}",
                         "{sample}.{target}.bowtie2_map.bam.bai"),
        dupmetric = join(WORKPATH, "{sample}", "alignment", "{target}",
                         "{sample}.{target}.bowtie2_map.dup_metrics.txt"),
    params:
        rname        = "mark_duplicates",
        extra_picard = config["parameters"]["alignment"]["picard_markdup"],
        skip_markdup = config["parameters"]["alignment"].get("skip_markduplicates", "true"),
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.mark_duplicates.log"),
    resources:
        partition = allocated("partition", "mark_duplicates", cluster),
        mem       = allocated("mem",       "mark_duplicates", cluster),
        time      = allocated("time",      "mark_duplicates", cluster),
    threads:
        int(allocated("threads", "mark_duplicates", cluster))
    container:
        config["images"]["picard"]
    shell: """
set -euo pipefail

if [ "{params.skip_markdup}" = "true" ]; then
    echo "skip_markduplicates=true; not running Picard MarkDuplicates" >> "{log}"
    # Hard-link rather than copy: same filesystem, instant, and removing the
    # temp input later leaves this name intact.
    ln -f "{input.bam}" "{output.bam}" 2>/dev/null || cp "{input.bam}" "{output.bam}"
    ln -f "{input.bai}" "{output.bai}" 2>/dev/null || cp "{input.bai}" "{output.bai}"
    : > "{output.dupmetric}"
else
    picard MarkDuplicates {params.extra_picard} \
        I="{input.bam}" O="{output.bam}" M="{output.dupmetric}" \
        >> "{log}" 2>&1
    picard BuildBamIndex I="{output.bam}" >> "{log}" 2>&1 \
        || samtools index "{output.bam}" >> "{log}" 2>&1
fi
"""


# ── flagstat_align ────────────────────────────────────────────────────────────

rule flagstat_align:
    """
    samtools flagstat on the final aligned BAM.

    @Input:  duplicate-marked BAM from bowtie2_map
    @Output: flagstat text file
    """
    input:
        bam = rules.mark_duplicates.output.bam,
    output:
        stat = join(WORKPATH, "{sample}", "alignment", "{target}",
                    "{sample}.{target}.bowtie2_map.flagstat"),
    params:
        rname = "flagstat_align",
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.flagstat_align.log"),
    resources:
        partition = allocated("partition", "flagstat_align", cluster),
        mem       = allocated("mem",       "flagstat_align", cluster),
        time      = allocated("time",      "flagstat_align", cluster),
    threads: 2
    container:
        config["images"]["samtools"]
    shell: """
set -euo pipefail
samtools flagstat --threads {threads} "{input.bam}" > "{output.stat}" 2>> "{log}"
"""


# ── idxstats_align ────────────────────────────────────────────────────────────

rule idxstats_align:
    """
    samtools idxstats on the final aligned BAM.

    @Input:  duplicate-marked BAM + BAI from bowtie2_map
    @Output: idxstats TSV
    """
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bai,
    output:
        stat = join(WORKPATH, "{sample}", "alignment", "{target}",
                    "{sample}.{target}.bowtie2_map.idxstats"),
    params:
        rname = "idxstats_align",
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.idxstats_align.log"),
    resources:
        partition = allocated("partition", "idxstats_align", cluster),
        mem       = allocated("mem",       "idxstats_align", cluster),
        time      = allocated("time",      "idxstats_align", cluster),
    threads: 1
    container:
        config["images"]["samtools"]
    shell: """
set -euo pipefail
samtools idxstats "{input.bam}" > "{output.stat}" 2>> "{log}"
"""


# ── mosdepth_align ───────────────────────────────────────────────────────────

rule mosdepth_align:
    """
    Per-base and windowed coverage depth with mosdepth.

    @Input:  duplicate-marked BAM from bowtie2_map
    @Output: mosdepth summary, per-base bed.gz, global dist
    """
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bai,
    output:
        summary  = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.mosdepth.summary.txt"),
        global_d = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.mosdepth.global.dist.txt"),
        per_base = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.per-base.bed.gz"),
        # --by 200 windows: MultiQC's mosdepth module renders these as the
        # coverage-across-genome plot, which is what nf-core produces with a
        # bespoke R script.
        regions  = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.regions.bed.gz"),
    params:
        rname  = "mosdepth_align",
        prefix = join(WORKPATH, "{sample}", "alignment", "{target}",
                      "mosdepth", "{sample}.{target}"),
        extra  = config["parameters"]["alignment"]["mosdepth"],
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.mosdepth_align.log"),
    resources:
        partition = allocated("partition", "mosdepth_align", cluster),
        mem       = allocated("mem",       "mosdepth_align", cluster),
        time      = allocated("time",      "mosdepth_align", cluster),
    threads:
        int(allocated("threads", "mosdepth_align", cluster))
    container:
        config["images"]["mosdepth"]
    shell: """
set -euo pipefail
mkdir -p "$(dirname "{output.summary}")"
mosdepth {params.extra} --threads {threads} \
    "{params.prefix}" "{input.bam}" >> "{log}" 2>&1
"""


# ── picard_collect_metrics ────────────────────────────────────────────────────

rule picard_collect_metrics:
    """
    Picard CollectMultipleMetrics: insert size, alignment summary, GC bias.

    @Input:  duplicate-marked BAM + reference FASTA
    @Output: insert-size metrics/histogram, alignment summary, GC-bias metrics
    """
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bai,
        fa  = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        insert_metrics = join(WORKPATH, "{sample}", "alignment", "{target}",
                              "picard", "{sample}.{target}.insert_size_metrics"),
        insert_hist    = join(WORKPATH, "{sample}", "alignment", "{target}",
                              "picard", "{sample}.{target}.insert_size_histogram.pdf"),
        aln_summary    = join(WORKPATH, "{sample}", "alignment", "{target}",
                              "picard", "{sample}.{target}.alignment_summary_metrics"),
    params:
        rname      = "picard_collect_metrics",
        prefix     = join(WORKPATH, "{sample}", "alignment", "{target}",
                          "picard", "{sample}.{target}"),
        extra      = config["parameters"]["alignment"]["picard_collect_metrics"],
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.picard_collect_metrics.log"),
    resources:
        partition = allocated("partition", "picard_collect_metrics", cluster),
        mem       = allocated("mem",       "picard_collect_metrics", cluster),
        time      = allocated("time",      "picard_collect_metrics", cluster),
    threads:
        int(allocated("threads", "picard_collect_metrics", cluster))
    container:
        config["images"]["picard"]
    shell: """
set -euo pipefail
mkdir -p "$(dirname "{output.insert_metrics}")"

# Picard's behaviour on a BAM with few or no aligned reads is not uniform:
# with no *paired* mapped reads it exits 0 but silently skips the insert-size
# outputs, and with no mapped reads at all it exits non-zero. Both happen
# legitimately in a multi-target run - a library barely maps to a divergent
# reference - so neither should take the run down. The declared outputs are
# touched either way, and the failure is recorded rather than propagated.
if ! picard CollectMultipleMetrics {params.extra} \
    I="{input.bam}" O="{params.prefix}" R="{input.fa}" \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    >> "{log}" 2>&1; then
    echo "Picard CollectMultipleMetrics found nothing to report for {wildcards.sample} / {wildcards.target}; see run_summary.tsv for the mapping rate." | tee -a "{log}"
fi

touch "{output.insert_metrics}" "{output.insert_hist}" "{output.aln_summary}"
"""


