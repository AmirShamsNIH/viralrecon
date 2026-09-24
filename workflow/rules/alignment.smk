# alignment.smk: map depleted reads to each target, mark duplicates and collect alignment QC.

from os.path import join
from scripts.common import allocated


rule bowtie2_map:
    """Align depleted reads to a target with Bowtie2, sort, and drop unmapped, secondary
    and supplementary records. The raw flagstat is taken before that filter."""
    input:
        r1 = join(WORKPATH, "{sample}", "pre_process",
                  "{sample}.kraken2_decon.R1.fastq.gz"),
        r2 = join(WORKPATH, "{sample}", "pre_process",
                  "{sample}.kraken2_decon.R2.fastq.gz") if PAIRED else [],
        bt2 = join(WORKPATH, "ref_db", "{target}", "{target}.1.bt2"),
    output:
        # Handed to mark_duplicates, which produces the canonical
        # {sample}.{target}.bowtie2_map.bam that everything downstream reads.
        bam = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                        "{sample}.{target}.aligned.bam")),
        bai = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                        "{sample}.{target}.aligned.bam.bai")),
        raw_stat = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "{sample}.{target}.bowtie2_map.raw.flagstat"),
    params:
        rname = "bowtie2_map",
        paired = PAIRED,
        extra_bt2 = config["parameters"]["alignment"]["bowtie2_map"],
        extra_picard = config["parameters"]["alignment"]["picard_markdup"],
        filt = config["parameters"]["alignment"].get("samtools_filter", "-F 0x004"),
        min_mapped = config["parameters"]["alignment"].get("min_mapped_reads", "1000"),
        idx = join(WORKPATH, "ref_db", "{target}", "{target}"),
        outdir = join(WORKPATH, "{sample}", "alignment", "{target}"),
        tmpbam = join(WORKPATH, "{sample}", "alignment", "{target}",
                      "{sample}.{target}.raw.bam"),
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.bowtie2_map.log"),
    resources:
        partition = allocated("partition", "bowtie2_map", cluster),
        mem = allocated("mem", "bowtie2_map", cluster),
        time = allocated("time", "bowtie2_map", cluster),
    threads:
        int(allocated("threads", "bowtie2_map", cluster))
    container:
        config["images"]["bowtie2"]
    shell: """
set -euo pipefail

mkdir -p "$(dirname "{output.bam}")"

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

# The pre-filter flagstat is the only view of the real mapping rate.
samtools flagstat --threads {threads} "{params.tmpbam}" > "{output.raw_stat}" 2>> "{log}"

# Too few mapped reads is a warning, never a reason to abort the run.
MAPPED=$(awk '/ mapped \(/ && !/primary/ {{print $1; exit}}' "{output.raw_stat}")
MAPPED=${{MAPPED:-0}}
if [ "$MAPPED" -lt "{params.min_mapped}" ]; then
    echo "QC WARNING: {wildcards.sample} / {wildcards.target}: $MAPPED mapped reads < min_mapped_reads={params.min_mapped}" \\
        | tee -a "{log}" > "{params.outdir}/{wildcards.sample}.{wildcards.target}.qc_warn"
fi

samtools view -b {params.filt} -@ {threads} \
    -o "{output.bam}" "{params.tmpbam}" >> "{log}" 2>&1
samtools index "{output.bam}" >> "{log}" 2>&1
rm -f "{params.tmpbam}" "{params.tmpbam}.bai"
"""


rule mark_duplicates:
    """Produce the canonical alignment BAM, marking duplicates only when
    skip_markduplicates=false (off by default: nearly every viral read is a duplicate)."""
    input:
        bam = rules.bowtie2_map.output.bam,
        bai = rules.bowtie2_map.output.bai,
    output:
        bam = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam"),
        bai = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam.bai"),
        dupmetric = join(WORKPATH, "{sample}", "alignment", "{target}",
                         "{sample}.{target}.bowtie2_map.dup_metrics.txt"),
    params:
        rname = "mark_duplicates",
        extra_picard = config["parameters"]["alignment"]["picard_markdup"],
        skip_markdup = config["parameters"]["alignment"].get("skip_markduplicates", "true"),
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.mark_duplicates.log"),
    resources:
        partition = allocated("partition", "mark_duplicates", cluster),
        mem = allocated("mem", "mark_duplicates", cluster),
        time = allocated("time", "mark_duplicates", cluster),
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
    picard BuildBamIndex I="{output.bam}" O="{output.bai}" >> "{log}" 2>&1 \
        || samtools index "{output.bam}" "{output.bai}" >> "{log}" 2>&1
fi
"""


rule flagstat_align:
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
        mem = allocated("mem", "flagstat_align", cluster),
        time = allocated("time", "flagstat_align", cluster),
    threads: 2
    container:
        config["images"]["samtools"]
    shell: """
set -euo pipefail
samtools flagstat --threads {threads} "{input.bam}" > "{output.stat}" 2>> "{log}"
"""


rule idxstats_align:
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
        mem = allocated("mem", "idxstats_align", cluster),
        time = allocated("time", "idxstats_align", cluster),
    threads: 1
    container:
        config["images"]["samtools"]
    shell: """
set -euo pipefail
samtools idxstats "{input.bam}" > "{output.stat}" 2>> "{log}"
"""


rule mosdepth_align:
    """Per-base and windowed coverage depth with mosdepth."""
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bai,
    output:
        summary = join(WORKPATH, "{sample}", "alignment", "{target}",
                       "mosdepth", "{sample}.{target}.mosdepth.summary.txt"),
        global_d = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.mosdepth.global.dist.txt"),
        per_base = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.per-base.bed.gz"),
        # --by 200 windows feed the MultiQC coverage-across-genome plot.
        regions = join(WORKPATH, "{sample}", "alignment", "{target}",
                       "mosdepth", "{sample}.{target}.regions.bed.gz"),
    params:
        rname = "mosdepth_align",
        prefix = join(WORKPATH, "{sample}", "alignment", "{target}",
                      "mosdepth", "{sample}.{target}"),
        extra = config["parameters"]["alignment"]["mosdepth"],
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.mosdepth_align.log"),
    resources:
        partition = allocated("partition", "mosdepth_align", cluster),
        mem = allocated("mem", "mosdepth_align", cluster),
        time = allocated("time", "mosdepth_align", cluster),
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


rule picard_collect_metrics:
    """Picard CollectMultipleMetrics: insert size, alignment summary, GC bias."""
    input:
        bam = rules.mark_duplicates.output.bam,
        bai = rules.mark_duplicates.output.bai,
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        insert_metrics = join(WORKPATH, "{sample}", "alignment", "{target}",
                              "picard", "{sample}.{target}.insert_size_metrics"),
        insert_hist = join(WORKPATH, "{sample}", "alignment", "{target}",
                           "picard", "{sample}.{target}.insert_size_histogram.pdf"),
        aln_summary = join(WORKPATH, "{sample}", "alignment", "{target}",
                           "picard", "{sample}.{target}.alignment_summary_metrics"),
    params:
        rname = "picard_collect_metrics",
        prefix = join(WORKPATH, "{sample}", "alignment", "{target}",
                      "picard", "{sample}.{target}"),
        extra = config["parameters"]["alignment"]["picard_collect_metrics"],
    log:
        join(WORKPATH, "logfiles", "alignment",
             "{sample}.{target}.picard_collect_metrics.log"),
    resources:
        partition = allocated("partition", "picard_collect_metrics", cluster),
        mem = allocated("mem", "picard_collect_metrics", cluster),
        time = allocated("time", "picard_collect_metrics", cluster),
    threads:
        int(allocated("threads", "picard_collect_metrics", cluster))
    container:
        config["images"]["picard"]
    shell: """
set -euo pipefail
mkdir -p "$(dirname "{output.insert_metrics}")"

# Picard skips or fails on BAMs with few or no mapped reads, which is normal in a
# multi-target run, so outputs are touched and the failure is logged, not raised.
if ! picard CollectMultipleMetrics {params.extra} \
    I="{input.bam}" O="{params.prefix}" R="{input.fa}" \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    >> "{log}" 2>&1; then
    echo "Picard CollectMultipleMetrics found nothing to report for {wildcards.sample} / {wildcards.target}; see run_summary.tsv for the mapping rate." | tee -a "{log}"
fi

touch "{output.insert_metrics}" "{output.insert_hist}" "{output.aln_summary}"
"""


