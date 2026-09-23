# report.smk: cross-sample aggregation written straight into final_report/{target}/.
# Helper scripts are declared as inputs, so editing one reruns the rule that uses it.

from os.path import join
from scripts.common import allocated

_FR = join(WORKPATH, "final_report")   # shorthand used throughout


# ── bcftools_merge / snpeff_aggregate / bgzip_aggregate /
#    gatk_aggregate_table / collect_mapping_summary ──────────────────────────

rule bcftools_merge:
    """Merge the per-sample annotated VCFs for a target, stripping FORMAT/GL and PL,
    which htsjdk cannot parse once the merge leaves them partially missing."""
    input:
        vcfs = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.freebayes_haplotypecaller.snpeff.vcf.gz"),
            sample=SAMPLES,
        ),
        tbis = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.freebayes_haplotypecaller.snpeff.vcf.gz.tbi"),
            sample=SAMPLES,
        ),
    output:
        agg_vcf = join(_FR, "{target}", "aggregate.{target}.vcf.gz"),
        agg_tbi = join(_FR, "{target}", "aggregate.{target}.vcf.gz.tbi"),
    params:
        rname       = "bcftools_merge",
        extra_merge = config["parameters"]["report"]["bcftools_merge"],
        outdir      = join(_FR, "{target}"),
    log:
        join(WORKPATH, "logfiles", "report", "bcftools_merge.{target}.log"),
    resources:
        partition = allocated("partition", "bcftools_merge", cluster),
        mem       = allocated("mem",       "bcftools_merge", cluster),
        time      = allocated("time",      "bcftools_merge", cluster),
    threads:
        int(allocated("threads", "bcftools_merge", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

NVCF=$(echo {input.vcfs} | wc -w)
if [ "$NVCF" -gt 1 ]; then
    # Indexed merge in bcftools 1.21 segfaults when no input has a record (a target
    # nothing mapped to); every VCF shares one reference order, so stream instead.
    bcftools merge --no-index {params.extra_merge} --threads {threads} \
        -O u {input.vcfs} 2>> "{log}" \
    | bcftools annotate -x FORMAT/GL,FORMAT/PL \
        -O z -o "{output.agg_vcf}" 2>> "{log}"
else
    bcftools annotate -x FORMAT/GL,FORMAT/PL \
        -O z -o "{output.agg_vcf}" {input.vcfs} >> "{log}" 2>&1
fi
tabix -p vcf "{output.agg_vcf}" >> "{log}" 2>&1
"""


rule snpeff_aggregate:
    """Re-annotate the merged VCF with SnpEff, emitting plain VCF for bgzip_aggregate."""
    input:
        vcf = rules.bcftools_merge.output.agg_vcf,
        snpeff_cfg = join(WORKPATH, "ref_db", "{target}", "snpEff.config"),
    output:
        vcf = temp(join(_FR, "{target}", "aggregate.{target}.snpeff.vcf")),
        snpeff_sum   = join(_FR, "{target}", "aggregate.{target}.snpEff_summary.html"),
        snpeff_genes = join(_FR, "{target}", "aggregate.{target}.snpEff_summary.genes.txt"),
    params:
        rname  = "snpeff_aggregate",
        target = "{target}",
        extra  = config["parameters"]["report"]["snpeff"],
    log:
        join(WORKPATH, "logfiles", "report", "snpeff_aggregate.{target}.log"),
    resources:
        partition = allocated("partition", "snpeff_aggregate", cluster),
        mem       = allocated("mem",       "snpeff_aggregate", cluster),
        time      = allocated("time",      "snpeff_aggregate", cluster),
    threads:
        int(allocated("threads", "snpeff_aggregate", cluster))
    container:
        config["images"]["snpeff"]
    shell: """
set -euo pipefail
snpEff ann {params.extra} \
    -config "{input.snpeff_cfg}" \
    -stats "{output.snpeff_sum}" \
    "{params.target}" "{input.vcf}" > "{output.vcf}" 2>> "{log}"
touch "{output.snpeff_sum}" "{output.snpeff_genes}"
"""


rule bgzip_aggregate:
    """Compress and index the annotated aggregate VCF."""
    input:
        vcf = rules.snpeff_aggregate.output.vcf,
    output:
        ann_vcf = join(_FR, "{target}", "aggregate.{target}.snpeff.vcf.gz"),
        ann_tbi = join(_FR, "{target}", "aggregate.{target}.snpeff.vcf.gz.tbi"),
    params:
        rname = "bgzip_aggregate",
    log:
        join(WORKPATH, "logfiles", "report", "bgzip_aggregate.{target}.log"),
    resources:
        partition = allocated("partition", "bgzip_aggregate", cluster),
        mem       = allocated("mem",       "bgzip_aggregate", cluster),
        time      = allocated("time",      "bgzip_aggregate", cluster),
    threads:
        int(allocated("threads", "bgzip_aggregate", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
bgzip -c "{input.vcf}" > "{output.ann_vcf}" 2>> "{log}"
tabix -p vcf "{output.ann_vcf}" >> "{log}" 2>&1
"""


rule gatk_aggregate_table:
    """Flatten the aggregate annotated VCF into a table. DP is the cross-sample sum;
    per-sample depth is in the .DP columns."""
    input:
        vcf = rules.bgzip_aggregate.output.ann_vcf,
        tbi = rules.bgzip_aggregate.output.ann_tbi,
        fa  = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        variants_tbl = join(_FR, "{target}", "aggregate.{target}.snpeff.variants.txt"),
    params:
        rname = "gatk_aggregate_table",
        extra = config["parameters"]["report"]["gatk_variants_to_table"],
    log:
        join(WORKPATH, "logfiles", "report", "gatk_aggregate_table.{target}.log"),
    resources:
        partition = allocated("partition", "gatk_aggregate_table", cluster),
        mem       = allocated("mem",       "gatk_aggregate_table", cluster),
        time      = allocated("time",      "gatk_aggregate_table", cluster),
    threads:
        int(allocated("threads", "gatk_aggregate_table", cluster))
    container:
        config["images"]["gatk4"]
    shell: """
set -euo pipefail
gatk VariantsToTable {params.extra} \
    -V "{input.vcf}" -R "{input.fa}" \
    -O "{output.variants_tbl}" >> "{log}" 2>&1
"""


rule collect_mapping_summary:
    """Per-sample mapping summary for a target, from the pre-filter flagstat (the
    filtered BAM always reads about 100% mapped)."""
    input:
        flagstats = expand(
            join(WORKPATH, "{sample}", "alignment", "{{target}}",
                 "{sample}.{{target}}.bowtie2_map.raw.flagstat"),
            sample=SAMPLES,
        ),
        script = join(WORKPATH, "workflow", "scripts", "collect_report.py"),
    output:
        summary_tsv = join(_FR, "{target}", "aggregate.{target}.mapping_summary.tsv"),
    params:
        rname  = "collect_mapping_summary",
        target = "{target}",
    log:
        join(WORKPATH, "logfiles", "report", "collect_mapping_summary.{target}.log"),
    resources:
        partition = allocated("partition", "collect_mapping_summary", cluster),
        mem       = allocated("mem",       "collect_mapping_summary", cluster),
        time      = allocated("time",      "collect_mapping_summary", cluster),
    threads:
        int(allocated("threads", "collect_mapping_summary", cluster))
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail
python3 "{input.script}" \
    --flagstats {input.flagstats} \
    --output "{output.summary_tsv}" \
    --target "{params.target}" >> "{log}" 2>&1
"""


# ── quast_consensus ───────────────────────────────────────────────────────────

rule quast_consensus:
    """QUAST assessment of per-sample consensus FASTAs against the target reference."""
    input:
        fastas = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.consensus.fa"),
            sample=SAMPLES,
        ),
        ref = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        html = join(_FR, "{target}", "quast", "report.html"),
        tsv  = join(_FR, "{target}", "quast", "report.tsv"),
    params:
        rname  = "quast_consensus",
        outdir = join(_FR, "{target}", "quast"),
        extra  = config["parameters"]["report"]["quast"],
    log:
        join(WORKPATH, "logfiles", "report", "quast_consensus.{target}.log"),
    resources:
        partition = allocated("partition", "quast_consensus", cluster),
        mem       = allocated("mem",       "quast_consensus", cluster),
        time      = allocated("time",      "quast_consensus", cluster),
    threads:
        int(allocated("threads", "quast_consensus", cluster))
    container:
        config["images"]["quast"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"
# Consensus sequences that are all N after masking make QUAST exit non-zero, a
# legitimate outcome for a weak target, so it is logged and stubbed.
INFORMATIVE=0
for fa in {input.fastas}; do
    if grep -v "^>" "$fa" | tr -d "\n" | tr -d "Nn" | grep -q .; then
        INFORMATIVE=1; break
    fi
done

if [ "$INFORMATIVE" -eq 0 ]; then
    echo "All consensus sequences for {wildcards.target} are entirely N; skipping QUAST." | tee -a "{log}"
    printf 'Assembly\tall_N\n#_contigs\t0\nNote\tno informative consensus for {wildcards.target}\n' > "{output.tsv}"
    printf '<html><body><h2>QUAST skipped: {wildcards.target}</h2><p>Every consensus sequence was entirely N after depth masking, so there was nothing to assess. See run_summary.tsv for mapping rate and coverage.</p></body></html>\n' > "{output.html}"
else
    quast.py {params.extra} \
        --threads {threads} \
        -r "{input.ref}" \
        -o "{params.outdir}" \
        {input.fastas} >> "{log}" 2>&1
fi
"""


# ── make_variants_long_table ──────────────────────────────────────────────────

rule make_variants_long_table:
    """Reshape the aggregate variants table into a long TSV, one row per sample x variant."""
    input:
        tbl = join(_FR, "{target}", "aggregate.{target}.snpeff.variants.txt"),
        script = join(WORKPATH, "workflow", "scripts", "make_variants_long_table.py"),
    output:
        long_tbl = join(_FR, "{target}", "aggregate.{target}.variants_long.tsv"),
    params:
        rname   = "make_variants_long_table",
        target  = "{target}",
        samples = SAMPLES,
    log:
        join(WORKPATH, "logfiles", "report", "make_variants_long_table.{target}.log"),
    resources:
        partition = allocated("partition", "make_variants_long_table", cluster),
        mem       = allocated("mem",       "make_variants_long_table", cluster),
        time      = allocated("time",      "make_variants_long_table", cluster),
    threads: 1
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail
python3 "{input.script}" \
    --input   "{input.tbl}" \
    --output  "{output.long_tbl}" \
    --target  "{params.target}" \
    --samples {params.samples} \
    >> "{log}" 2>&1
"""


# ── make_variants_matrix ──────────────────────────────────────────────────────

rule make_variants_matrix:
    """Variant x sample matrix (AD and percentage per sample), built for the raw and
    filtered sets. Needs -F TYPE -F ANN in gatk_variants_to_table."""
    input:
        tbl = lambda wc: join(
            _FR, wc.target,
            "aggregate.%s.snpeff.variants.txt" % wc.target if wc.vset == "raw"
            else "aggregate.%s.filtered.variants.txt" % wc.target),
        script = join(WORKPATH, "workflow", "scripts", "make_variants_matrix.py"),
    output:
        matrix = join(_FR, "{target}",
                      "aggregate.{target}.variants_matrix.{vset}.tsv"),
    params:
        rname   = "make_variants_matrix",
        samples = SAMPLES,
    log:
        join(WORKPATH, "logfiles", "report",
             "make_variants_matrix.{target}.{vset}.log"),
    resources:
        partition = allocated("partition", "make_variants_matrix", cluster),
        mem       = allocated("mem",       "make_variants_matrix", cluster),
        time      = allocated("time",      "make_variants_matrix", cluster),
    threads: 1
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail
python3 "{input.script}" \
    --input   "{input.tbl}" \
    --output  "{output.matrix}" \
    --samples {params.samples} \
    >> "{log}" 2>&1
"""


# ── bcftools_merge_filtered ───────────────────────────────────────────────────

rule bcftools_merge_filtered:
    """Merge the per-sample SnpSift-filtered VCFs into one aggregate, stripping
    FORMAT/GL and PL as bcftools_merge does."""
    input:
        vcfs = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.snpsift_filtered.vcf.gz"),
            sample=SAMPLES,
        ),
        tbis = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.snpsift_filtered.vcf.gz.tbi"),
            sample=SAMPLES,
        ),
    output:
        vcf = join(_FR, "{target}", "aggregate.{target}.filtered.vcf.gz"),
        tbi = join(_FR, "{target}", "aggregate.{target}.filtered.vcf.gz.tbi"),
    params:
        rname       = "bcftools_merge_filtered",
        outdir      = join(_FR, "{target}"),
        extra_merge = config["parameters"]["report"]["bcftools_merge"],
    log:
        join(WORKPATH, "logfiles", "report", "bcftools_merge_filtered.{target}.log"),
    resources:
        partition = allocated("partition", "bcftools_merge_filtered", cluster),
        mem       = allocated("mem",       "bcftools_merge_filtered", cluster),
        time      = allocated("time",      "bcftools_merge_filtered", cluster),
    threads:
        int(allocated("threads", "bcftools_merge_filtered", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

NVCF=$(echo {input.vcfs} | wc -w)
if [ "$NVCF" -gt 1 ]; then
    # Indexed merge in bcftools 1.21 segfaults when no input has a record (a target
    # nothing mapped to); every VCF shares one reference order, so stream instead.
    bcftools merge --no-index {params.extra_merge} --threads {threads} \
        -O u {input.vcfs} 2>> "{log}" \
    | bcftools annotate -x FORMAT/GL,FORMAT/PL \
        -O z -o "{output.vcf}" 2>> "{log}"
else
    bcftools annotate -x FORMAT/GL,FORMAT/PL \
        -O z -o "{output.vcf}" {input.vcfs} >> "{log}" 2>&1
fi
tabix -p vcf "{output.vcf}" >> "{log}" 2>&1
"""


# ── gatk_aggregate_table_filtered ─────────────────────────────────────────────

rule gatk_aggregate_table_filtered:
    """Flatten the aggregate filtered VCF into a table with the same columns as
    gatk_aggregate_table."""
    input:
        vcf = rules.bcftools_merge_filtered.output.vcf,
        tbi = rules.bcftools_merge_filtered.output.tbi,
        fa  = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        tbl = join(_FR, "{target}", "aggregate.{target}.filtered.variants.txt"),
    params:
        rname = "gatk_aggregate_table_filtered",
        extra = config["parameters"]["report"]["gatk_variants_to_table"],
    log:
        join(WORKPATH, "logfiles", "report",
             "gatk_aggregate_table_filtered.{target}.log"),
    resources:
        partition = allocated("partition", "gatk_aggregate_table_filtered", cluster),
        mem       = allocated("mem",       "gatk_aggregate_table_filtered", cluster),
        time      = allocated("time",      "gatk_aggregate_table_filtered", cluster),
    threads:
        int(allocated("threads", "gatk_aggregate_table_filtered", cluster))
    container:
        config["images"]["gatk4"]
    shell: """
set -euo pipefail
gatk VariantsToTable {params.extra} \
    -V "{input.vcf}" -R "{input.fa}" \
    -O "{output.tbl}" >> "{log}" 2>&1
"""


# ── plot_report_figures ───────────────────────────────────────────────────────

rule plot_report_figures:
    """Static per-target figures for readers who will not open a genome browser, drawn
    only from files the pipeline already wrote."""
    input:
        matrix = join(_FR, "{target}", "aggregate.{target}.variants_matrix.filtered.tsv"),
        comps  = expand(
            join(WORKPATH, "{sample}", "pre_process", "kraken2",
                 "{sample}.kraken2_decon.composition.tsv"),
            sample=SAMPLES,
        ),
        script = join(WORKPATH, "workflow", "scripts", "plot_report_figures.py"),
    output:
        done = join(_FR, "{target}", "figures", ".figures_done"),
    params:
        rname     = "plot_report_figures",
        outdir    = join(_FR, "{target}", "figures"),
        workpath  = WORKPATH,
        samples   = SAMPLES,
        target    = "{target}",
        min_depth = config["parameters"]["variant_calling"].get("consensus_min_depth", "10"),
        taxids    = TARGET_TAXIDS,
    log:
        join(WORKPATH, "logfiles", "report", "plot_report_figures.{target}.log"),
    resources:
        partition = allocated("partition", "plot_report_figures", cluster),
        mem       = allocated("mem",       "plot_report_figures", cluster),
        time      = allocated("time",      "plot_report_figures", cluster),
    threads: 1
    container:
        config["images"]["plotting"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"
python3 "{input.script}" \
    --workpath "{params.workpath}" \
    --outdir   "{params.outdir}" \
    --target   "{params.target}" \
    --samples  {params.samples} \
    --matrix   "{input.matrix}" \
    --target-taxids {params.taxids} \
    --consensus-min-depth {params.min_depth} >> "{log}" 2>&1
touch "{output.done}"
"""


# ── igv_downsample_bam ────────────────────────────────────────────────────────

rule igv_downsample_bam:
    """Downsample a sample's alignment to igv_target_depth for the IGV report only, so
    report size depends on variant count rather than sequencing depth."""
    input:
        bam = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam"),
        summ = join(WORKPATH, "{sample}", "alignment", "{target}", "mosdepth",
                    "{sample}.{target}.mosdepth.summary.txt"),
    output:
        bam = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                        "{sample}.{target}.igvshallow.bam")),
        bai = temp(join(WORKPATH, "{sample}", "alignment", "{target}",
                        "{sample}.{target}.igvshallow.bam.bai")),
    params:
        rname  = "igv_downsample_bam",
        target_depth = config["parameters"]["report"].get("igv_target_depth", "200"),
    log:
        join(WORKPATH, "logfiles", "report",
             "igv_downsample_bam.{sample}.{target}.log"),
    resources:
        partition = allocated("partition", "igv_downsample_bam", cluster),
        mem       = allocated("mem",       "igv_downsample_bam", cluster),
        time      = allocated("time",      "igv_downsample_bam", cluster),
    threads:
        int(allocated("threads", "igv_downsample_bam", cluster))
    container:
        config["images"]["samtools"]
    shell: """
set -euo pipefail

# Mean depth over the whole reference, from the mosdepth summary already built.
MEAN=$(awk -F'\\t' '$1=="total" {{print $4}}' "{input.summ}" | head -1)
[ -n "$MEAN" ] || MEAN=$(awk -F'\\t' 'NR==2 {{print $4}}' "{input.summ}")
[ -n "$MEAN" ] || MEAN=0

# Keep every read when the sample is already at or below the target depth.
FRAC=$(awk -v m="$MEAN" -v t={params.target_depth} \
    'BEGIN {{ if (m <= t || m <= 0) print 1; else printf "%.6f", t/m }}')
echo "mean depth $MEAN, target {params.target_depth}, keeping fraction $FRAC" >> "{log}"

if [ "$FRAC" = "1" ]; then
    cp "{input.bam}" "{output.bam}"
else
    samtools view -b -s "$FRAC" --threads {threads} \
        -o "{output.bam}" "{input.bam}" 2>> "{log}"
fi
samtools index -@ {threads} "{output.bam}" >> "{log}" 2>&1
"""


# ── igv_report ────────────────────────────────────────────────────────────────

rule igv_report:
    """Self-contained HTML variant report with embedded igv.js, viewable without IGV or
    access to the run directory."""
    input:
        vcf  = join(_FR, "{target}", "aggregate.{target}.filtered.vcf.gz"),
        tbi  = join(_FR, "{target}", "aggregate.{target}.filtered.vcf.gz.tbi"),
        fa   = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
        bams = expand(
            join(WORKPATH, "{sample}", "alignment", "{{target}}",
                 "{sample}.{{target}}.igvshallow.bam"),
            sample=SAMPLES,
        ),
        bais = expand(
            join(WORKPATH, "{sample}", "alignment", "{{target}}",
                 "{sample}.{{target}}.igvshallow.bam.bai"),
            sample=SAMPLES,
        ),
    output:
        html = join(_FR, "{target}", "igv_report.{target}.html"),
    params:
        rname     = "igv_report",
        flanking  = config["parameters"]["report"].get("igv_report_flanking", "150"),
        subsample = config["parameters"]["report"].get("igv_report_subsample", "100"),
    log:
        join(WORKPATH, "logfiles", "report", "igv_report.{target}.log"),
    resources:
        partition = allocated("partition", "igv_report", cluster),
        mem       = allocated("mem",       "igv_report", cluster),
        time      = allocated("time",      "igv_report", cluster),
    threads: 1
    container:
        config["images"]["igvreports"]
    shell: """
set -euo pipefail

# An empty variant set still gets a report, so a missing report means a failure.
if ! create_report "{input.vcf}" \
        --fasta "{input.fa}" \
        --flanking {params.flanking} \
        --tracks "{input.vcf}" {input.bams} \
        --output "{output.html}" >> "{log}" 2>&1; then
    echo "igv-reports failed for {wildcards.target}; see {log}" | tee -a "{log}"
    printf '<html><body><h3>IGV report unavailable for %s</h3><p>See logfiles/report/igv_report.%s.log</p></body></html>\\n' \
        "{wildcards.target}" "{wildcards.target}" > "{output.html}"
fi
"""


# ── multiqc_report ────────────────────────────────────────────────────────────

rule multiqc_report:
    """Project-wide MultiQC across all stages and samples. Inputs only order it after
    every stage, since MultiQC scans the run directory rather than a file list."""
    input:
        fastp_json = expand(
            join(WORKPATH, "{sample}", "pre_process",
                 "{sample}.fastp_trim.json"),
            sample=SAMPLES,
        ),
        fastqc_zip = expand(
            join(WORKPATH, "{sample}", "pre_process", "fastqc",
                 "{sample}.kraken2_decon.R1_fastqc.zip"),
            sample=SAMPLES,
        ),
        flagstats = expand(
            join(WORKPATH, "{sample}", "alignment", "{target}",
                 "{sample}.{target}.bowtie2_map.flagstat"),
            sample=SAMPLES, target=TARGETS,
        ),
        idxstats = expand(
            join(WORKPATH, "{sample}", "alignment", "{target}",
                 "{sample}.{target}.bowtie2_map.idxstats"),
            sample=SAMPLES, target=TARGETS,
        ),
        snpeff_sums = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{target}",
                 "{sample}.{target}.snpEff_summary.html"),
            sample=SAMPLES, target=TARGETS,
        ),
        bcf_stats = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{target}",
                 "{sample}.{target}.bcftools_stats.txt"),
            sample=SAMPLES, target=TARGETS,
        ),
        snpsift = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{target}",
                 "{sample}.{target}.snpsift_filtered.variants.txt"),
            sample=SAMPLES, target=TARGETS,
        ),
        nextclade = expand(
            join(WORKPATH, "{sample}", "lineage", "{target}",
                 "{sample}.{target}.nextclade.csv"),
            sample=SAMPLES, target=NEXTCLADE_TARGETS,
        ),
    output:
        html = join(_FR, "multiqc", "project_multiqc_report.html"),
    params:
        rname  = "multiqc_report",
        indir  = WORKPATH,
        outdir = join(_FR, "multiqc"),
        extra  = config["parameters"]["report"]["multiqc"],
        mqc_config = join(WORKPATH, "resources", "multiqc_config.yaml"),
    log:
        join(WORKPATH, "logfiles", "report", "multiqc_report.log"),
    resources:
        partition = allocated("partition", "multiqc_report", cluster),
        mem       = allocated("mem",       "multiqc_report", cluster),
        time      = allocated("time",      "multiqc_report", cluster),
    threads:
        int(allocated("threads", "multiqc_report", cluster))
    container:
        config["images"]["multiqc"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"
multiqc {params.extra} --config "{params.mqc_config}" --force \\
    --outdir "{params.outdir}" \\
    "{params.indir}" >> "{log}" 2>&1
mv "{params.outdir}/multiqc_report.html" "{output.html}"
"""


# ── make_igv_session ──────────────────────────────────────────────────────────

rule make_igv_session:
    """IGV XML session for a target: every sample BAM, annotated VCF and consensus FASTA."""
    input:
        bams = expand(
            join(WORKPATH, "{sample}", "alignment", "{{target}}",
                 "{sample}.{{target}}.bowtie2_map.bam"),
            sample=SAMPLES,
        ),
        vcfs = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.freebayes_haplotypecaller.snpeff.vcf.gz"),
            sample=SAMPLES,
        ),
        fastas = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{{target}}",
                 "{sample}.{{target}}.consensus.fa"),
            sample=SAMPLES,
        ),
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
        script = join(WORKPATH, "workflow", "scripts", "make_igv_session.py"),
    output:
        xml = join(_FR, "{target}", "igv_session.{target}.xml"),
    params:
        rname  = "make_igv_session",
        target = "{target}",
    log:
        join(WORKPATH, "logfiles", "report", "make_igv_session.{target}.log"),
    resources:
        partition = allocated("partition", "make_igv_session", cluster),
        mem       = allocated("mem",       "make_igv_session", cluster),
        time      = allocated("time",      "make_igv_session", cluster),
    threads: 1
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail
python3 "{input.script}" \\
    --target  "{params.target}" \\
    --ref     "{input.fa}" \\
    --bams    {input.bams} \\
    --vcfs    {input.vcfs} \\
    --fastas  {input.fastas} \\
    --output  "{output.xml}" >> "{log}" 2>&1
"""


# ── final_report ──────────────────────────────────────────────────────────────

rule final_report:
    """Assemble final_report/, the directory to open first: sort aggregates into
    consensus/, lineage/, variants/ and qc/, and write run_summary.tsv."""
    input:
        multiqc_html = join(_FR, "multiqc", "project_multiqc_report.html"),
        consensus = expand(
            join(WORKPATH, "{sample}", "variant_calling", "{target}",
                 "{sample}.{target}.consensus.fa"),
            sample=SAMPLES, target=TARGETS,
        ),
        composition = expand(
            join(WORKPATH, "{sample}", "pre_process", "kraken2",
                 "{sample}.kraken2_decon.composition.tsv"),
            sample=SAMPLES,
        ),
        nextclade = expand(
            join(WORKPATH, "{sample}", "lineage", "{target}",
                 "{sample}.{target}.nextclade.tsv"),
            sample=SAMPLES, target=NEXTCLADE_TARGETS,
        ),
        agg_vcfs = expand(
            join(_FR, "{target}", "aggregate.{target}.vcf.gz"),
            target=TARGETS,
        ),
        long_tbls = expand(
            join(_FR, "{target}", "aggregate.{target}.variants_long.tsv"),
            target=TARGETS,
        ),
        matrices = expand(
            join(_FR, "{target}", "aggregate.{target}.variants_matrix.{vset}.tsv"),
            target=TARGETS, vset=["raw", "filtered"],
        ),
        # Produced by its own rule since the bcftools_merge split; nothing else
        # names it, so final_report is where it has to be requested.
        mapping_summaries = expand(
            join(_FR, "{target}", "aggregate.{target}.mapping_summary.tsv"),
            target=TARGETS,
        ),
        quast_htmls = expand(
            join(_FR, "{target}", "quast", "report.html"),
            target=TARGETS,
        ),
        igv_xmls = expand(
            join(_FR, "{target}", "igv_session.{target}.xml"),
            target=TARGETS,
        ),
        figures = expand(
            join(_FR, "{target}", "figures", ".figures_done"),
            target=TARGETS,
        ),
        igv_htmls = expand(
            join(_FR, "{target}", "igv_report.{target}.html"),
            target=TARGETS,
        ),
        script = join(WORKPATH, "workflow", "scripts", "collect_final_report.py"),
    output:
        flag    = join(_FR, ".done"),
        summary = join(_FR, "run_summary.tsv"),
    params:
        rname    = "final_report",
        outdir   = _FR,
        workpath = WORKPATH,
        samples  = SAMPLES,
        targets  = TARGETS,
        lineage  = LINEAGE_TARGETS if LINEAGE_TARGETS else [],
        taxids   = ["%s=%s" % (t, (_TARGET_REFS.get(t) or {}).get("taxid", ""))
                    for t in TARGETS if (_TARGET_REFS.get(t) or {}).get("taxid")],
        min_cov  = config["parameters"]["variant_calling"].get("min_genome_coverage", "0.80"),
        min_depth= config["parameters"]["variant_calling"].get("consensus_min_depth", "10"),
    log:
        join(WORKPATH, "logfiles", "report", "final_report.log"),
    resources:
        partition = allocated("partition", "final_report", cluster),
        mem       = allocated("mem",       "final_report", cluster),
        time      = allocated("time",      "final_report", cluster),
    threads: 1
    container:
        config["images"]["python3"]
    shell: """
set -euo pipefail

OUTDIR="{params.outdir}"
WPATH="{params.workpath}"

# ── Sort aggregates, pull in per-sample results, write run_summary.tsv ─────
python3 "{input.script}" \
    --workpath "$WPATH" \
    --outdir   "$OUTDIR" \
    --samples  {params.samples} \
    --targets  {params.targets} \
    --lineage-targets {params.lineage} \
    --target-taxids {params.taxids} \
    --min-genome-coverage {params.min_cov} \
    --consensus-min-depth {params.min_depth} \
    >> "{log}" 2>&1

touch "{output.flag}"
echo "final_report complete: $OUTDIR" >> "{log}" 2>&1
"""
