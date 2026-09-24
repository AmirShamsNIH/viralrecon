# variant_calling.smk: FreeBayes calling, normalisation, SnpEff/SnpSift annotation, consensus.

from os.path import join
from scripts.common import allocated


# SnpEff and SnpSift images lack htslib, so bcftools-image rules compress their VCFs.

rule freebayes_call:
    """Call variants with FreeBayes per sample and target, downsampling ultra-deep
    libraries to max_reads so calling stays tractable."""
    input:
        bam = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam"),
        bai = join(WORKPATH, "{sample}", "alignment", "{target}",
                   "{sample}.{target}.bowtie2_map.bam.bai"),
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
        fai = join(WORKPATH, "ref_db", "{target}", "{target}.fa.fai"),
    output:
        vcf_raw = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.freebayes_haplotypecaller.vcf.gz"),
        tbi_raw = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.freebayes_haplotypecaller.vcf.gz.tbi"),
    params:
        rname = "freebayes_call",
        extra_fb = config["parameters"]["variant_calling"]["freebayes"],
        max_reads = int(config["parameters"]["variant_calling"].get(
            "freebayes_max_reads", "5000000")),
        outdir = join(WORKPATH, "{sample}", "variant_calling", "{target}"),
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.freebayes_call.log"),
    resources:
        partition = allocated("partition", "freebayes_call", cluster),
        mem = allocated("mem", "freebayes_call", cluster),
        time = allocated("time", "freebayes_call", cluster),
    threads:
        int(allocated("threads", "freebayes_call", cluster))
    container:
        config["images"]["freebayes"]
    shell: """
set -euo pipefail
mkdir -p "{params.outdir}"

MAPPED=$(samtools view -c -F 260 "{input.bam}" 2>> "{log}")
FRAC=$(python3 -c '
import sys
m, mx = int(sys.argv[1]), int(sys.argv[2])
f = min(1.0, mx / max(m, 1))
# samtools -s takes SEED.FRACTION; format the fraction with enough digits that
# very deep samples (f can be <1e-4) do not round down to 0 and drop all reads.
print("NO" if f >= 1.0 else "42" + ("%.6f" % max(f, 0.000001))[1:])
' "$MAPPED" {params.max_reads})
echo "Mapped: $MAPPED  max_reads: {params.max_reads}  downsample_frac: $FRAC" >> "{log}"

if [ "$FRAC" = "NO" ]; then
    freebayes {params.extra_fb} -f "{input.fa}" "{input.bam}" \
        | bgzip -c > "{output.vcf_raw}" 2>> "{log}"
else
    samtools view -h -b -s "$FRAC" "{input.bam}" \
        | freebayes {params.extra_fb} -f "{input.fa}" /dev/stdin \
        | bgzip -c > "{output.vcf_raw}" 2>> "{log}"
fi
tabix -p vcf "{output.vcf_raw}" >> "{log}" 2>&1
"""


rule bcftools_norm:
    """Left-align indels and split multi-allelic records."""
    input:
        vcf = rules.freebayes_call.output.vcf_raw,
        tbi = rules.freebayes_call.output.tbi_raw,
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        vcf_norm = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.freebayes_haplotypecaller.norm.vcf.gz"),
        tbi_norm = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.freebayes_haplotypecaller.norm.vcf.gz.tbi"),
    params:
        rname = "bcftools_norm",
        extra = config["parameters"]["variant_calling"]["bcftools_norm"],
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.bcftools_norm.log"),
    resources:
        partition = allocated("partition", "bcftools_norm", cluster),
        mem = allocated("mem", "bcftools_norm", cluster),
        time = allocated("time", "bcftools_norm", cluster),
    threads:
        int(allocated("threads", "bcftools_norm", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
bcftools norm {params.extra} \
    --fasta-ref "{input.fa}" \
    --output-type z --output "{output.vcf_norm}" \
    "{input.vcf}" >> "{log}" 2>&1
tabix -p vcf "{output.vcf_norm}" >> "{log}" 2>&1
"""


rule snpeff_annotate:
    """Annotate variant effects with SnpEff, emitting plain VCF for bgzip_snpeff."""
    input:
        vcf = rules.bcftools_norm.output.vcf_norm,
        snpeff_cfg = join(WORKPATH, "ref_db", "{target}", "snpEff.config"),
    output:
        vcf = temp(join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.freebayes_haplotypecaller.snpeff.vcf")),
        snpeff_sum = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                          "{sample}.{target}.snpEff_summary.html"),
    params:
        rname = "snpeff_annotate",
        target = "{target}",
        extra = config["parameters"]["variant_calling"]["snpeff"],
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.snpeff_annotate.log"),
    resources:
        partition = allocated("partition", "snpeff_annotate", cluster),
        mem = allocated("mem", "snpeff_annotate", cluster),
        time = allocated("time", "snpeff_annotate", cluster),
    threads:
        int(allocated("threads", "snpeff_annotate", cluster))
    container:
        config["images"]["snpeff"]
    shell: """
set -euo pipefail
snpEff ann {params.extra} \
    -config "{input.snpeff_cfg}" \
    -stats "{output.snpeff_sum}" \
    "{params.target}" "{input.vcf}" > "{output.vcf}" 2>> "{log}"
# snpEff skips its stats files when the VCF holds no variants
touch "{output.snpeff_sum}"
"""


rule bgzip_snpeff:
    """Compress and index the SnpEff output, the canonical annotated VCF."""
    input:
        vcf = rules.snpeff_annotate.output.vcf,
    output:
        vcf_ann = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.freebayes_haplotypecaller.snpeff.vcf.gz"),
        tbi_ann = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                       "{sample}.{target}.freebayes_haplotypecaller.snpeff.vcf.gz.tbi"),
    params:
        rname = "bgzip_snpeff",
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.bgzip_snpeff.log"),
    resources:
        partition = allocated("partition", "bgzip_snpeff", cluster),
        mem = allocated("mem", "bgzip_snpeff", cluster),
        time = allocated("time", "bgzip_snpeff", cluster),
    threads:
        int(allocated("threads", "bgzip_snpeff", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
bgzip -c "{input.vcf}" > "{output.vcf_ann}" 2>> "{log}"
tabix -p vcf "{output.vcf_ann}" >> "{log}" 2>&1
"""


rule gatk_variants_table:
    input:
        vcf = rules.bgzip_snpeff.output.vcf_ann,
        tbi = rules.bgzip_snpeff.output.tbi_ann,
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        tbl = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                   "{sample}.{target}.freebayes_haplotypecaller.snpeff.variants.txt"),
    params:
        rname = "gatk_variants_table",
        extra = config["parameters"]["variant_calling"]["gatk_variants_to_table"],
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.gatk_variants_table.log"),
    resources:
        partition = allocated("partition", "gatk_variants_table", cluster),
        mem = allocated("mem", "gatk_variants_table", cluster),
        time = allocated("time", "gatk_variants_table", cluster),
    threads:
        int(allocated("threads", "gatk_variants_table", cluster))
    container:
        config["images"]["gatk4"]
    shell: """
set -euo pipefail
gatk VariantsToTable {params.extra} \
    -V "{input.vcf}" -R "{input.fa}" \
    -O "{output.tbl}" >> "{log}" 2>&1
"""


rule bcftools_stats:
    input:
        vcf = rules.bgzip_snpeff.output.vcf_ann,
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
    output:
        stats = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                     "{sample}.{target}.bcftools_stats.txt"),
    params:
        rname = "bcftools_stats",
        extra = config["parameters"]["variant_calling"]["bcftools_stats"],
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.bcftools_stats.log"),
    resources:
        partition = allocated("partition", "bcftools_stats", cluster),
        mem = allocated("mem", "bcftools_stats", cluster),
        time = allocated("time", "bcftools_stats", cluster),
    threads:
        int(allocated("threads", "bcftools_stats", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
bcftools stats {params.extra} --fasta-ref "{input.fa}" \
    "{input.vcf}" > "{output.stats}" 2>> "{log}"
"""


rule snpsift_annotate:
    """SnpSift filter on the annotated VCF plus a flat field extraction, as plain VCF."""
    input:
        vcf = rules.bgzip_snpeff.output.vcf_ann,
        tbi = rules.bgzip_snpeff.output.tbi_ann,
    output:
        vcf_filt = temp(join(WORKPATH, "{sample}", "variant_calling", "{target}",
                             "{sample}.{target}.snpsift_filtered.vcf")),
        txt_filt = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.snpsift_filtered.variants.txt"),
    params:
        rname = "snpsift_annotate",
        filter_exp = config["parameters"]["variant_calling"]["snpsift_filter"],
        fields = config["parameters"]["variant_calling"]["snpsift_fields"],
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.snpsift_annotate.log"),
    resources:
        partition = allocated("partition", "snpsift_annotate", cluster),
        mem = allocated("mem", "snpsift_annotate", cluster),
        time = allocated("time", "snpsift_annotate", cluster),
    threads:
        int(allocated("threads", "snpsift_annotate", cluster))
    container:
        config["images"]["snpsift"]
    shell: """
set -euo pipefail
SnpSift filter "{params.filter_exp}" "{input.vcf}" > "{output.vcf_filt}" 2>> "{log}"
SnpSift extractFields -s "," -e "." \
    "{output.vcf_filt}" {params.fields} \
    > "{output.txt_filt}" 2>> "{log}"
"""


rule bgzip_snpsift:
    """Compress and index the SnpSift-filtered VCF for bcftools consensus."""
    input:
        vcf = rules.snpsift_annotate.output.vcf_filt,
    output:
        vcf_filt = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.snpsift_filtered.vcf.gz"),
        tbi_filt = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                        "{sample}.{target}.snpsift_filtered.vcf.gz.tbi"),
    params:
        rname = "bgzip_snpsift",
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.bgzip_snpsift.log"),
    resources:
        partition = allocated("partition", "bgzip_snpsift", cluster),
        mem = allocated("mem", "bgzip_snpsift", cluster),
        time = allocated("time", "bgzip_snpsift", cluster),
    threads:
        int(allocated("threads", "bgzip_snpsift", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
bgzip -c "{input.vcf}" > "{output.vcf_filt}" 2>> "{log}"
tabix -p vcf "{output.vcf_filt}" >> "{log}" 2>&1
"""


rule bcftools_consensus:
    """Per-sample consensus FASTA from the filtered variants, masking positions below
    consensus_min_depth to N so uncovered sites are not called reference."""
    input:
        vcf = rules.bgzip_snpsift.output.vcf_filt,
        tbi = rules.bgzip_snpsift.output.tbi_filt,
        fa = join(WORKPATH, "ref_db", "{target}", "{target}.fa"),
        per_base = join(WORKPATH, "{sample}", "alignment", "{target}",
                        "mosdepth", "{sample}.{target}.per-base.bed.gz"),
    output:
        consensus = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                         "{sample}.{target}.consensus.fa"),
        mask = join(WORKPATH, "{sample}", "variant_calling", "{target}",
                    "{sample}.{target}.lowcov_mask.bed"),
    params:
        rname = "bcftools_consensus",
        sample = "{sample}",
        target = "{target}",
        extra = config["parameters"]["variant_calling"]["bcftools_consensus"],
        min_depth = config["parameters"]["variant_calling"].get(
            "consensus_min_depth", "10"),
    log:
        join(WORKPATH, "logfiles", "variant_calling",
             "{sample}.{target}.bcftools_consensus.log"),
    resources:
        partition = allocated("partition", "bcftools_consensus", cluster),
        mem = allocated("mem", "bcftools_consensus", cluster),
        time = allocated("time", "bcftools_consensus", cluster),
    threads:
        int(allocated("threads", "bcftools_consensus", cluster))
    container:
        config["images"]["bcftools"]
    shell: """
set -euo pipefail
zcat "{input.per_base}" \
    | awk -v d={params.min_depth} 'BEGIN{{OFS="\t"}} $4 < d {{print $1,$2,$3}}' \
    > "{output.mask}" 2>> "{log}"
echo "masked $(awk '{{s+=$3-$2}} END{{print s+0}}' "{output.mask}") bases below {params.min_depth}x" >> "{log}"

bcftools consensus {params.extra} \
    --fasta-ref "{input.fa}" \
    --sample "{params.sample}" \
    --mask "{output.mask}" --mask-with N \
    --output "{output.consensus}" \
    "{input.vcf}" >> "{log}" 2>&1

# Header is sample.target, the key other tools derive for MultiQC, so each sample is
# one row. Snakemake formats shell comments too, so avoid braces here.
HEADER="{wildcards.sample}.{wildcards.target}"
sed -i "s/^>.*/>$HEADER/" "{output.consensus}" >> "{log}" 2>&1
"""


