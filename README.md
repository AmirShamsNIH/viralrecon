<div align="center">

<h1>viralrecon 🦠</h1>
<b>Viral Variation Analysis Pipeline</b><br/>
<i>NIH Biowulf · GRS BigSky · <a href="https://github.com/OpenOmics/baseline">OpenOmics/baseline</a> structure</i>

<br/><br/>

<a href="#"><img alt="Version" src="https://img.shields.io/badge/version-0.1.0-blue"></a>
<a href="#22-dependencies"><img alt="Snakemake" src="https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen"></a>
<a href="#23-containers"><img alt="Singularity" src="https://img.shields.io/badge/singularity-only-%23663399"></a>
<a href=".github/workflows/main.yaml"><img alt="CI" src="https://img.shields.io/badge/CI-dry--run%20DAG-lightgrey"></a>
<a href="#"><img alt="Platform" src="https://img.shields.io/badge/platform-Biowulf%20%7C%20BigSky-orange"></a>

<h3>
<a href="#1-introduction">Introduction</a> ·
<a href="#2-overview-of-the-pipeline">Overview</a> ·
<a href="#3-run-the-pipeline">Run</a> ·
<a href="#4-output">Output</a> ·
<a href="#5-contribute">Contribute</a> ·
<a href="#6-references">References</a>
</h3>

</div>

---

## 1. Introduction

Sequencing a virus directly from a clinical or environmental sample answers three
different questions at once: **what is in the sample**, **what does its genome look
like**, and **which known lineage is it**. Answering them separately, with ad-hoc
scripts per project, is where reproducibility is usually lost.

**viralrecon** answers all three in one pass. It takes raw Illumina FASTQ files and a
set of reference accessions, and produces per-sample consensus genomes, annotated
variant calls (including sub-consensus, intra-host variants), contamination profiles,
and lineage assignments — with a single directory a user is meant to open first.

Two design commitments make its results reproducible:

- **Every external tool comes from a pinned Singularity image.** No `module load`
  appears anywhere in the workflow, or in the reference-building CLI. A cluster module
  can be upgraded or removed underneath a pipeline and its version is not recorded in
  the run directory; a pinned `.sif` path is. This matters most for tools that bundle a
  database — pangolin, nextclade, freyja — where a floating version silently changes
  lineage assignments.
- **The pipeline is virus-agnostic.** Nothing outside the optional lineage stage assumes
  SARS-CoV-2. Targets are selected by accession, never by name, and a run may carry
  several unrelated viruses at once; every post-alignment file is namespaced
  `{sample}.{target}.*` so two targets cannot mix.

Orchestration is [Snakemake][1], which handles the job DAG, SLURM submission, restarts,
and container invocation.

---

## 2. Overview of the pipeline

### 2.1 Stages

Stages run in the order below. `lineage`, `amplicon`, and `assembly` are optional.

| # | Stage | Purpose | Key tools |
|---|-------|---------|-----------|
| 1 | **build_environment** | Per-run reference setup | samtools |
| 2 | **pre_process** | Repair, trim, profile, deplete | BBTools · [fastp][2] · [Kraken2][3] · [Krona][4] · [FastQC][5] |
| 3 | **alignment** | Map to each target, QC the mapping | [Bowtie2][6] · [samtools][7] · [Picard][8] · [mosdepth][9] |
| 4 | **variant_calling** | Call, normalise, annotate, build consensus | [FreeBayes][10] · [bcftools][7] · [SnpEff/SnpSift][11] · [GATK4][12] |
| 5 | **report** | Aggregate across samples, assemble `final_report/` | bcftools · GATK4 · [QUAST][13] · [MultiQC][14] |
| — | **lineage** *(opt.)* | Lineage / clade / mixture calls | [pangolin][15] · [Nextclade][16] · [Freyja][17] |

#### Quality control and decontamination

- **BBTools `reformat`** repairs the raw pairs first — junk bases, IUPAC codes, broken
  reads, and out-of-range quality scores are fixed or dropped before anything reads them.
  It does *not* tolerate a truncated `.gz`; a failure here means the input transfer, not
  the data.
- **fastp** trims adapters (auto-detected per pair), low-complexity reads, and low-quality
  tails.
- **Kraken2** classifies the trimmed reads against the standard database **before**
  depletion, so contamination is visible rather than silently discarded. **Krona** renders
  that profile as an interactive HTML chart.
- **Depletion is subtractive by taxon**, not "keep everything unclassified": reads
  assigned anywhere in the subtree of `9606` (human), `10847` (phiX), or `2` (bacteria)
  are removed, and everything else — including unclassified reads, which is where a novel
  or divergent virus lives — is kept. Taxids are configurable.
- **FastQC** runs on the depleted reads, and a single project-wide **MultiQC** report
  aggregates every stage across every sample into `final_report/multiqc/`.

#### Alignment and variant calling

- **Bowtie2** (`--local`) maps the depleted reads to each target independently. A target
  whose mapping falls below `min_mapped_reads` raises a QC warning and its downstream
  stages are skipped — **it does not fail the run**, so one weak reference cannot sink a
  multi-target analysis.
- **samtools** filters unmapped, secondary, and supplementary alignments; **Picard**
  collects insert-size and alignment metrics; **mosdepth** produces per-base and windowed
  depth.
- **FreeBayes** runs in `--pooled-continuous` mode. This is the deliberate choice for viral
  work: a sample is a *population*, not a diploid individual, so variants are called on
  allele fraction rather than genotype. The default floor of `--min-alternate-fraction
  0.02` surfaces intra-host minor variants that a genotype-based caller would discard.
- **bcftools norm** splits multi-allelic records; **SnpEff** and **SnpSift** annotate;
  **GATK4** flattens calls into a table, which is emitted in two shapes: a long table
  (one row per sample x variant, for scripting) and a **variant x sample matrix** — one
  row per variant with an allele-depth and percentage column per sample, plus parsed
  effect, severity, gene and amino-acid change. The matrix is what shows at a glance
  whether a variant is fixed across the run or a minor allele in a few samples; an empty
  cell means not called in that sample, and is left empty rather than zeroed.
- **Both a raw and a filtered set are carried to the aggregate**, as
  `variants_matrix.raw.tsv` and `variants_matrix.filtered.tsv`. The filtered set applies
  the SnpSift expression and is the one to report; the raw set is what FreeBayes called
  under its own thresholds, and is what answers "was this variant never called, or called
  and then filtered out?" — a question the filtered table alone cannot distinguish.
- **Consensus is depth-masked.** Positions below `consensus_min_depth` are written as `N`
  rather than inheriting the reference base — without this, a consensus fabricates
  reference sequence in regions with no evidence.

#### Lineage

Each of the three callers is gated by what it can actually describe, because they are not
equivalent:

- **pangolin** — Pango nomenclature exists only for SARS-CoV-2 and its database ships
  inside the image, so `pangolin_targets` is simply a list of accessions.
- **Nextclade** — dataset-driven, and the dataset belongs to the reference. Fetch it at
  build time with `--nextclade-dataset`; the path is recorded in `genome.json` and
  nextclade runs on any target that has one. A target without a dataset is skipped rather
  than described against another virus's dataset, which would return confident nonsense
  instead of failing.

- **Freyja** (mixture deconvolution, which detects co-infection or multi-strain samples a
  single consensus hides) is gated by `freyja_targets`, and that is a coordinate
  constraint rather than a preference: freyja reads a BAM in *reference* coordinates and
  matches it against barcodes keyed to Wuhan-Hu-1, so against any other reference its
  answer is meaningless even when it does not crash.

Agreement between the callers is the point; disagreement is a finding.

### 2.2 Dependencies

| Requirement | Version |
|---|---|
| Snakemake | ≥ 7.0 |
| Singularity / Apptainer | ≥ 3.5 |
| Python | ≥ 3.8 |

Input is paired or single-end Illumina FASTQ; the layout is detected from whether any
R2 file is present, not configured. Single-end has been verified against the same
samples run paired: identical lineage calls and fixed markers, and sub-consensus
frequencies within about a point, at proportionally lower depth.

Nothing else is needed locally — every tool is containerised.

### 2.3 Containers

Images are resolved from `config/containers.json`, drawn from two locations in this order
of preference:

1. **`/data/OpenOmics/SIFs`** — the shared lab library, used whenever it carries the tool
   at a version at least as new as the one pinned here.
2. **`/data/RTB_GRS/references/singularity`** — built by us where the shared library has
   no image, or only an older one.

A Snakemake rule may declare only **one** container, so a step chaining tools from
different images is split into one rule per tool. When adding a rule, give it a
`container:` directive — never a module.

### 2.4 Installation

```bash
git clone https://github.com/OpenOmics/viralrecon.git
cd viralrecon
./viralrecon --version
```

---

### 2.5 Platform profile: Biowulf

The pipeline is developed and validated on NIH Biowulf. Everything below has to
be in place before a run; most of it already is, and the paths are recorded in
`config/containers.json` and `config/config.json` rather than discovered at
runtime.

**Access**

| need | detail |
|---|---|
| Biowulf account | with access to `/data/RTB_GRS` |
| SLURM partition | `norm` — the only partition `cluster.json` uses |
| Shared references | read access to `/data/OpenOmics/SIFs` |

**Software on the submitting shell**

```bash
module load python/3.10        # provides snakemake 7.30.1
```

Snakemake is the one tool that does *not* come from a container — it is the
thing that launches the containers. The master job inherits the submitting
environment, so if `snakemake` is not on your `PATH` when you submit, the master
job fails immediately. Everything else is loaded by the job itself: the master
script runs `module load singularity` and no rule loads anything at all.

**Data that must exist on disk**

| path | what it is |
|---|---|
| `/fdb/kraken/20260226_standard_kraken2` | Kraken2 standard database, ~98 GB. Biowulf shared, nothing to install |
| `/data/RTB_GRS/references/krona/taxonomy` | Krona taxonomy, staged once. The Krona image ships a placeholder that produces empty charts |
| `/data/OpenOmics/SIFs` | shared lab image library |
| `/data/RTB_GRS/references/singularity` | images built or pulled for this pipeline |

**Which node does what** — this trips people up, because the two capabilities
live on opposite hosts:

| task | where | why |
|---|---|---|
| `viralrecon build` | **compute node** (`sbatch`) | needs `singularity`, which the login node does not have |
| `viralrecon run` | either — it only submits | the master job it creates runs on a compute node |
| `git push` | **login node** | compute nodes cannot resolve external hostnames at all |

Compute nodes have no direct internet. `viralrecon build` downloads a reference
from NCBI and a Nextclade dataset, so a build job must set the session proxy:

```bash
export http_proxy=http://dtn20-e0:3128
export https_proxy=http://dtn20-e0:3128
```

That covers HTTP and HTTPS only. SSH is not proxied, which is why a `git push`
over an SSH remote has to run from the login node.

**A complete run**, as a single batch script:

```bash
#!/usr/bin/env bash
#SBATCH --cpus-per-task=4 --mem=16g --time=2:00:00

VIRALRECON=/data/RTB_GRS/internal/pipeline/viralrecon/viralrecon
BASE=/data/RTB_GRS/IDSS_Projects/<project>
REF=$BASE/target_reference
OUT=$BASE/viralrecon_execution

module load singularity
export http_proxy=http://dtn20-e0:3128
export https_proxy=http://dtn20-e0:3128

$VIRALRECON build --virus SARS --accession NC_045512.2 --output "$REF" \
    --nextclade-dataset sars-cov-2

$VIRALRECON run --input /path/to/reads/*_R[12]_001.fastq.gz \
    --output "$OUT" --genome "$REF/genome.json" \
    --targets SARS_NC_045512.2
```

Submit it with `sbatch`, having run `module load python/3.10` first so the
`run` step can find Snakemake.

**Scratch.** The pipeline keeps its temporary files in `$OUTDIR/tmp` rather than
`lscratch`, so no `--gres=lscratch` allocation is needed. One exception is
deliberate: `pangolin_lineage` redirects `TMPDIR` to a node-local path, because
scorpio opens a Unix domain socket for `multiprocessing` and those do not work
on GPFS.

## 3. Run the pipeline

### 3.1 Build a reference

Run once per accession. References accumulate in a single `genome.json`: a new accession
is **added** to an existing build, and an accession already present is verified
file-by-file and skipped if complete, or completed if not.

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs
```

Fetch a Nextclade dataset into the reference so that target can be clade-called. The
dataset lands in `<reference>/nextclade/` and its path is recorded in `genome.json`:

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs \
                   --nextclade-dataset sars-cov-2
```

Names come from `nextclade dataset list` — `sars-cov-2`, `mpox`, `rsv_a`,
`flu_h1n1pdm_ha` and others. This step needs the proxy set.

Record curated knowledge about a reference — the kind nothing can be derived from the
files themselves. The note is stored in `genome.json`, echoed whenever a run selects that
target, and carried forward across rebuilds unless given again:

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs \
                   --notes "contains long N runs; do not use for consensus"
```

Provide your own sequence instead of downloading from NCBI — mutually exclusive with a
bare `--accession` download:

```bash
./viralrecon build --virus MYVIRUS --accession MYV_001 \
                   --fasta /local/myvirus.fa \
                   --annotation /local/myvirus.gff \
                   --output /data/refs
```

### 3.2 Run on Biowulf

```bash
./viralrecon run -i /data/fastq/*_R[12]_001.fastq.gz \
                 -o /data/out \
                 --genome /data/refs/genome.json \
                 --targets SARS_NC_045512.2
```

Preview the DAG without submitting anything:

```bash
./viralrecon run -i ... -o ... --genome ... --dry-run
```

After editing a config template in the repo, refresh it into an existing output
directory — `config/` is copied once and does not otherwise refresh:

```bash
./viralrecon run -i ... -o ... --genome ... --overwrite-pipeline-template
```

Other subcommands: `unlock` (release a stale Snakemake lock), `cache` (pre-pull images),
`install` (fetch reference data).

### 3.3 Configuration

| File | Contents |
|---|---|
| `config/config.json` | Pipeline options and every tool parameter |
| `config/cluster.json` | Per-rule SLURM resources |
| `config/genome.json` | Reference paths, written by `viralrecon build` |
| `config/containers.json` | Singularity image paths |

Frequently adjusted keys:

| Key | Default | Meaning |
|---|---|---|
| `kraken2_decon_taxids` | `9606 10847 2` | Subtrees removed during depletion |
| `consensus_min_depth` | `10` | Below this, consensus is masked `N` |
| `min_genome_coverage` | `0.80` | Coverage below this raises a QC warning |
| `min_mapped_reads` | `1000` | Below this, a target's downstream stages are skipped |
| `pangolin_targets` | `["NC_045512", "PP115423"]` | Targets pangolin runs on |
| `freyja_targets` | `["NC_045512"]` | Targets freyja runs on — Wuhan coordinates only |

Nextclade has no list: it runs on whichever targets carry a `nextclade_dataset` path in
`genome.json`, written at build time.

---

## 4. Output

Everything a reader needs is in `final_report/`. The rest of the run directory is
intermediate — kept for debugging, not for reading.

```
final_report/
├── run_summary.tsv                  ← START HERE
├── multiqc/project_multiqc_report.html
└── {target}/
    ├── figures/                     ← LOOK HERE FIRST if you want pictures
    ├── consensus/                   ← the genomes
    ├── variants/                    ← what differs from the reference
    ├── lineage/                     ← what strain it is
    ├── qc/                          ← whether to believe the above
    ├── igv_report.{target}.html     ← interactive, no IGV needed
    └── igv_session.{target}.xml     ← for people who do have IGV
```

Outputs fall into five categories. Read them in this order.

### 4.1 Start here — the one-page answer

**`run_summary.tsv`** — one row per sample × target, and the only file most runs need.
Twenty-three columns grouped as:

| group | columns | reads as |
|---|---|---|
| Input | `input_read_pairs` | how much data went in |
| Composition | `pct_viral`, `pct_sars_cov_2`, `pct_human`, `pct_depleted` | what the library was made of, *before* depletion |
| Mapping | `mapped_reads`, `mapped_pct`, `reads_used` | how much of it hit this target |
| Assembly | `mean_depth`, `genome_length`, `consensus_masked_bases`, `pct_genome_covered` | how good the genome is |
| Variants | `n_variants` | how much it differs from the reference |
| Lineage | `pangolin_lineage`, `pangolin_qc`, `nextclade_clade`, `nextclade_qc`, `nextclade_coverage`, `freyja_lineages`, `freyja_abundances` | what strain it is, from three independent callers |
| Verdict | `qc_status` | **`pass`, or `WARN:` plus reasons** |

**`qc_status` is the column to scan first.** `WARN:LOW_GENOME_COVERAGE`,
`WARN:LOW_MAPPED_READS`, `WARN:FREYJA_FAILED` and friends say exactly which check failed,
and several reasons combine with `+`. A stage that was *skipped* rather than *failed*
never appears here — see §4.6.

Two units differ and it matters: `input_read_pairs` counts **pairs** (Kraken2 classifies a
pair as one unit) while `mapped_reads` comes from flagstat and counts **individual reads**,
so expect roughly double on paired data. `mapped_reads`/`mapped_pct` are measured *before*
unmapped records are filtered, so they describe the library; `reads_used` is what actually
reached variant calling.

### 4.2 Pictures — `{target}/figures/` and `igv_report.{target}.html`

For readers who will not open a genome browser. Four static figures per target:

| file | answers |
|---|---|
| `variant_heatmap.png` | **the whole run on one page** — variants down, samples across, cell = alt %, **grey = not called** (never 0) |
| `genome_overview.{sample}.png` | one sample across the genome: depth (log), masked regions shaded red, variant needles by allele fraction and coloured by snpEff impact |
| `coverage_comparison.png` | every sample's depth on one axis, so a weak library stands out |
| `composition.png` | what each library was made of — target virus / other viral / human / other |

`composition.png` is usually the one that explains a bad sample: low coverage is far more
often a library-composition problem than an alignment problem.

**`igv_report.{target}.html`** is a self-contained interactive report — an embedded
`igv.js` browser with the reference, variants and the read pileups around each site baked
into the file. Unlike `igv_session.{target}.xml`, it needs neither an IGV installation
nor access to the run directory, so it can simply be emailed.

### 4.3 The genomes — `{target}/consensus/`

| file | what it is |
|---|---|
| `{sample}.{target}.consensus.fa` | per-sample consensus, **depth-masked** |
| `all_samples.{target}.consensus.fa` | all samples concatenated, ready for a tree or an upload |

Positions below `consensus_min_depth` are written as `N` rather than inheriting the
reference base. An `N` therefore means *no evidence here*, not *matches the reference* —
the distinction that keeps a consensus from fabricating sequence. `consensus_masked_bases`
in the summary counts them.

### 4.4 What differs from the reference — `{target}/variants/`

Same calls, four shapes, for four different questions:

| file | answers |
|---|---|
| `aggregate.{target}.variants_matrix.filtered.tsv` | **"is this variant fixed everywhere, or a minor allele in a few samples?"** |
| `aggregate.{target}.variants_matrix.raw.tsv` | "was this variant never called, or called and then filtered out?" |
| `aggregate.{target}.variants_long.tsv` | one row per sample × variant, for scripting and plotting |
| `aggregate.{target}.snpeff.vcf.gz` (+ `.filtered.vcf.gz`) | the VCFs themselves, for any VCF-aware tool |

**The matrix is the one to open.** One row per variant, one column pair per sample:

```
POS    TYPE   GENE    AA_CHANGE  EFFECT            S1.AD    S1.PCT   S2.AD    S2.PCT
18962  SNP    ORF1ab  Q6249H     missense_variant  0,135    100.0    0,212    100.0
26212  INDEL  E                  ...inframe_del    187,52   21.8     416,99   19.2
```

`AD` is `ref,alt` read counts; `PCT` is the alt percentage. **An empty cell means the
variant was not called in that sample** — deliberately left blank rather than `0`, because
`0` would claim the site was examined and found reference. Annotation columns (`EFFECT`,
`SEVERITY`, `GENE`, `HGVS_C`, `HGVS_P`, `AA_CHANGE`) come from snpEff's highest-ranked
annotation; `AA_CHANGE` gives the short form (`Q6249H`) for substitutions and is blank for
indels and frameshifts, which defer to `HGVS_P`.

Raw versus filtered is not a formality: the SnpSift filter typically removes most calls,
and the pair is what lets you tell a genuinely absent variant from a filtered one.

### 4.5 What strain it is — `{target}/lineage/`

Present only for targets the callers can describe (see §4.7).

| file | caller | gives |
|---|---|---|
| `lineage_summary.tsv` | all three | **side-by-side table, read this one** |
| `{sample}.{target}.pangolin_lineage.csv` | pangolin | Pango lineage + QC |
| `{sample}.{target}.nextclade.tsv` | Nextclade | clade, QC, mutation list |
| `{sample}.{target}.freyja.demix` | Freyja | lineage *mixture* abundances |
| `{sample}.{target}.freyja.boot_lineages.csv` | Freyja | bootstrap confidence intervals |

Three independent methods on purpose. **Agreement is the result; disagreement is a
finding**, usually meaning a low-quality or genuinely mixed sample. Freyja is the one that
detects co-infection or multi-strain samples, which a single consensus hides by
construction.

### 4.6 Whether to believe it — `{target}/qc/` and `multiqc/`

| file | tells you |
|---|---|
| `{sample}.kraken2_decon.krona.html` | **what was in the library** — interactive, open in a browser |
| `{sample}.kraken2_decon.composition.tsv` | the same as a table; a low on-target fraction explains low coverage |
| `{sample}.{target}.mosdepth.summary.txt` | depth across the reference |
| `{sample}.{target}.lowcov_mask.bed` | exactly which positions were masked `N` |
| `{sample}.{target}.bowtie2_map.raw.flagstat` | mapping rate before filtering |
| `aggregate.{target}.mapping_summary.tsv` | all samples' mapping in one table |
| `multiqc/project_multiqc_report.html` | every QC tool, every stage, every sample, in one page |

When a sample looks wrong, the Krona chart usually explains it: low coverage is far more
often a library-composition problem than an alignment problem.

`igv_session.{target}.xml` opens all BAMs, VCFs and consensus sequences against the
reference in IGV — the fastest way to eyeball a specific variant.

### 4.7 Skipped is not failed

An absent output can mean a stage did not apply, which is not a warning:

- **No `lineage/` directory** — the target is not a lineage target. The Pango and
  Nextclade databases are SARS-CoV-2 specific.
- **No freyja files** — freyja is gated separately on Wuhan-Hu-1 coordinates
  (`freyja_targets`), because it matches a BAM against barcodes keyed to that genome and
  its answer is meaningless against any other reference.
- **A target missing downstream stages entirely** — mapping fell below `min_mapped_reads`,
  and `qc_status` says `LOW_MAPPED_READS`. One weak reference does not stop the others.

A stage that *ran and failed* always shows in `qc_status` and in the log tally at the end
of `final_report.log`. Empty directories are avoided precisely because they read as
failure.

### 4.8 Run provenance — the rest of the run directory

| path | purpose |
|---|---|
| `RUNNING` / `COMPLETED` / `FAILED` | state as a sentinel file, so a check never parses a log |
| `job_information_<ts>.tsv` | requested vs peak CPU and memory per job (`jobby`) — the basis for tuning `cluster.json` |
| `failed_jobs_<ts>.tsv` | just the failures from the above |
| `logfiles/` | per-rule logs; `master.err` is the Snakemake driver |
| `workflow/`, `config/` | the exact code and parameters that produced this run |

`workflow/` and `config/` are copied in deliberately: a run is self-describing, and
re-running later cannot silently pick up a changed repo.

---

## 5. Contribute

1. Fork the repository.
2. Create a feature branch.
3. Make your changes — give any new rule a `container:` directive, never `module load`.
4. Run `.tests/dryrun.sh` to confirm the DAG still resolves. **This is not optional:** a
   container split can silently orphan a rule, and the DAG omits it with no error.
5. Open a pull request.

---

## 6. References

<sup>1.</sup> Köster, J. & Rahmann, S. *Snakemake — a scalable bioinformatics workflow engine.* Bioinformatics (2012).
<sup>2.</sup> Chen, S. *et al.* *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics (2018).
<sup>3.</sup> Wood, D. E. *et al.* *Improved metagenomic analysis with Kraken 2.* Genome Biology (2019).
<sup>4.</sup> Ondov, B. D. *et al.* *Interactive metagenomic visualization in a web browser.* BMC Bioinformatics (2011).
<sup>5.</sup> Andrews, S. *FastQC: a quality control tool for high throughput sequence data.* (2010).
<sup>6.</sup> Langmead, B. & Salzberg, S. L. *Fast gapped-read alignment with Bowtie 2.* Nature Methods (2012).
<sup>7.</sup> Danecek, P. *et al.* *Twelve years of SAMtools and BCFtools.* GigaScience (2021).
<sup>8.</sup> Broad Institute. *Picard Toolkit.* (2019).
<sup>9.</sup> Pedersen, B. S. & Quinlan, A. R. *mosdepth: quick coverage calculation.* Bioinformatics (2018).
<sup>10.</sup> Garrison, E. & Marth, G. *Haplotype-based variant detection from short-read sequencing.* arXiv (2012).
<sup>11.</sup> Cingolani, P. *et al.* *A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff.* Fly (2012).
<sup>12.</sup> McKenna, A. *et al.* *The Genome Analysis Toolkit.* Genome Research (2010).
<sup>13.</sup> Gurevich, A. *et al.* *QUAST: quality assessment tool for genome assemblies.* Bioinformatics (2013).
<sup>14.</sup> Ewels, P. *et al.* *MultiQC: summarize analysis results for multiple tools and samples.* Bioinformatics (2016).
<sup>15.</sup> O'Toole, Á. *et al.* *Assignment of epidemiological lineages in an emerging pandemic using the pangolin tool.* Virus Evolution (2021).
<sup>16.</sup> Aksamentov, I. *et al.* *Nextclade: clade assignment, mutation calling and quality control.* JOSS (2021).
<sup>17.</sup> Karthikeyan, S. *et al.* *Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission.* Nature (2022).

[1]: https://snakemake.readthedocs.io
[2]: https://github.com/OpenGene/fastp
[3]: https://github.com/DerrickWood/kraken2
[4]: https://github.com/marbl/Krona
[5]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[6]: https://bowtie-bio.sourceforge.net/bowtie2/
[7]: https://www.htslib.org/
[8]: https://broadinstitute.github.io/picard/
[9]: https://github.com/brentp/mosdepth
[10]: https://github.com/freebayes/freebayes
[11]: https://pcingola.github.io/SnpEff/
[12]: https://gatk.broadinstitute.org/
[13]: https://quast.sourceforge.net/
[14]: https://multiqc.info/
[15]: https://github.com/cov-lineages/pangolin
[16]: https://clades.nextstrain.org/
[17]: https://github.com/andersen-lab/Freyja

---

<div align="center"><a href="#viralrecon-">Back to Top</a></div>
