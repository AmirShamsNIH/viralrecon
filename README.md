<div align="center">

<h1>viralrecon</h1>
<b>Viral Variation Analysis Pipeline</b><br/>
<i>NIH Biowulf · GRS BigSky · NIAID Skyline · <a href="https://github.com/OpenOmics/baseline">OpenOmics/baseline</a> structure</i>

<br/><br/>

<a href="#"><img alt="Version" src="https://img.shields.io/badge/version-0.1.0-blue"></a>
<a href="#22-dependencies"><img alt="Snakemake" src="https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen"></a>
<a href="#23-containers"><img alt="Singularity" src="https://img.shields.io/badge/singularity-only-%23663399"></a>
<a href="#"><img alt="Platform" src="https://img.shields.io/badge/platform-Biowulf%20%7C%20BigSky%20%7C%20Skyline-orange"></a>

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

viralrecon takes raw Illumina FASTQ files from clinical or environmental samples and a
set of reference accessions. It reports what is in each sample, builds per-sample
consensus genomes, calls and annotates variants (including sub-consensus, intra-host
variants), and assigns clades. Results for the reader are collected in `final_report/`.

Two design choices:

- **Every external tool comes from a pinned Singularity image.** No `module load`
  appears in the workflow or in the reference-building CLI. A cluster module can be
  upgraded or removed underneath a pipeline and its version is not recorded in the run
  directory; a pinned `.sif` path is. For tools that bundle a dataset, like Nextclade, a
  floating version changes clade assignments.
- **The pipeline is virus-agnostic.** No rule assumes a particular virus. Everything that
  differs between viruses (the FASTA, the annotation, the taxid, the Nextclade dataset)
  is recorded on the reference, so adding a virus needs a `viralrecon build`, not a code
  change. Targets are selected by accession, never by name. A run may carry several
  unrelated viruses at once; every post-alignment file is namespaced `{sample}.{target}.*`
  so targets cannot mix.

Orchestration is [Snakemake][1], which handles the job DAG, SLURM submission, restarts,
and container invocation.

---

## 2. Overview of the pipeline

Every sample is processed against every selected reference. Steps run in this order,
scheduled as SLURM jobs by Snakemake.

<p align="center"><img src="docs/workflow.svg" width="100%" alt="viralrecon workflow: 1 Reference build, 2 Read clean-up, 3 Composition, 4 Alignment, 5 Variant calling, 6 Consensus, 7 Clade, 8 Reporting"></p>

The dashed step runs only for references that carry a Nextclade dataset.

### 2.1 Stages

The rules are grouped into these stages; `lineage` is optional.

| # | Stage | Purpose | Key tools |
|---|-------|---------|-----------|
| 1 | **build_environment** | Per-run reference setup | samtools |
| 2 | **pre_process** | Repair, trim, profile, deplete | BBTools · [fastp][2] · [Kraken2][3] · [Krona][4] · [FastQC][5] |
| 3 | **alignment** | Map to each target, QC the mapping | [Bowtie2][6] · [samtools][7] · [Picard][8] · [mosdepth][9] |
| 4 | **variant_calling** | Call, normalise, annotate, build consensus | [FreeBayes][10] · [bcftools][7] · [SnpEff/SnpSift][11] · [GATK4][12] |
| 5 | **report** | Aggregate across samples, assemble `final_report/` | bcftools · GATK4 · [QUAST][13] · [MultiQC][14] |
| - | **lineage** *(opt.)* | Clade assignment | [Nextclade][15] |

#### Quality control and decontamination

- **BBTools `reformat`** repairs the raw pairs first: junk bases, IUPAC codes, broken
  reads, and out-of-range quality scores are fixed or dropped before anything reads them.
  A truncated `.gz` fails here, which points to the input transfer rather than the data.
- **fastp** trims adapters (auto-detected per pair), low-complexity reads, and low-quality
  tails.
- **Kraken2** classifies the trimmed reads against the standard database before
  depletion, so contamination is reported rather than discarded unseen. **Krona** renders
  that profile as an interactive HTML chart.
- **Depletion removes taxon subtrees.** Reads assigned anywhere in the subtree of `9606`
  (human), `10847` (phiX), or `2` (bacteria) are removed. Everything else is kept,
  including unclassified reads, which is where a novel or divergent virus appears. Taxids
  are configurable.
- **FastQC** runs on the depleted reads, and a single project-wide **MultiQC** report
  aggregates every stage across every sample into `final_report/multiqc/`.

#### Alignment and variant calling

- **Bowtie2** (`--local`) maps the depleted reads to each target independently. A target
  whose mapping falls below `min_mapped_reads` raises a QC warning and its downstream
  stages are skipped; the run continues with the other targets.
- **samtools** filters unmapped, secondary, and supplementary alignments; **Picard**
  collects insert-size and alignment metrics; **mosdepth** produces per-base and windowed
  depth.
- **FreeBayes** runs in `--pooled-continuous` mode. A viral sample is a population, not a
  diploid individual, so variants are called on allele fraction instead of genotype. The
  default floor of `--min-alternate-fraction 0.02` reports intra-host minor variants that
  a genotype-based caller would discard.
- **bcftools norm** splits multi-allelic records; **SnpEff** and **SnpSift** annotate;
  **GATK4** flattens calls into a table, which is emitted in two shapes: a long table
  (one row per sample × variant, for scripting) and a variant × sample matrix: one
  row per variant with an allele-depth and percentage column per sample, plus parsed
  effect, severity, gene and amino-acid change. The matrix shows whether a variant is
  fixed across the run or a minor allele in a few samples; an empty cell means not called
  in that sample, and is left empty rather than zeroed.
- **Raw and filtered sets** both reach the aggregate, as `variants_matrix.raw.tsv` and
  `variants_matrix.filtered.tsv`. The filtered set applies the SnpSift expression and is
  the one to report. The raw set is what FreeBayes called under its own thresholds;
  comparing the two separates a variant that was never called from one that was called
  and then filtered out.
- **Consensus is depth-masked.** Positions below `consensus_min_depth` are written as `N`
  instead of the reference base, so the consensus contains no reference sequence where
  there is no evidence.

#### Lineage

Clade assignment is done by **[Nextclade][15]**, and it works for any virus that has a
Nextclade dataset: SARS-CoV-2, influenza, RSV, mpox, Ebola, Marburg, dengue, measles and
about a hundred more. `resources/nextclade_datasets.tsv` lists them all.

- **The dataset is matched automatically at build time.** `viralrecon build` runs
  `nextclade sort` on the reference FASTA, which compares it against every published
  dataset. When all records match one dataset, it lands in `<reference>/nextclade/` and
  its path is recorded in `genome.json`. When nothing matches, lineage is skipped.
  `--no-nextclade` opts out; there is no other switch.
- **The dataset belongs to the reference, not the run.** Nextclade runs on every target
  that carries one.
- **A target without a dataset is skipped.** Using another virus's dataset would return
  a wrong clade instead of failing.
- **Nextclade aligns to the dataset's own reference.** The clade does not depend on which
  reference the reads were mapped to, so two targets for the same virus should agree.
- **Pick the dataset that matches what was sequenced.** QC is scored against the
  dataset's reference tree. In validation, six JN.1-lineage samples (clade 24A) scored
  `mediocre`/`bad` against the Wuhan-rooted `sars-cov-2` dataset and `good` against the
  BA.2.86-rooted one, with the same clade call from both. Ancestral samples go the other
  way: clade 19B reads as `outgroup` against a BA.2.86-rooted tree.

### 2.2 Dependencies

| Requirement | Version |
|---|---|
| Snakemake | ≥ 7.0 |
| Singularity / Apptainer | ≥ 3.5 |
| Python | ≥ 3.8 |

Input is paired or single-end Illumina FASTQ; the layout is detected from whether any
R2 file is present, not configured. Single-end has been verified against the same
samples run paired: identical clade calls and fixed markers, and sub-consensus
frequencies within about a point, at proportionally lower depth.

Nothing else is needed locally; every tool is containerised.

### 2.3 Containers

Images are resolved from `config/containers.json`, which names two **roots** per
platform and writes every image against one of them:

| root | what it holds |
|---|---|
| `{shared}` | a read-only library maintained by someone else |
| `{ours}` | images we pull or build for this pipeline |

Only the roots change between platforms; the image file names, which are the version
pins of record, are written once and cannot drift apart per platform. On Biowulf the
roots are `/data/OpenOmics/SIFs` and `/data/RTB_GRS/references/singularity`. A root may
contain `{repo_parent}`, which expands to the directory holding the clone, so a cluster
with no institutional reference tree keeps its images beside the checkout.

`src/containers.py` resolves them for the platform given to
`viralrecon build/run --platform`. The same idea covers databases the pipeline does not
build: `config/config.json` holds a `paths` block keyed by platform for the Kraken2
index and the Krona taxonomy, and `config/cluster.json` holds `__partition__` for the
SLURM queue name. All three are resolved before Snakemake starts, so no rule ever asks
which cluster it is on.

A Snakemake rule may declare only **one** container, so a step chaining tools from
different images is split into one rule per tool. When adding a rule, give it a
`container:` directive, never a module.

### 2.4 Installation

```bash
git clone https://github.com/AmirShamsNIH/viralrecon.git
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
| SLURM partition | `norm`, the only partition `cluster.json` uses |
| Shared references | read access to `/data/OpenOmics/SIFs` |

**Software on the submitting shell**

```bash
module load python/3.10        # provides snakemake 7.30.1
```

Snakemake is the one tool not run from a container, because it launches the
containers. The master job inherits the submitting environment, so if `snakemake` is not on your `PATH` when you submit, the master
job fails immediately. Everything else is loaded by the job itself: the master
script runs `module load singularity` and no rule loads anything at all.

**Data that must exist on disk**

| path | what it is |
|---|---|
| `/fdb/kraken/20260226_standard_kraken2` | Kraken2 standard database, ~98 GB. Biowulf shared, nothing to install |
| `/data/RTB_GRS/references/krona/taxonomy` | Krona taxonomy, staged once. The Krona image ships a placeholder that produces empty charts |
| `/data/OpenOmics/SIFs` | shared lab image library |
| `/data/RTB_GRS/references/singularity` | images built or pulled for this pipeline |

**Which node does what.** Singularity and internet access are on different hosts:

| task | where | why |
|---|---|---|
| `viralrecon build` | **compute node** (`sbatch`) | needs `singularity`, which the login node does not have |
| `viralrecon run` | either; it only submits | the master job it creates runs on a compute node |
| `git push` | **login node** | compute nodes cannot resolve external hostnames |

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

$VIRALRECON build --virus SARS --accession NC_045512.2 --output "$REF"

$VIRALRECON run --input /path/to/reads/*_R[12]_001.fastq.gz \
    --output "$OUT" --genome "$REF/genome.json" \
    --targets SARS_NC_045512.2
```

Submit it with `sbatch`, having run `module load python/3.10` first so the
`run` step can find Snakemake.

**Scratch.** The pipeline keeps its temporary files in `$OUTDIR/tmp` rather than
`lscratch`, so no `--gres=lscratch` allocation is needed.

### 2.6 Platform profile: BigSky

Pass `--platform BIGSKY` to both `build` and `run`. Everything that differs from
Biowulf is config, not code: image roots in `config/containers.json`, database
paths in `config/config.json`, and the queue name in `config/cluster.json`.

**The layout.** BigSky has no institutional reference tree, so the deployment is
one directory with the references beside the checkout:

```
/data/rml_ngs/viralrecon/
├── viralrecon/     the clone
├── singularity/    the images  ({repo_parent}/singularity)
└── references/     `viralrecon build --output` goes here
```

**What differs**

| | Biowulf | BigSky |
|---|---|---|
| SLURM partition | `norm` | `all` (also `himem`, 4 TB × 2, and `gpu`) |
| `singularity` | compute nodes only | submit and compute nodes |
| Internet | login node only; proxy for HTTP | everywhere, including compute nodes; no proxy needed |
| Snakemake | `module load python/3.10` | no module (see below) |
| Node-local scratch | `/lscratch/$SLURM_JOB_ID` | none; `/tmp` is mounted noexec |
| Kraken2 database | `/fdb/kraken/20260226_standard_kraken2` | `/data/rml_ngs/kraken_db/K2/k2_standard` |
| Krona taxonomy | `/data/RTB_GRS/references/krona/taxonomy` | `/data/rml_ngs/ngs_dbs/krona/taxonomy` |
| `git` | on `PATH` | `module load git` |

Because compute nodes reach the internet directly, `viralrecon build` needs no
proxy exports on BigSky, and there is no login/compute split: one node can build
images, build references, and push to GitHub.

**Snakemake.** There is no snakemake module. Create it once:

```bash
PY=$(module load python/3.11.9-4z43o4e >/dev/null 2>&1; command -v python3)
$PY -m venv /data/rml_ngs/viralrecon/sm_venv
/data/rml_ngs/viralrecon/sm_venv/bin/pip install "snakemake==7.30.1" "pulp<2.8"
```

The `pulp<2.8` pin is required. Snakemake 7 calls `pulp.list_solvers`, which
pulp 3 renamed, and an unpinned install fails on `--version` with
`AttributeError: module 'pulp' has no attribute 'list_solvers'` before it reads a
single rule.

**Images.** The pinned set is ~6 GB. `/data/openomics/SIFs` exists on BigSky but
carries only 1 of the 20 images, so both roots point at `singularity/` beside the
clone and everything is staged there. To pull an image directly, set two temporary
directories first:

```bash
export APPTAINER_TMPDIR=/data/rml_ngs/viralrecon/tmp
export PROOT_TMP_DIR=/data/rml_ngs/viralrecon/tmp
apptainer pull image.sif docker://quay.io/biocontainers/<tool>:<tag>
```

Without them the pull fails inside proot with
`mksquashfs: No such file or directory`; the real cause is that `/tmp` is
mounted `noexec`, and apptainer falls back to proot because there is no setuid
starter and no `/etc/subuid` entry (so `--fakeroot` is unavailable, exactly as on
Biowulf).

**Apptainer version string.** Snakemake 7 parses `singularity --version` with
`packaging.Version`, and BigSky's apptainer answers `apptainer version
1.5.0-1.el9`. The RPM release suffix is not PEP 440, so the master job dies in
28 seconds with `InvalidVersion: '1.5.0-1.el9'` before running a single rule.
A wrapper script earlier on `PATH` fixes it:

```bash
cat > /data/rml_ngs/viralrecon/bin/singularity <<'SH'
#!/usr/bin/env bash
if [ "$1" = "--version" ]; then echo "apptainer version 1.5.0"; exit 0; fi
exec /usr/bin/singularity "$@"
SH
chmod +x /data/rml_ngs/viralrecon/bin/singularity
export PATH=/data/rml_ngs/viralrecon/bin:/data/rml_ngs/viralrecon/sm_venv/bin:$PATH
```

**Scratch.** `$OUTDIR/tmp` on GPFS is used for everything, as on Biowulf, so no
`--gres=lscratch` is needed. BigSky has no `/lscratch` in any case.

### 2.7 Platform profile: Skyline

Pass `--platform SKYLINE` to both `build` and `run`. Skyline (NIAID, submit host
`ai-hpcsubmit1.niaid.nih.gov`) has its own `/data`, separate from BigSky's, so none of
the BigSky paths exist there.

**The layout.** The shared install follows the other pipelines in `/data/openomics/prod`:
one directory per release, a `latest` link, and the images and Krona taxonomy beside
them, shared by every release:

```
/data/openomics/prod/viralrecon/
├── v0.1.0/                 the clone
├── latest -> v0.1.0
└── viralrecon_resources/
    ├── singularity/         the images  ({repo_parent}/viralrecon_resources/singularity)
    └── krona/taxonomy/      Krona taxonomy, copied from Biowulf
```

**What differs**

| | Biowulf | Skyline |
|---|---|---|
| SLURM partition | `norm` | `all` (also `gpu`) |
| `singularity` | compute nodes only | apptainer 1.5.0 on submit and compute nodes |
| Internet | login node only; proxy for HTTP | everywhere; no proxy needed |
| Snakemake | `module load python/3.10` | 7.22.0 in a user miniconda, or `module load snakemake/7.22.0-ufanewz` |
| Node-local scratch | `/lscratch/$SLURM_JOB_ID` | none; `/tmp` is mounted noexec |
| Kraken2 database | `/fdb/kraken/20260226_standard_kraken2` | `/data/bio_db/kraken_db/plus_PFV_Oct2025` |
| Krona taxonomy | `/data/RTB_GRS/references/krona/taxonomy` | `viralrecon_resources/krona/taxonomy` |
| Kraken2 loading | `--memory-mapping` | whole database into RAM |

**Kraken2.** Skyline has no standard database. `plus_PFV_Oct2025` is the closest: the
standard libraries (archaea, bacteria, viral, plasmid, human, UniVec_Core) plus protozoa
and fungi. Depletion removes the same taxa, but composition percentages will not match
a Biowulf run exactly. Kraken2 runs without `--memory-mapping` here: paging the 101 GB
database in at random from GPFS left a 20 MB sample unfinished after an hour, while one
sequential load fits the rule's 150 GB request.

**Snakemake.** Stay on 7.x; the `snakemake/8.18.2` modules will not run this
workflow. Snakemake 7.22 parses the apptainer version with `LooseVersion`, so the
BigSky version shim is not needed.

**Images.** The pinned set (~6 GB) was copied from Biowulf. `/data/openomics/SIFs`
exists on Skyline but holds almost none of the pins, so both roots point at
`viralrecon_resources/singularity`.

**Scratch.** `$OUTDIR/tmp` is used for everything, as on the other clusters.

## 3. Run the pipeline

### 3.0 Start from a template

`execution/` holds one batch script per cluster (`viralrecon_skyline.sh`, `viralrecon_biowulf.sh`) that builds
the references and runs the pipeline. Copy the one for your cluster into your project,
set `WORKDIR`, `FASTQ_DIR` and one `build` line per reference, and submit it from there:

```bash
cp /data/openomics/prod/viralrecon/latest/execution/viralrecon_skyline.sh /data/<group>/<project>/
cd /data/<group>/<project> && sbatch viralrecon_skyline.sh
```

Without `--targets`, `run` uses every reference in `genome.json`. The job log lands next
to the script as `viralrecon.<jobid>.log`; the run itself reports under `$OUT` (§4).

### 3.1 Build a reference

Run once per accession. References accumulate in a single `genome.json`: a new accession
is added to an existing build, and an accession already present is verified
file-by-file and skipped if complete, or completed if not.

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs
```

A Nextclade dataset is matched to the reference automatically (§2.1, *Lineage*). A
virus with no dataset, or a segmented reference whose segments match different datasets,
is skipped with a message. To build a reference without one:

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs --no-nextclade
```

On an existing reference, `--no-nextclade` removes the dataset, and a plain rebuild of a
reference without one matches it then. On Biowulf this step needs the proxy set (§2.5);
BigSky and Skyline reach the internet directly. `resources/nextclade_datasets.tsv` lists
the catalogue; refresh it with `nextclade dataset list --json`.

`--notes` records curated knowledge about a reference that cannot be derived from its
files. The note is stored in `genome.json`, echoed whenever a run selects that
target, and carried forward across rebuilds unless given again:

```bash
./viralrecon build --virus SARS --accession NC_045512.2 --output /data/refs \
                   --notes "contains long N runs; do not use for consensus"
```

To provide your own sequence instead of downloading it from NCBI (the two are mutually
exclusive):

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
directory; `config/` is copied once and does not otherwise refresh:

```bash
./viralrecon run -i ... -o ... --genome ... --overwrite-pipeline-template
```

Pass `--platform BIGSKY` or `--platform SKYLINE` to both `build` and `run` on those
clusters; Biowulf is the default.
`unlock` releases a stale Snakemake lock on an output directory.

### 3.3 Configuration

| File | Contents |
|---|---|
| `config/config.json` | Pipeline options and every tool parameter |
| `config/cluster.json` | Per-rule SLURM resources |
| `config/genome.json` | Template only; the real registry is `<build --output>/genome.json` |
| `config/containers.json` | Singularity image paths |

Frequently adjusted keys:

| Key | Default | Meaning |
|---|---|---|
| `kraken2_decon_taxids` | `9606 10847 2` | Subtrees removed during depletion |
| `consensus_min_depth` | `10` | Below this, consensus is masked `N` |
| `min_genome_coverage` | `0.80` | Coverage below this raises a QC warning |
| `min_mapped_reads` | `1000` | Below this, a target's downstream stages are skipped |

Nextclade has no list: it runs on whichever targets carry a `nextclade_dataset` path in
`genome.json`, matched at build time.

---

## 4. Output

Everything a reader needs is in `final_report/`. The rest of the run directory is
intermediate, kept for debugging, not for reading.

```
final_report/
├── run_summary.tsv                  one row per sample × target; read first
├── multiqc/project_multiqc_report.html
└── {target}/
    ├── figures/                     static plots
    ├── consensus/                   consensus genomes
    ├── variants/                    variant tables and VCFs
    ├── lineage/                     Nextclade clades
    ├── qc/                          composition, depth and mapping QC
    ├── igv_report.{target}.html     interactive, no IGV needed
    └── igv_session.{target}.xml     session file for IGV
```

The sections below follow the order in which to read the outputs.

### 4.1 Run summary: `run_summary.tsv`

One row per sample × target; for most runs this is the only file needed.
Nineteen columns grouped as:

| group | columns | reads as |
|---|---|---|
| Input | `input_read_pairs` | how much data went in |
| Composition | `pct_viral`, `pct_target`, `pct_human`, `pct_depleted` | what the library was made of, *before* depletion |
| Mapping | `mapped_reads`, `mapped_pct`, `reads_used` | how much of it hit this target |
| Assembly | `mean_depth`, `genome_length`, `consensus_masked_bases`, `pct_genome_covered` | how good the genome is |
| Variants | `n_variants` | how much it differs from the reference |
| Lineage | `nextclade_clade`, `nextclade_qc`, `nextclade_coverage` | what clade it is |
| Verdict | `qc_status` | `pass`, or `WARN:` plus reasons |

Check `qc_status` first. `WARN:LOW_GENOME_COVERAGE`, `WARN:LOW_MAPPED_READS`,
`WARN:NEXTCLADE_FAILED` and the other codes name the check that failed, and several
reasons combine with `+`. A stage that was skipped rather than failed never appears here
(see §4.7).

The read counts use different units: `input_read_pairs` counts pairs (Kraken2 classifies a
pair as one unit) while `mapped_reads` comes from flagstat and counts individual reads,
so expect roughly double on paired data. `mapped_reads`/`mapped_pct` are measured *before*
unmapped records are filtered, so they describe the library; `reads_used` is what actually
reached variant calling.

### 4.2 Figures: `{target}/figures/` and `igv_report.{target}.html`

Four static figures per target, for readers who do not use a genome browser:

| file | answers |
|---|---|
| `variant_heatmap.png` | the whole run on one page: variants down, samples across, cell = alt %, grey = not called (never 0) |
| `genome_overview.{sample}.png` | one sample across the genome: depth (log), masked regions shaded red, variant needles by allele fraction and coloured by snpEff impact |
| `coverage_comparison.png` | every sample's depth on one axis, so a weak library stands out |
| `composition.png` | what each library was made of: target virus / other viral / human / other |

`composition.png` is usually the one that explains a bad sample: low coverage is far more
often a library-composition problem than an alignment problem.

**`igv_report.{target}.html`** is a self-contained interactive report, an embedded
`igv.js` browser with the reference, variants and the read pileups around each site baked
into the file. Unlike `igv_session.{target}.xml`, it needs neither an IGV installation
nor access to the run directory, so it can be emailed.

### 4.3 Consensus genomes: `{target}/consensus/`

| file | what it is |
|---|---|
| `{sample}.{target}.consensus.fa` | per-sample consensus, depth-masked |
| `all_samples.{target}.consensus.fa` | all samples concatenated, ready for a tree or an upload |

Positions below `consensus_min_depth` are written as `N` rather than inheriting the
reference base, so an `N` means no evidence at that position, not a match to the
reference. `consensus_masked_bases` in the summary counts them.

### 4.4 Variants: `{target}/variants/`

The same calls in four forms:

| file | answers |
|---|---|
| `aggregate.{target}.variants_matrix.filtered.tsv` | is this variant fixed everywhere, or a minor allele in a few samples? |
| `aggregate.{target}.variants_matrix.raw.tsv` | was this variant never called, or called and then filtered out? |
| `aggregate.{target}.variants_long.tsv` | one row per sample × variant, for scripting and plotting |
| `aggregate.{target}.snpeff.vcf.gz` (+ `.filtered.vcf.gz`) | the VCFs themselves, for any VCF-aware tool |

The filtered matrix is the main table: one row per variant, one column pair per sample:

```
POS    TYPE   GENE    AA_CHANGE  EFFECT            S1.AD    S1.PCT   S2.AD    S2.PCT
18962  SNP    ORF1ab  Q6249H     missense_variant  0,135    100.0    0,212    100.0
26212  INDEL  E                  ...inframe_del    187,52   21.8     416,99   19.2
```

`AD` is `ref,alt` read counts; `PCT` is the alt percentage. An empty cell means the
variant was not called in that sample. It is left blank rather than `0`, because `0`
would claim the site was examined and found reference. Annotation columns (`EFFECT`,
`SEVERITY`, `GENE`, `HGVS_C`, `HGVS_P`, `AA_CHANGE`) come from snpEff's highest-ranked
annotation; `AA_CHANGE` gives the short form (`Q6249H`) for substitutions and is blank for
indels and frameshifts, which defer to `HGVS_P`.

The SnpSift filter typically removes most calls, so compare the raw and filtered
matrices to tell an absent variant from a filtered one.

### 4.5 Clades: `{target}/lineage/`

Present only for targets whose reference carries a Nextclade dataset (see §4.7).

| file | gives |
|---|---|
| `lineage_summary.tsv` | clade, QC and coverage per sample; the main table |
| `{sample}.{target}.nextclade.tsv` | full Nextclade output: clade, QC breakdown, mutation list |

Read `nextclade_qc` beside the clade. On a well-covered sample, `mediocre` or `bad` most
often means the dataset does not match what was sequenced rather than a problem with the
sample; check the clade against a dataset rooted closer to it before distrusting the
data (see *Lineage* in §2.1). `nextclade_coverage` is the fraction of the dataset's
reference the consensus covers.

### 4.6 QC: `{target}/qc/` and `multiqc/`

| file | tells you |
|---|---|
| `{sample}.kraken2_decon.krona.html` | what was in the library, as an interactive chart for a browser |
| `{sample}.kraken2_decon.composition.tsv` | the same as a table; a low on-target fraction explains low coverage |
| `{sample}.{target}.mosdepth.summary.txt` | depth across the reference |
| `{sample}.{target}.lowcov_mask.bed` | exactly which positions were masked `N` |
| `{sample}.{target}.bowtie2_map.raw.flagstat` | mapping rate before filtering |
| `aggregate.{target}.mapping_summary.tsv` | all samples' mapping in one table |
| `multiqc/project_multiqc_report.html` | every QC tool, every stage, every sample, in one page |

When a sample looks wrong, start with the Krona chart or `composition.png` (§4.2).

`igv_session.{target}.xml` opens all BAMs, VCFs and consensus sequences against the
reference in IGV, for inspecting a specific variant.

### 4.7 Skipped stages

An absent output can mean a stage did not apply, which is not a warning:

- **No `lineage/` directory**: no Nextclade dataset matched the reference at build time,
  or it was built with `--no-nextclade`.
- **A target missing downstream stages entirely**: mapping fell below `min_mapped_reads`,
  and `qc_status` says `LOW_MAPPED_READS`. One weak reference does not stop the others.

A stage that *ran and failed* always shows in `qc_status` and in the log tally at the end
of `final_report.log`. A skipped stage leaves no empty directory, since one would read
as a failure.

### 4.8 Run provenance: the rest of the run directory

| path | purpose |
|---|---|
| `RUNNING` / `COMPLETED` / `FAILED` | state as a sentinel file, so a check never parses a log |
| `job_information_<ts>.tsv` | requested vs peak CPU and memory per job (`jobby`), the basis for tuning `cluster.json` |
| `failed_jobs_<ts>.tsv` | just the failures from the above |
| `logfiles/` | per-rule logs; `master.log` is the full Snakemake driver record, `master.err` holds failures only and stays empty on success |
| `workflow/`, `config/` | the exact code and parameters that produced this run |

`workflow/` and `config/` are copied in so a run is self-describing, but they refresh
differently. `workflow/` is re-copied from the repo on every launch, so code fixes reach an
existing run directory; a changed helper script reruns the rules that use it. `config/` is
copied once so per-run edits survive, which means a new config key does not reach an
existing run directory until you pass `--overwrite-pipeline-template`.

---

## 5. Contribute

1. Fork the repository.
2. Create a feature branch.
3. Make your changes: give any new rule a `container:` directive, never `module load`.
4. Run `viralrecon run ... --dry-run` on test data to confirm the DAG still resolves.
   This is required: a container split can orphan a rule, and the DAG omits it
   without an error.
5. Open a pull request.

---

## 6. References

<sup>1.</sup> Köster, J. & Rahmann, S. *Snakemake: a scalable bioinformatics workflow engine.* Bioinformatics (2012).
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
<sup>15.</sup> Aksamentov, I. *et al.* *Nextclade: clade assignment, mutation calling and quality control.* JOSS (2021).

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
[15]: https://clades.nextstrain.org/

---

<div align="center"><a href="#viralrecon">Back to Top</a></div>
