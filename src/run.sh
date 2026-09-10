#!/usr/bin/env bash
# viralrecon SLURM master-job submission script.
# Called by src/run.py runner() when --mode slurm is used.
# Usage: run.sh slurm -j <jobname> -b <bindpaths> -o <outdir> -c <cache> -t <tmpdir>
set -euo pipefail

MODE="$1"; shift

JOBNAME="viralrecon"
# Paths that must be visible inside every container: the run itself, the
# reference/image tree, and node-local scratch.
BINDPATHS="/data/RTB_GRS,/fdb"
OUTDIR=""
CACHEDIR=""
TMPDIR=""            # defaults to ${OUTDIR}/tmp once OUTDIR is known

while getopts "j:b:o:c:t:" opt; do
    case "$opt" in
        j) JOBNAME="$OPTARG" ;;
        b) BINDPATHS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        c) CACHEDIR="$OPTARG" ;;
        t) TMPDIR="$OPTARG" ;;
        *) echo "Unknown option: $opt" >&2; exit 1 ;;
    esac
done

LOGDIR="${OUTDIR}/logfiles"
SLURM_DIR="${LOGDIR}/slurmfiles"
# Temporary files live inside the run directory, not in node-local /lscratch.
# Keeping them here means they are on the same filesystem as their outputs, they
# survive a job for inspection when something fails, and the run directory is
# self-contained: nothing it produced is stranded on a compute node.
: "${TMPDIR:=${OUTDIR}/tmp}"
mkdir -p "${LOGDIR}" "${SLURM_DIR}" "${TMPDIR}"

# Clear the previous run's status markers at SUBMIT time, not in Snakemake's
# onstart hook. onstart only fires once the master job starts, so between
# sbatch and that moment a stale COMPLETED sits over a job that is queued or
# running - which reads as a finished, successful run to anything checking the
# sentinel. Nothing downstream can tell the difference.
rm -f "${OUTDIR}/COMPLETED" "${OUTDIR}/FAILED" "${OUTDIR}/RUNNING"

# Write a self-contained kickoff.sh so the run can be re-submitted
# directly (sbatch kickoff.sh) without going through the Python CLI.
cat > "${OUTDIR}/kickoff.sh" << EOF
#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=6g
#SBATCH --time=0-12:00:00
#SBATCH --job-name="${JOBNAME}"
#SBATCH --output="${LOGDIR}/master.log"
#SBATCH --error="${LOGDIR}/master.log"
#SBATCH --open-mode=append
set -euo pipefail
module load singularity/4.3.7 || module load singularity || true
export SINGULARITY_CACHEDIR="${CACHEDIR}"
export TMPDIR="${TMPDIR}"
export SINGULARITY_TMPDIR="${TMPDIR}"
export APPTAINER_TMPDIR="${TMPDIR}"
mkdir -p "${TMPDIR}"
cd "${OUTDIR}"

# Rules declare their own resources via allocated() in the Snakefile;
# cluster.json is loaded inside the Snakefile, not via --cluster-config.
#
# --keep-going matters with several targets: the per-target branches are
# independent, so a reference that a library barely maps to should not abandon
# work already queued for the others. Snakemake still exits non-zero and the
# FAILED sentinel is still written, so a partial run is never mistaken for a
# clean one.
CLUSTER_OPTS="sbatch --cpus-per-task {threads} --mem {resources.mem} --time {resources.time} --partition {resources.partition} --job-name viralrecon.{rule} --output ${SLURM_DIR}/slurm-%j_{rule}.out --error ${SLURM_DIR}/slurm-%j_{rule}.out"

snakemake -pr \\
    --rerun-incomplete \\
    --keep-going \\
    --use-singularity \\
    --singularity-args "--bind ${BINDPATHS}" \\
    --cluster "\$CLUSTER_OPTS" \\
    --jobs=20 \\
    --max-jobs-per-second=1 \\
    --max-status-checks-per-second=0.01 \\
    --latency-wait=120 \\
    --keep-incomplete \\
    --rerun-triggers mtime \\
    --configfile=config.json \\
    -s workflow/Snakefile
EOF

chmod +x "${OUTDIR}/kickoff.sh"

JOBID=$(sbatch --parsable "${OUTDIR}/kickoff.sh")
echo "${JOBID}" > "${LOGDIR}/mjobid.log"
echo "${JOBID}"
