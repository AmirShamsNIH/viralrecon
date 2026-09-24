#!/usr/bin/env bash
# SLURM master-job submission, called by src/run.py for --mode slurm.
# Usage: run.sh slurm -j <jobname> -b <bindpaths> -o <outdir> -c <cache> -t <tmpdir>
set -euo pipefail

MODE="$1"
shift

JOBNAME="viralrecon"
# Fallback binds when -b is not given: the reference/image tree and /fdb.
BINDPATHS="/data/RTB_GRS,/fdb"
OUTDIR=""
CACHEDIR=""
TMPDIR=""  # defaults to ${OUTDIR}/tmp once OUTDIR is known

while getopts "j:b:o:c:t:" opt; do
    case "$opt" in
        j) JOBNAME="$OPTARG" ;;
        b) BINDPATHS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        c) CACHEDIR="$OPTARG" ;;
        t) TMPDIR="$OPTARG" ;;
        *)
            echo "Unknown option: $opt" >&2
            exit 1
            ;;
    esac
done

LOGDIR="${OUTDIR}/logfiles"
SLURM_DIR="${LOGDIR}/slurmfiles"
# Temporary files stay in the run directory beside their outputs, so nothing is
# stranded on a compute node and a failure can be inspected.
: "${TMPDIR:=${OUTDIR}/tmp}"
mkdir -p "${LOGDIR}" "${SLURM_DIR}" "${TMPDIR}"

# Clear old status markers at submit time: onstart fires only once the master job
# starts, so a stale COMPLETED would otherwise sit over a queued run.
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

# Rules take resources from allocated() in the Snakefile, not --cluster-config.
# --keep-going lets independent targets finish; any failure still writes FAILED.
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
