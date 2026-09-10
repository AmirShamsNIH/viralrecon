# ############################################################################
# hooks.smk — pipeline lifecycle hooks
#
# Sentinel files
# --------------
# Snakemake's exit status is only visible in the master job's log, which is
# awkward to check from a script, a cron job, or another pipeline. These hooks
# keep exactly one marker file in the output directory root at all times:
#
#   RUNNING    the master job is alive
#   COMPLETED  every rule finished, exit 0
#   FAILED     at least one rule failed
#
# So `[ -f COMPLETED ]` is a complete answer to "did this run work", with no
# log parsing. The previous marker is always removed first, so a rerun after a
# failure cannot leave both FAILED and COMPLETED behind.
#
# final_report/.done is a different thing and both are worth having: .done says
# the report was assembled, COMPLETED says the whole DAG succeeded.
#
# Job accounting
# --------------
# On exit, `jobby` (workflow/scripts/jobby, from OpenOmics/RNA-seek) is run
# over every SLURM job id the master submitted, producing:
#
#   job_information_<timestamp>.tsv   every job: state, requested cpus/mem and
#                                     time, plus the cpu_max / mem_max actually
#                                     used, node, and stdout/stderr paths
#   failed_jobs_<timestamp>.tsv       the FAILED subset, for triage
#
# mem_max and cpu_max are the point of this: they turn cluster.json from a set
# of guesses into something measurable. It runs on both success and failure -
# a failed run is exactly when the resource numbers matter most.
# ############################################################################

# Both master streams go to master.log, opened in append mode so a re-run does
# not destroy the previous run's driver log. Snakemake writes nearly everything
# to stderr, including ordinary progress, so master.log is the full record and
# master.err is reserved for failures only: empty when a run succeeds, and the
# first thing to read when the FAILED sentinel appears.
#
# Because master.log now spans runs, anything derived from it has to look only
# at the current run. onstart writes a marker line and _CURRENT_SLICE keeps
# just the text after the last one - otherwise jobby would report SLURM ids
# from previous runs as if they belonged to this one.
_MASTER_LOG = join("logfiles", "master.log")
_MASTER_ERR = join("logfiles", "master.err")
_RUN_MARKER = "=== viralrecon run started"

_CURRENT_SLICE = r"""
CURRENT=$(mktemp)
awk -v m="{marker}" 'index($0, m)==1 {{buf=""}} {{buf = buf $0 ORS}} END {{printf "%s", buf}}'     "{log}" > "$CURRENT" 2>/dev/null || : > "$CURRENT"
""".replace("{marker}", _RUN_MARKER)

_JOBBY = r"""
sleep 15
rm -f COMPLETED FAILED RUNNING
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# Collect the SLURM ids this run submitted. Nothing is submitted when the DAG
# is already satisfied, so an empty list is normal, not an error.
# The `|| true` is load-bearing, not defensive noise: grep exits 1 when it
# matches nothing, and under `set -e` a VAR=$(...) assignment takes the
# command's status, so a run whose DAG was already satisfied - submitting no
# jobs at all - aborted this hook and reported FAILED despite succeeding.
JOBIDS=$(grep --color=never "^Submitted .* external jobid" "$CURRENT" 2>/dev/null \
    | awk '{{print $NF}}' | sed "s/['.]//g" | sort -u | tr "\n" " " || true)

if [ -n "$JOBIDS" ] && [ -x workflow/scripts/jobby ]; then
    ./workflow/scripts/jobby $JOBIDS > "job_information_${{timestamp}}.tsv" \
        2>> "{log}" || : > "job_information_${{timestamp}}.tsv"
else
    : > "job_information_${{timestamp}}.tsv"
fi

grep --color=never '^jobid\|FAILED' "job_information_${{timestamp}}.tsv" \
    > "failed_jobs_${{timestamp}}.tsv" 2>/dev/null \
    || : > "failed_jobs_${{timestamp}}.tsv"
"""


onstart:
    shell(
        'rm -f COMPLETED FAILED RUNNING\n'
        ': > "%s"\n'
        'echo "%s $(date +\'%%Y-%%m-%%d %%H:%%M:%%S\')'
        '  SLURM job ${{SLURM_JOB_ID:-unknown}} ===" >> "%s"\n'
        'touch RUNNING\n' % (_MASTER_ERR, _RUN_MARKER, _MASTER_LOG)
    )
    print("\n▶  viralrecon pipeline started.\n")


onsuccess:
    shell(_CURRENT_SLICE.replace("{log}", _MASTER_LOG)
          + _JOBBY.replace("{log}", _MASTER_LOG)
          + "\ntouch COMPLETED\n")
    print("\n✓  viralrecon pipeline completed successfully.\n"
          "   Sentinel:  COMPLETED\n"
          "   Resources: job_information_<timestamp>.tsv\n")


# On failure, master.err gets the rule-level errors Snakemake reported for this
# run, with the per-rule log path each one names, plus the SLURM failures. That
# is the whole point of keeping it empty otherwise: FAILED sentinel -> open
# master.err -> see what broke and which log to read, without scrolling a
# driver log that spans several runs.
_ERR_REPORT = r"""
{
  echo "viralrecon FAILED  $(date +'%Y-%m-%d %H:%M:%S')"
  echo "master SLURM job: ${{SLURM_JOB_ID:-unknown}}"
  echo
  if grep -q "^Error in rule" "$CURRENT" 2>/dev/null; then
      echo "Rules that failed (per-rule SLURM ids are in failed_jobs_*.tsv):"
      echo
      awk '/^Error in rule/{{p=1}} p{{print "  " $0}} /^$/{{p=0}}' "$CURRENT"
  else
      echo "No rule-level error was reported; the driver itself failed."
      echo "See logfiles/master.log for this run."
  fi
  echo
  echo "SLURM failures : failed_jobs_${{timestamp}}.tsv"
  echo "Resources      : job_information_${{timestamp}}.tsv"
  echo "Full driver log: logfiles/master.log"
} >> "{err}" 2>/dev/null || true
"""


onerror:
    shell(_CURRENT_SLICE.replace("{log}", _MASTER_LOG)
          + _JOBBY.replace("{log}", _MASTER_LOG)
          + _ERR_REPORT.replace("{err}", _MASTER_ERR)
          + "\ntouch FAILED\n")
    print("\n✗  viralrecon pipeline failed.\n"
          "   Sentinel:  FAILED\n"
          "   What broke: logfiles/master.err\n"
          "   Failures:  failed_jobs_<timestamp>.tsv\n"
          "   Resources: job_information_<timestamp>.tsv\n")
