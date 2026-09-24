# hooks.smk: keep one sentinel (RUNNING, COMPLETED or FAILED) in the run directory and
# run jobby over this run's SLURM jobs on exit, whether it succeeded or failed.

# master.log is appended across runs and master.err holds failures only.
# _CURRENT_SLICE keeps the text after onstart's marker, so jobby sees this run only.
_MASTER_LOG = join("logfiles", "master.log")
_MASTER_ERR = join("logfiles", "master.err")
_RUN_MARKER = "=== viralrecon run started"

_CURRENT_SLICE = r"""
CURRENT=$(mktemp)
awk -v m="{marker}" 'index($0, m)==1 {{buf=""}} {{buf = buf $0 ORS}} END {{printf "%s", buf}}' "{log}" > "$CURRENT" 2>/dev/null || : > "$CURRENT"
""".replace("{marker}", _RUN_MARKER)

_JOBBY = r"""
sleep 15
rm -f COMPLETED FAILED RUNNING
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# SLURM ids this run submitted; none is normal when the DAG was already satisfied.
# `|| true` stops grep's exit 1 on no match from failing the hook under set -e.
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
    print("\nviralrecon pipeline started.\n")


onsuccess:
    shell(_CURRENT_SLICE.replace("{log}", _MASTER_LOG)
          + _JOBBY.replace("{log}", _MASTER_LOG)
          + "\ntouch COMPLETED\n")
    print("\nviralrecon pipeline completed successfully.\n"
          "   Sentinel:  COMPLETED\n"
          "   Resources: job_information_<timestamp>.tsv\n")


# On failure, master.err gets this run's rule errors with their log paths plus the
# SLURM failures, so FAILED leads straight to what broke.
_ERR_REPORT = r"""
{{
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
}} >> "{err}" 2>/dev/null || true
"""


onerror:
    shell(_CURRENT_SLICE.replace("{log}", _MASTER_LOG)
          + _JOBBY.replace("{log}", _MASTER_LOG)
          + _ERR_REPORT.replace("{err}", _MASTER_ERR)
          + "\ntouch FAILED\n")
    print("\nviralrecon pipeline failed.\n"
          "   Sentinel:  FAILED\n"
          "   What broke: logfiles/master.err\n"
          "   Failures:  failed_jobs_<timestamp>.tsv\n"
          "   Resources: job_information_<timestamp>.tsv\n")
