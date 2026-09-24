#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Parse samtools flagstat outputs and write a mapping summary TSV for one target."""

import argparse
import os
import re
import sys


def parse_flagstat(path):
    metrics = {
        "total_reads": 0,
        "mapped_reads": 0,
        "mapped_pct": "0.00%",
        "properly_paired": 0,
    }
    try:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip()
                m = re.match(r"(\d+) \+ \d+ in total", line)
                if m:
                    metrics["total_reads"] = int(m.group(1))
                m = re.match(r"(\d+) \+ \d+ mapped \((.+?)\)", line)
                if m:
                    metrics["mapped_reads"] = int(m.group(1))
                    metrics["mapped_pct"] = m.group(2).split(":")[0].strip()
                m = re.match(r"(\d+) \+ \d+ properly paired", line)
                if m:
                    metrics["properly_paired"] = int(m.group(1))
    except FileNotFoundError:
        pass
    return metrics


def sample_name(flagstat_path, target):
    """Sample name from .../{sample}.{target}.bowtie2_map.flagstat. The target is
    stripped by exact match because target names contain dots."""
    base = os.path.basename(flagstat_path)
    for suf in (".flagstat", ".bowtie2_map"):
        base = base.rsplit(suf, 1)[0]
    if base.endswith("." + target):
        base = base[: -(len(target) + 1)]
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--flagstats", nargs="+", required=True,
                    metavar="FILE", help="samtools flagstat files")
    ap.add_argument("--output", required=True,
                    metavar="TSV", help="output TSV path")
    ap.add_argument("--target", required=True,
                    metavar="NAME", help="viral target name")
    args = ap.parse_args()

    rows = []
    for fs in sorted(args.flagstats):
        sample = sample_name(fs, args.target)
        metrics = parse_flagstat(fs)
        rows.append({
            "sample": sample,
            "target": args.target,
            "total_reads": metrics["total_reads"],
            "mapped_reads": metrics["mapped_reads"],
            "mapped_pct": metrics["mapped_pct"],
            "properly_paired": metrics["properly_paired"],
        })

    header = ["sample", "target", "total_reads",
              "mapped_reads", "mapped_pct", "properly_paired"]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    print(f"[collect_report] Wrote {len(rows)} rows -> {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
