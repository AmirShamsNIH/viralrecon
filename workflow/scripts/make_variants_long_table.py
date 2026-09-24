#!/usr/bin/env python3
"""Reshape a wide GATK VariantsToTable output into a long TSV with one row per
sample x variant."""

import argparse
import csv
import os
import sys


def parse_args():
    ap = argparse.ArgumentParser(
        description="Reshape wide GATK VariantsToTable -> long-format TSV"
    )
    ap.add_argument("--input", required=True, help="GATK VariantsToTable TSV")
    ap.add_argument("--output", required=True, help="Output long-format TSV")
    ap.add_argument("--target", required=True, help="Target/reference name")
    ap.add_argument("--samples", nargs="+", required=True, help="Sample names")
    return ap.parse_args()


def _shared_cols(header, sample_names):
    """Return column names that are NOT per-sample (no 'SAMPLE.' prefix)."""
    sample_prefixes = tuple("{}.".format(s) for s in sample_names)
    return [c for c in header if not c.startswith(sample_prefixes)]


def _sample_cols(header, sample):
    """Return (header_name, local_name) pairs for a given sample."""
    prefix = "{}.".format(sample)
    return [(c, c[len(prefix):]) for c in header if c.startswith(prefix)]


def reshape(input_path, output_path, target, samples):
    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)

    shared = _shared_cols(header, samples)

    out_header = ["TARGET", "SAMPLE"] + shared
    local_col_sets = [set(lc for _, lc in _sample_cols(header, s)) for s in samples]
    all_local = []
    seen = set()
    for s in samples:
        for _, lc in _sample_cols(header, s):
            if lc not in seen:
                all_local.append(lc)
                seen.add(lc)
    out_header += all_local

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=out_header, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()

        for row in rows:
            for sample in samples:
                s_cols = _sample_cols(header, sample)
                # Skip rows where all sample-specific fields are missing / '.'
                if all(row.get(hc, ".") in (".", "", "NA") for hc, _ in s_cols):
                    continue

                out_row = {"TARGET": target, "SAMPLE": sample}
                for col in shared:
                    out_row[col] = row.get(col, ".")
                for hc, lc in s_cols:
                    out_row[lc] = row.get(hc, ".")
                writer.writerow(out_row)

    print("Wrote long-format table -> {}".format(output_path), file=sys.stderr)


def main():
    args = parse_args()
    reshape(args.input, args.output, args.target, args.samples)


if __name__ == "__main__":
    main()
