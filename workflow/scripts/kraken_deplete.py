#!/usr/bin/env python3
"""Deplete host and contaminant reads from a FASTQ pair using Kraken2 output. Only reads
inside --deplete-taxids subtrees are dropped; unclassified and viral reads are kept."""

import argparse
import gzip
import sys


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def subtree_taxids(report_path, roots):
    """Collect every taxid at or below each root. Report nesting is two leading spaces
    per level, so a subtree ends at the first row back at or above its indent."""
    roots = set(roots)
    keep, active_depth = set(), None

    with open_maybe_gzip(report_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            taxid, name = cols[4].strip(), cols[5]
            depth = (len(name) - len(name.lstrip(" "))) // 2

            if active_depth is not None and depth <= active_depth:
                active_depth = None          # left the subtree

            if taxid in roots:
                active_depth = depth
                keep.add(taxid)
            elif active_depth is not None:
                keep.add(taxid)

    missing = roots - keep
    if missing:
        print("WARNING: taxid(s) not present in this Kraken2 report, nothing "
              "will be depleted for them: {}".format(", ".join(sorted(missing))),
              file=sys.stderr)
    return keep


# Clades reported for a viral assay. "Unclassified" is its own row, not a viral
# proxy: the target virus is in the database and comes back classified.
PROFILE_CLADES = [
    ("0",     "unclassified"),
    ("10239", "Viruses"),
    ("9606",  "Homo sapiens"),
    ("2",     "Bacteria"),
    ("2759",  "Eukaryota"),
    ("2157",  "Archaea"),
]


def clade_counts(report_path):
    """Map taxid -> clade-level read count (column 2 of the Kraken2 report)."""
    counts = {}
    with open_maybe_gzip(report_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            try:
                counts[cols[4].strip()] = int(cols[1])
            except ValueError:
                continue
    return counts


def write_profile(report_path, out_path, extra_taxids=()):
    """Write a compositional breakdown of the library before depletion."""
    counts = clade_counts(report_path)
    total = sum(counts.get(t, 0) for t in ("0", "1")) or 1

    rows = list(PROFILE_CLADES)
    for t in extra_taxids:
        if all(t != taxid for taxid, _ in rows):
            rows.append((t, "taxid {}".format(t)))

    with open(out_path, "w") as fh:
        fh.write("taxid\tclade\treads\tpct_of_total\n")
        for taxid, name in rows:
            n = counts.get(taxid, 0)
            fh.write("{}\t{}\t{}\t{:.4f}\n".format(
                taxid, name, n, 100.0 * n / total))


def fastq_records(handle):
    """Yield raw 4-line FASTQ records."""
    while True:
        header = handle.readline()
        if not header:
            return
        yield header + handle.readline() + handle.readline() + handle.readline()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kraken-output", required=True,
                    help="Kraken2 per-read output (may be gzipped)")
    ap.add_argument("--report", required=True, help="Kraken2 report")
    ap.add_argument("--in1", required=True)
    ap.add_argument("--in2", help="R2 (omit for single-end)")
    ap.add_argument("--out1", required=True)
    ap.add_argument("--out2", help="depleted R2")
    ap.add_argument("--deplete-taxids", nargs="+",
                    default=["9606", "10847", "2"],
                    help="root taxids to remove, with descendants "
                         "(default: 9606 Homo sapiens, 10847 phiX174, "
                         "2 Bacteria)")
    ap.add_argument("--summary", help="write a TSV of counts here")
    ap.add_argument("--profile", help="write a pre-depletion compositional "
                                      "breakdown (contamination profile) here")
    ap.add_argument("--profile-taxids", nargs="*", default=[],
                    help="extra taxids to include as rows in --profile, "
                         "e.g. the target virus")
    args = ap.parse_args()

    paired = bool(args.in2)
    if paired and not args.out2:
        ap.error("--out2 is required when --in2 is given")

    # Profile the library as Kraken2 saw it, i.e. before anything is removed.
    if args.profile:
        write_profile(args.report, args.profile, args.profile_taxids)

    drop_taxids = subtree_taxids(args.report, args.deplete_taxids)
    print("depleting {} taxids under {}".format(
        len(drop_taxids), " ".join(args.deplete_taxids)), file=sys.stderr)

    total = dropped = 0
    kh = open_maybe_gzip(args.kraken_output)
    i1 = open_maybe_gzip(args.in1)
    o1 = gzip.open(args.out1, "wt", compresslevel=4)
    i2 = open_maybe_gzip(args.in2) if paired else None
    o2 = gzip.open(args.out2, "wt", compresslevel=4) if paired else None

    try:
        r1 = fastq_records(i1)
        r2 = fastq_records(i2) if paired else None
        for line in kh:
            cols = line.split("\t", 3)
            if len(cols) < 3:
                continue
            # Column 3 is the taxid, unless kraken2 ran with --use-names, in
            # which case it reads "Homo sapiens (taxid 9606)".  Accept both.
            field = cols[2].strip()
            if field.endswith(")") and "(taxid " in field:
                field = field.rsplit("(taxid ", 1)[1][:-1].strip()
            rec1 = next(r1, None)
            rec2 = next(r2, None) if paired else None
            if rec1 is None:
                break
            total += 1
            if field in drop_taxids:
                dropped += 1
                continue
            o1.write(rec1)
            if paired:
                o2.write(rec2)
    finally:
        for h in (kh, i1, o1, i2, o2):
            if h:
                h.close()

    kept = total - dropped
    pct = (100.0 * dropped / total) if total else 0.0
    msg = "total={} depleted={} ({:.4f}%) kept={}".format(total, dropped, pct, kept)
    print(msg, file=sys.stderr)

    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("total_reads\tdepleted_reads\tdepleted_pct\tkept_reads\n")
            fh.write("{}\t{}\t{:.4f}\t{}\n".format(total, dropped, pct, kept))


if __name__ == "__main__":
    main()
