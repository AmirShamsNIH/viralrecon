#!/usr/bin/env python3
"""
make_variants_matrix.py
───────────────────────
Reshape a GATK VariantsToTable output into a variant x sample matrix: one row
per variant, one column pair per sample carrying allele depth and the alt
fraction as a percentage.

Why this exists alongside the long table
----------------------------------------
The long table is one row per sample x variant, which is the right shape for
scripts and plots but the wrong shape for the question people actually ask of
a multi-sample viral run: "is this variant fixed in every sample, or is it a
minor allele in three of them?" Answering that from long format means pivoting
it by hand. Here it is one row, read left to right.

An empty cell is meaningful and is deliberately left empty rather than filled
with 0 or NA: the variant was not called in that sample. A "0.0" would claim
the site was examined and found reference, which is a stronger statement than
the VCF supports.

The layout follows the RML covid-seq final_filtered_variants.txt convention,
which this reimplements from a merged VCF rather than a per-sample cat.

Input columns (GATK VariantsToTable, requires -F ANN -F TYPE):
    CHROM  POS  TYPE  REF  ALT  ANN  SAMPLE1.AD  SAMPLE1.DP  SAMPLE2.AD ...

Output columns:
    CHROM POS TYPE REF ALT EFFECT SEVERITY GENE HGVS_C HGVS_P AA_CHANGE
    <sample>.AD  <sample>.PCT   (one pair per sample)
"""

import argparse
import csv
import os
import sys

# snpEff ANN subfields, pipe-delimited, in the order the spec defines them.
# Only the ones that end up as columns are named here.
_ANN_ALLELE, _ANN_EFFECT, _ANN_IMPACT, _ANN_GENE = 0, 1, 2, 3
_ANN_HGVS_C, _ANN_HGVS_P = 9, 10

# Three-letter to one-letter amino acid codes, so p.Gln6249His also appears as
# Q6249H -- the form people search for and cite.
_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O",
}

MISSING = (".", "", "NA", "./.")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Reshape wide GATK VariantsToTable -> variant x sample matrix"
    )
    ap.add_argument("--input",   required=True, help="GATK VariantsToTable TSV")
    ap.add_argument("--output",  required=True, help="Output matrix TSV")
    ap.add_argument("--samples", nargs="+", required=True, help="Sample names")
    ap.add_argument("--min-pct", type=float, default=0.0,
                    help="Blank out a sample's cell below this alt percentage "
                         "(default 0.0, i.e. report whatever was called)")
    return ap.parse_args()


def aa_one_letter(hgvs_p):
    """
    Convert p.Gln6249His to Q6249H. Returns "" when there is no protein change
    or the notation is not a simple substitution -- deletions and frameshifts
    are left to the HGVS_P column rather than mangled into a shorter form that
    would misrepresent them.
    """
    if not hgvs_p or hgvs_p in MISSING:
        return ""
    s = hgvs_p[2:] if hgvs_p.startswith("p.") else hgvs_p
    # Only a plain <AA3><pos><AA3> substitution converts cleanly.
    if len(s) < 7:
        return ""
    ref3, rest = s[:3], s[3:]
    if ref3 not in _AA3TO1:
        return ""
    digits = ""
    i = 0
    while i < len(rest) and rest[i].isdigit():
        digits += rest[i]
        i += 1
    alt3 = rest[i:]
    if not digits or alt3 not in _AA3TO1:
        return ""
    return "%s%s%s" % (_AA3TO1[ref3], digits, _AA3TO1[alt3])


def parse_ann(ann_field):
    """
    Pull the first (highest-ranked) annotation out of a snpEff ANN field.

    snpEff orders annotations by predicted severity, so the first entry is the
    one a reviewer wants to see. The rest stay in the VCF; duplicating every
    transcript here would defeat the point of a one-row-per-variant table.
    """
    blank = {"EFFECT": "", "SEVERITY": "", "GENE": "",
             "HGVS_C": "", "HGVS_P": "", "AA_CHANGE": ""}
    if not ann_field or ann_field in MISSING:
        return blank
    # GATK renders a multi-value INFO field as a comma-separated list.
    first = ann_field.split(",")[0]
    parts = first.split("|")

    def get(i):
        return parts[i].strip() if i < len(parts) else ""

    hgvs_p = get(_ANN_HGVS_P)
    return {
        "EFFECT":    get(_ANN_EFFECT),
        "SEVERITY":  get(_ANN_IMPACT),
        "GENE":      get(_ANN_GENE),
        "HGVS_C":    get(_ANN_HGVS_C),
        "HGVS_P":    hgvs_p,
        "AA_CHANGE": aa_one_letter(hgvs_p),
    }


def alt_percent(ad):
    """
    Alt fraction as a percentage from an AD field of the form "ref,alt[,alt2]".

    Multi-allelic records are split upstream by `bcftools norm -m -any`, so a
    third value is not expected; if one appears, every non-reference count is
    summed rather than silently dropped.

    Returns None when there is no usable depth -- including total depth 0,
    which is a site with no reads rather than a site with no alt reads.
    """
    if not ad or ad in MISSING:
        return None
    try:
        counts = [int(x) for x in ad.split(",")]
    except ValueError:
        return None
    if len(counts) < 2:
        return None
    total = sum(counts)
    if total <= 0:
        return None
    return 100.0 * sum(counts[1:]) / total


def build(input_path, output_path, samples, min_pct):
    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)

    shared = ["CHROM", "POS", "TYPE", "REF", "ALT"]
    ann_cols = ["EFFECT", "SEVERITY", "GENE", "HGVS_C", "HGVS_P", "AA_CHANGE"]

    out_header = list(shared) + ann_cols
    for s in samples:
        out_header += ["%s.AD" % s, "%s.PCT" % s]

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    written = 0
    with open(output_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()

        for row in rows:
            rec = {c: row.get(c, "") for c in shared}
            # TYPE is only present when VariantsToTable was given -F TYPE.
            if not rec.get("TYPE"):
                rec["TYPE"] = ""
            rec.update(parse_ann(row.get("ANN", "")))

            any_called = False
            for s in samples:
                ad = row.get("%s.AD" % s, "")
                pct = alt_percent(ad)
                if pct is None or pct < min_pct:
                    rec["%s.AD" % s] = ""
                    rec["%s.PCT" % s] = ""
                    continue
                rec["%s.AD" % s] = ad
                rec["%s.PCT" % s] = "%.1f" % pct
                any_called = True

            # A variant present in the merged VCF but called in no sample after
            # the percentage floor carries no information; drop the row rather
            # than emit one that is blank across every sample column.
            if not any_called:
                continue
            w.writerow(rec)
            written += 1

    print("Wrote variant matrix: %d variants x %d samples -> %s"
          % (written, len(samples), output_path), file=sys.stderr)


def main():
    a = parse_args()
    build(a.input, a.output, a.samples, a.min_pct)


if __name__ == "__main__":
    main()
