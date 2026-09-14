#!/usr/bin/env python3
"""Assemble final_report/, the one directory to open after a run, from outputs that
already exist elsewhere in the run tree."""

import argparse
import csv
import glob
import gzip
import os
import shutil
from os.path import basename, exists, join


# ── small readers ────────────────────────────────────────────────────────────

def _copy(src, dstdir, newname=None):
    if not exists(src):
        return None
    os.makedirs(dstdir, exist_ok=True)
    dst = join(dstdir, newname or basename(src))
    shutil.copy2(src, dst)
    return dst


def read_tsv_rows(path):
    if not exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kraken_composition(path):
    """taxid -> pct, from the kraken2_decon composition TSV."""
    out = {}
    for r in read_tsv_rows(path):
        out[r.get("taxid", "")] = r.get("pct_of_total", "")
    return out


def mosdepth_mean(path):
    for r in read_tsv_rows(path):
        if r.get("chrom") == "total":
            return r.get("mean", "")
    rows = read_tsv_rows(path)
    return rows[0].get("mean", "") if rows else ""


def count_ns(fasta):
    if not exists(fasta):
        return ""
    n = total = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            s = line.strip()
            total += len(s)
            n += s.upper().count("N")
    return str(n)


def genome_length(fai_path):
    """Reference length from the .fa.fai (sum of all sequence lengths)."""
    total = 0
    for r in read_tsv_rows(fai_path):
        pass  # .fai has no header; parsed manually below
    if not exists(fai_path):
        return 0
    with open(fai_path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 1:
                try:
                    total += int(f[1])
                except ValueError:
                    pass
    return total


def stub_reason(path, marker):
    """True when a report file is one of our stubs rather than a real result."""
    if not exists(path):
        return False
    try:
        with open(path, errors="ignore") as fh:
            return marker in fh.read(4096)
    except Exception:
        return False


def lineage_failures(next_):
    """Name the lineage caller if it ran and left a stub for this sample x target.
    A missing file is a deliberate skip (no dataset), not a failure."""
    failed = []
    if stub_reason(next_, "nextclade_error"):
        failed.append("NEXTCLADE_FAILED")
    return failed


def target_failures(FR, target):
    """Name the target-level artefacts that ran and left a stub, such as an OOM-killed
    IGV report, so run_summary does not read "pass" over a failure."""
    failed = []
    igv = join(FR, target, "igv_report.%s.html" % target)
    if stub_reason(igv, "IGV report unavailable"):
        failed.append("IGV_REPORT_FAILED")
    figdir = join(FR, target, "figures")
    if exists(figdir) and not glob.glob(join(figdir, "*.png")):
        failed.append("FIGURES_MISSING")
    return failed


def nextclade_call(path):
    rows = read_tsv_rows(path)
    if not rows:
        return "", "", ""
    r = rows[0]
    return (r.get("clade", ""), r.get("qc.overallStatus", ""),
            r.get("coverage", ""))



def flagstat_mapped(path):
    """(mapped_reads, mapped_pct) from samtools flagstat."""
    if not exists(path):
        return "", ""
    with open(path) as fh:
        for line in fh:
            if " mapped (" in line and "primary" not in line:
                n = line.split()[0]
                pct = line.split("(")[1].split("%")[0]
                return n, pct
    return "", ""


def count_variants(vcf_gz, sample):
    """Variants with a non-reference genotype for this sample."""
    if not exists(vcf_gz):
        return ""
    try:
        with gzip.open(vcf_gz, "rt") as fh:
            idx, n = None, 0
            for line in fh:
                if line.startswith("##"):
                    continue
                f = line.rstrip("\n").split("\t")
                if line.startswith("#CHROM"):
                    idx = f.index(sample) if sample in f else None
                    if idx is None:
                        return ""
                    continue
                gt = f[idx].split(":")[0]
                if gt not in ("0/0", "0|0", "./.", ".|.", "."):
                    n += 1
        return str(n)
    except Exception:
        return ""


# ── main ─────────────────────────────────────────────────────────────────────

COLUMNS = [
    "sample", "target",
    "input_read_pairs", "pct_viral", "pct_target", "pct_human", "pct_depleted",
    "mapped_reads", "mapped_pct", "reads_used", "mean_depth", "qc_status",
    "genome_length", "consensus_masked_bases", "pct_genome_covered",
    "n_variants",
    "nextclade_clade", "nextclade_qc", "nextclade_coverage",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workpath", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--lineage-targets", nargs="*", default=[],
                    help="targets eligible for lineage calling")
    ap.add_argument("--target-taxids", nargs="*", default=[], metavar="NAME=TAXID",
                    help="NCBI taxid per target, from genome.json. pct_target is "
                         "read from the Kraken2 report at that taxid, so the "
                         "column means the same thing for every virus; a target "
                         "registered without a taxid leaves it blank rather than "
                         "reporting some other organism's percentage.")
    ap.add_argument("--min-genome-coverage", type=float, default=0.80,
                    help="fraction of the reference that must reach the "
                         "consensus depth floor before a target is trusted")
    ap.add_argument("--consensus-min-depth", default="10")
    a = ap.parse_args()

    W, FR = a.workpath, a.outdir
    lineage_targets = set(a.lineage_targets)
    target_taxid = {}
    for pair in a.target_taxids:
        name, _, taxid = pair.partition("=")
        if name and taxid:
            target_taxid[name] = taxid
    os.makedirs(FR, exist_ok=True)
    rows = []

    for target in a.targets:
        tdir = join(FR, target)
        tgt_failures = target_failures(FR, target)
        cons_d, lin_d = join(tdir, "consensus"), join(tdir, "lineage")
        var_d,  qc_d  = join(tdir, "variants"),  join(tdir, "qc")
        for d in (cons_d, var_d, qc_d):
            os.makedirs(d, exist_ok=True)
        # lineage/ only for targets the stage applies to: an empty directory reads as
        # a failure, and run_summary.tsv says why it is absent.
        if target in lineage_targets:
            os.makedirs(lin_d, exist_ok=True)

        # Aggregates are copied, never moved: they are Snakemake outputs, and moving
        # them would make the next run redo bcftools_merge.
        for pat in ("aggregate.%s.vcf.gz*", "aggregate.%s.snpeff.vcf.gz*",
                    "aggregate.%s.snpeff.variants.txt",
                    "aggregate.%s.variants_long.tsv",
                    "aggregate.%s.variants_matrix.raw.tsv",
                    "aggregate.%s.variants_matrix.filtered.tsv",
                    "aggregate.%s.filtered.vcf.gz*",
                    "aggregate.%s.filtered.variants.txt",
                    "aggregate.%s.snpEff_summary.html",
                    "aggregate.%s.snpEff_summary.genes.txt"):
            for f in glob.glob(join(tdir, pat % target)):
                if os.path.isfile(f):
                    _copy(f, var_d)
        _copy(join(tdir, "aggregate.%s.mapping_summary.tsv" % target), qc_d)

        agg_vcf = join(var_d, "aggregate.%s.snpeff.vcf.gz" % target)
        combined = join(cons_d, "all_samples.%s.consensus.fa" % target)
        cons_parts = []

        for s in a.samples:
            vc = join(W, s, "variant_calling", target)
            al = join(W, s, "alignment", target)
            k2 = join(W, s, "pre_process", "kraken2")
            ln = join(W, s, "lineage", target)

            cons = join(vc, "%s.%s.consensus.fa" % (s, target))
            if _copy(cons, cons_d):
                cons_parts.append(cons)
            _copy(join(vc, "%s.%s.lowcov_mask.bed" % (s, target)), qc_d)

            comp = join(k2, "%s.kraken2_decon.composition.tsv" % s)
            _copy(comp, qc_d)
            _copy(join(k2, "%s.kraken2_decon.krona.html" % s), qc_d)
            _copy(join(al, "mosdepth", "%s.%s.mosdepth.summary.txt" % (s, target)), qc_d)
            _copy(join(al, "%s.%s.bowtie2_map.raw.flagstat" % (s, target)), qc_d)

            next_ = join(ln, "%s.%s.nextclade.tsv" % (s, target))
            if target in lineage_targets:
                _copy(next_, lin_d)

            k = kraken_composition(comp)
            summ = read_tsv_rows(join(k2, "%s.kraken2_decon.summary.tsv" % s))
            # Mapping rate comes from the PRE-filter flagstat: after unmapped
            # records are dropped every later flagstat reads ~100% mapped.
            mapped, mpct = flagstat_mapped(
                join(al, "%s.%s.bowtie2_map.raw.flagstat" % (s, target)))
            used, _ = flagstat_mapped(
                join(al, "%s.%s.bowtie2_map.flagstat" % (s, target)))
            # Breadth, not just depth: cross-mapping reads can pile into conserved
            # regions of a divergent reference while most of it stays uncovered.
            glen   = genome_length(join(W, "ref_db", target, "%s.fa.fai" % target))
            masked = count_ns(cons)
            covered = ""
            reasons = []
            if exists(join(al, "%s.%s.qc_warn" % (s, target))):
                reasons.append("LOW_MAPPED_READS")
            if glen and masked != "":
                frac = 1.0 - (int(masked) / float(glen))
                covered = "%.4f" % frac
                if frac < a.min_genome_coverage:
                    reasons.append("LOW_GENOME_COVERAGE")
            reasons.extend(lineage_failures(next_))
            reasons.extend(tgt_failures)
            qc = "pass" if not reasons else "WARN:" + "+".join(reasons)
            nc, nq, ncov = nextclade_call(next_)

            rows.append({
                "sample": s, "target": target,
                "input_read_pairs": summ[0].get("total_reads", "") if summ else "",
                "pct_viral":      k.get("10239", ""),
                "pct_target":     k.get(target_taxid.get(target, ""), ""),
                "pct_human":      k.get("9606", ""),
                "pct_depleted":   summ[0].get("depleted_pct", "") if summ else "",
                "mapped_reads": mapped, "mapped_pct": mpct,
                "reads_used": used, "qc_status": qc,
                "mean_depth": mosdepth_mean(
                    join(al, "mosdepth", "%s.%s.mosdepth.summary.txt" % (s, target))),
                "genome_length": str(glen) if glen else "",
                "consensus_masked_bases": masked,
                "pct_genome_covered": covered,
                "n_variants": count_variants(agg_vcf, s),
                "nextclade_clade": nc, "nextclade_qc": nq,
                "nextclade_coverage": ncov,
            })

        if cons_parts:
            with open(combined, "w") as out:
                for f in cons_parts:
                    text = open(f).read()
                    out.write(text if text.endswith("\n") else text + "\n")

        if target in lineage_targets:
            lin_summary = join(lin_d, "lineage_summary.tsv")
            with open(lin_summary, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["sample", "nextclade_clade", "nextclade_qc",
                            "nextclade_coverage"])
                for r in rows:
                    if r["target"] == target:
                        w.writerow([r["sample"], r["nextclade_clade"],
                                    r["nextclade_qc"], r["nextclade_coverage"]])

    with open(join(FR, "run_summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print("final_report assembled: %d rows across %d target(s)"
          % (len(rows), len(a.targets)))

    # Log what qc_status records: a stubbed caller still ends COMPLETED, so the log
    # has to carry the warning the exit code no longer does.
    warned = [r for r in rows if r["qc_status"] != "pass"]
    if warned:
        tally = {}
        for r in warned:
            for reason in r["qc_status"].split(":", 1)[1].split("+"):
                tally[reason] = tally.get(reason, 0) + 1
        print("WARNING: %d of %d sample x target rows are not 'pass':"
              % (len(warned), len(rows)))
        for reason, n in sorted(tally.items(), key=lambda kv: -kv[1]):
            print("  %-24s %d" % (reason, n))
        for r in warned:
            print("  %s / %s -> %s" % (r["sample"], r["target"], r["qc_status"]))
        print("  See run_summary.tsv column qc_status, and logfiles/ for detail.")
    else:
        print("All %d sample x target rows passed QC." % len(rows))


if __name__ == "__main__":
    main()
