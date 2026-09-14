#!/usr/bin/env python3
"""Static final_report figures per target: genome overview, variant heatmap, coverage
comparison and composition, drawn only from files the pipeline already wrote."""

import argparse
import csv
import gzip
import os
import sys
from os.path import exists, join

import matplotlib
matplotlib.use("Agg")            # no display on a compute node
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np

# snpEff impact -> colour. Ordered worst-first so legends read sensibly.
IMPACT_COLOUR = {
    "HIGH":     "#c0392b",
    "MODERATE": "#e67e22",
    "LOW":      "#27ae60",
    "MODIFIER": "#7f8c8d",
    "":         "#7f8c8d",
}

DPI = 150


def read_tsv(path):
    """Rows of a TSV as dicts; empty list when the file is absent or empty."""
    if not exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def mosdepth_per_base(path):
    """(positions, depths) from a mosdepth per-base BED, one point per interval
    midpoint."""
    if not exists(path):
        return np.array([]), np.array([])
    xs, ys = [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                start, end, depth = int(f[1]), int(f[2]), float(f[3])
            except ValueError:
                continue
            xs.append((start + end) / 2.0)
            ys.append(depth)
    return np.array(xs), np.array(ys)


def lowcov_intervals(path):
    """Masked regions from the low-coverage BED, for shading."""
    out = []
    if not exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                try:
                    out.append((int(f[1]), int(f[2])))
                except ValueError:
                    pass
    return out


def matrix_rows(path, samples):
    """Variant matrix rows as (pos, gene, aa, impact, {sample: pct}). Blank cells mean
    not called and stay absent rather than zero."""
    rows = []
    for r in read_tsv(path):
        try:
            pos = int(r.get("POS", ""))
        except ValueError:
            continue
        pct = {}
        for s in samples:
            v = (r.get("%s.PCT" % s) or "").strip()
            if v:
                try:
                    pct[s] = float(v)
                except ValueError:
                    pass
        rows.append((pos, r.get("GENE", ""), r.get("AA_CHANGE", ""),
                     (r.get("SEVERITY", "") or "").upper(), pct))
    rows.sort(key=lambda t: t[0])
    return rows


def genome_length(fai):
    if not exists(fai):
        return 0
    with open(fai) as fh:
        first = fh.readline().split("\t")
    return int(first[1]) if len(first) > 1 else 0


# ── figure 1: per-sample genome overview ────────────────────────────────────

def fig_genome_overview(sample, target, work, outdir, rows, glen, min_depth):
    al = join(work, sample, "alignment", target)
    vc = join(work, sample, "variant_calling", target)
    xs, ys = mosdepth_per_base(
        join(al, "mosdepth", "%s.%s.per-base.bed.gz" % (sample, target)))
    if xs.size == 0:
        return None

    fig, (ax, axv) = plt.subplots(
        2, 1, figsize=(13, 5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

    floor = max(1.0, min_depth / 2.0 if min_depth else 1.0)
    ax.fill_between(xs, floor, np.maximum(ys, floor), color="#3498db",
                    alpha=0.55, linewidth=0)
    ax.set_yscale("log")
    # Start near the depth threshold, not 1, so the axis covers depths actually reached.
    ax.set_ylim(floor, max(np.max(ys) * 1.6, floor * 10))
    ax.set_ylabel("depth (log)")
    ax.set_title("%s  vs  %s" % (sample, target), loc="left", fontsize=11)

    if min_depth > 0:
        ax.axhline(min_depth, color="#c0392b", lw=1, ls="--", alpha=0.8)
        ax.text(0, min_depth, " consensus min depth %d" % min_depth,
                color="#c0392b", fontsize=7, va="bottom")

    # Regions written as N in the consensus. Shading them is the point: a
    # reader should see the gaps, not infer them from a coverage number.
    masked = lowcov_intervals(
        join(vc, "%s.%s.lowcov_mask.bed" % (sample, target)))
    for a, b in masked:
        ax.axvspan(a, b, color="#e74c3c", alpha=0.18, linewidth=0)

    # Variant needles, height = allele fraction, colour = predicted impact.
    for pos, gene, aa, impact, pct in rows:
        if sample not in pct:
            continue
        c = IMPACT_COLOUR.get(impact, IMPACT_COLOUR[""])
        axv.vlines(pos, 0, pct[sample], color=c, lw=1.4)
        axv.plot([pos], [pct[sample]], "o", color=c, ms=3.5)
        if pct[sample] >= 50 and (aa or gene):
            axv.annotate(aa or gene, (pos, pct[sample]), fontsize=6,
                         rotation=90, ha="center", va="bottom",
                         xytext=(0, 3), textcoords="offset points")

    axv.set_ylim(0, 155)   # headroom so rotated AA labels are not clipped
    axv.set_ylabel("alt %")
    axv.set_xlabel("position (bp)")
    if glen:
        axv.set_xlim(0, glen)
    axv.grid(axis="y", alpha=0.25, lw=0.5)

    present = {i for _, _, _, i, p in rows if sample in p}
    handles = [Patch(facecolor=IMPACT_COLOUR[i], label=i.title())
               for i in ("HIGH", "MODERATE", "LOW", "MODIFIER") if i in present]
    # Only claim a masked region exists when one does. A legend entry for
    # something absent from the plot invites the reader to look for it.
    if masked:
        handles.append(Patch(facecolor="#e74c3c", alpha=0.18,
                             label="masked (below min depth)"))
    if handles:
        # Below the axes, so it cannot sit on top of a variant needle.
        axv.legend(handles=handles, fontsize=7, ncol=len(handles),
                   loc="upper center", bbox_to_anchor=(0.5, -0.42),
                   frameon=False)

    out = join(outdir, "genome_overview.%s.png" % sample)
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


# ── figure 2: variant x sample heatmap ──────────────────────────────────────

def fig_variant_heatmap(samples, target, outdir, rows):
    if not rows:
        return None
    # Starts at white so a low percentage reads as near-empty; "not called" is grey so
    # it cannot be mistaken for a low percentage.
    cmap = LinearSegmentedColormap.from_list(
        "af", ["#ffffff", "#9ecae1", "#3182bd", "#08306b"])

    grid = np.full((len(rows), len(samples)), np.nan)
    for i, (_, _, _, _, pct) in enumerate(rows):
        for j, s in enumerate(samples):
            if s in pct:
                grid[i, j] = pct[s]

    h = max(3.0, 0.32 * len(rows) + 1.6)
    fig, ax = plt.subplots(figsize=(max(6.0, 1.5 + 1.1 * len(samples)), h))
    cmap.set_bad("#b0bec5")           # not called: grey, never 0
    im = ax.imshow(np.ma.masked_invalid(grid), aspect="auto", cmap=cmap,
                   vmin=0, vmax=100)

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(
        ["%d  %s %s" % (p, g, a) if (g or a) else str(p)
         for p, g, a, _, _ in rows], fontsize=8)

    for i in range(len(rows)):
        for j in range(len(samples)):
            if not np.isnan(grid[i, j]):
                ax.text(j, i, "%.0f" % grid[i, j], ha="center", va="center",
                        fontsize=7,
                        color="white" if grid[i, j] > 55 else "#263238")

    ax.set_title("%s - filtered variants (%% alt reads; grey = not called)"
                 % target, fontsize=10, loc="left")
    fig.colorbar(im, ax=ax, shrink=0.7, label="alt %")
    out = join(outdir, "variant_heatmap.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


# ── figure 3: coverage of every sample together ─────────────────────────────

def fig_coverage_comparison(samples, target, work, outdir, glen, min_depth):
    fig, ax = plt.subplots(figsize=(13, 4.2))
    drawn = 0
    for s in samples:
        xs, ys = mosdepth_per_base(
            join(work, s, "alignment", target, "mosdepth",
                 "%s.%s.per-base.bed.gz" % (s, target)))
        if xs.size:
            ax.plot(xs, np.maximum(ys, 1), lw=0.9, alpha=0.85, label=s)
            drawn += 1
    if not drawn:
        plt.close(fig)
        return None

    if min_depth > 0:
        ax.axhline(min_depth, color="#c0392b", lw=1, ls="--", alpha=0.8,
                   label="min depth %d" % min_depth)
    ax.set_yscale("log")
    ax.set_xlabel("position (bp)")
    ax.set_ylabel("depth (log)")
    if glen:
        ax.set_xlim(0, glen)
    ax.set_title("%s - depth, all samples" % target, fontsize=10, loc="left")
    ax.legend(fontsize=7, ncol=max(1, min(6, drawn)), loc="lower center")
    ax.grid(alpha=0.25, lw=0.5)
    out = join(outdir, "coverage_comparison.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


# ── figure 4: library composition ───────────────────────────────────────────

def fig_composition(samples, work, outdir, target_taxids):
    """Library composition per sample from the Kraken2 composition TSVs (run_summary
    would make the DAG cyclic). The target is split out of Viruses, not counted twice."""
    per = {}
    for s_ in samples:
        comp = join(work, s_, "pre_process", "kraken2",
                    "%s.kraken2_decon.composition.tsv" % s_)
        rows = read_tsv(comp)
        if rows:
            per[s_] = {r.get("taxid", ""): r.get("pct_of_total", "")
                       for r in rows}
    order = [s_ for s_ in samples if s_ in per]
    if not order:
        return None

    def num(d, taxid):
        try:
            return max(0.0, float(d.get(taxid) or 0))
        except ValueError:
            return 0.0

    target, other_viral, human, other = [], [], [], []
    for s_ in order:
        d = per[s_]
        # Each target's own taxid, so this is not hardcoded to one virus.
        tv = sum(num(d, t) for t in target_taxids) if target_taxids else 0.0
        vir = num(d, "10239")              # Viruses superkingdom: a superset
        hu = num(d, "9606")
        ov = max(0.0, vir - tv)
        rest = max(0.0, 100.0 - tv - ov - hu)
        target.append(tv); other_viral.append(ov); human.append(hu); other.append(rest)

    if not any(target) and not any(other_viral) and not any(human):
        return None

    fig, ax = plt.subplots(figsize=(max(6.0, 1.3 * len(order)), 4.4))
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for vals, label, colour in (
            (target,      "target virus",        "#2ecc71"),
            (other_viral, "other viral",         "#16a085"),
            (human,       "human",               "#e67e22"),
            (other,       "other / unclassified", "#bdc3c7")):
        vals = np.array(vals)
        ax.bar(x, vals, bottom=bottom, label=label, color=colour, width=0.72)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of classified reads")
    ax.set_title("library composition, before depletion", fontsize=10, loc="left")
    ax.legend(fontsize=7, loc="upper right", ncol=2, framealpha=0.9)
    ax.grid(axis="y", alpha=0.25, lw=0.5)
    out = join(outdir, "composition.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description="Static figures for final_report")
    ap.add_argument("--workpath", required=True)
    ap.add_argument("--outdir", required=True, help="final_report/<target>/figures")
    ap.add_argument("--target", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--matrix", required=True, help="filtered variants matrix")
    ap.add_argument("--target-taxids", nargs="*", default=[],
                    help="taxids of this target, for the composition figure")
    ap.add_argument("--consensus-min-depth", type=int, default=10)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    rows = matrix_rows(a.matrix, a.samples)
    glen = genome_length(join(a.workpath, "ref_db", a.target,
                              "%s.fa.fai" % a.target))
    made = []
    for s in a.samples:
        p = fig_genome_overview(s, a.target, a.workpath, a.outdir, rows,
                                glen, a.consensus_min_depth)
        if p:
            made.append(p)
    for p in (fig_variant_heatmap(a.samples, a.target, a.outdir, rows),
              fig_coverage_comparison(a.samples, a.target, a.workpath,
                                      a.outdir, glen, a.consensus_min_depth),
              fig_composition(a.samples, a.workpath, a.outdir,
                              a.target_taxids)):
        if p:
            made.append(p)

    for p in made:
        print("  figure: %s" % os.path.basename(p))
    print("wrote %d figure(s) to %s" % (len(made), a.outdir))
    if not made:
        print("no figures produced; inputs were missing or empty",
              file=sys.stderr)


if __name__ == "__main__":
    main()
