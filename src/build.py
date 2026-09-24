#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Download, index and register a viral reference genome in genome.json."""

import csv
import json
import os
import shutil
import subprocess
import sys

try:
    from .utils import err, fatal, which
    from . import containers
except ImportError:
    from utils import err, fatal, which
    import containers

# Built-in virus to accession presets

VIRUS_PRESETS = {
    "SARS":           "NC_045512.2",   # SARS-CoV-2 Wuhan-Hu-1
    "SARS2":          "NC_045512.2",   # alias
    "MPOX":           "NC_063383.1",   # Mpox virus
    "HIV1":           "NC_001802.1",   # HIV-1
    "INFLUENZA_H1N1": "NC_026433.1",   # Influenza A H1N1 HA segment
    "EBV":            "NC_007605.1",   # Epstein-Barr virus
    "RSV":            "NC_038235.1",   # RSV-A
    "DENGUE1":        "NC_001477.1",   # Dengue virus 1
}


def _canonical_name(virus, accession):
    """VIRUS_ACCESSION, the target name used throughout the pipeline."""
    return "{}_{}".format(virus.upper().replace(" ", "_"), accession)


# Download helpers

def _curl_or_wget(url, out):
    if which("curl"):
        subprocess.check_call(["curl", "-fsSL", url, "-o", out])
    elif which("wget"):
        subprocess.check_call(["wget", "-q", url, "-O", out])
    else:
        fatal("Neither curl nor wget found: cannot download reference files.")


def _efetch_fasta(accession, out_path):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    url  = "{}/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text".format(
        base, accession)
    print("  Downloading FASTA for {} …".format(accession))
    _curl_or_wget(url, out_path)
    with open(out_path) as fh:
        if fh.read(1) != ">":
            os.remove(out_path)
            fatal("FASTA download for {} produced unexpected content. "
                  "Check the accession and try again.".format(accession))


def _efetch_taxid(accession):
    """NCBI taxid for an accession, or None. Recorded so the Kraken2 profile can
    report each target's own share of the library."""
    import tempfile
    import time

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    url  = "{}/esummary.fcgi?db=nuccore&id={}&retmode=json".format(base, accession)

    # NCBI throttles anonymous callers to about 3 requests/s and a build has just
    # made two efetch calls, so retry with backoff rather than lose the taxid.
    for attempt in range(3):
        tmp = None
        try:
            fd, tmp = tempfile.mkstemp(suffix=".json")
            os.close(fd)
            _curl_or_wget(url, tmp)
            with open(tmp) as fh:
                data = json.load(fh)
            result = data.get("result", {})
            for uid in result.get("uids", []):
                taxid = result.get(uid, {}).get("taxid")
                if taxid:
                    return str(taxid)
        except Exception:
            pass
        finally:
            if tmp and os.path.isfile(tmp):
                os.remove(tmp)
        time.sleep(1 + attempt)
    return None


def _efetch_gff(accession, out_path):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    url  = "{}/efetch.fcgi?db=nuccore&id={}&rettype=gff3&retmode=text".format(
        base, accession)
    print("  Downloading GFF3 for {} …".format(accession))
    _curl_or_wget(url, out_path)
    with open(out_path) as fh:
        if fh.read(1) != "#":
            os.remove(out_path)
            fatal("GFF3 download for {} produced unexpected content. "
                  "Check the accession and try again.".format(accession))


def _download_reference(accession, genome_dir, canonical_name):
    """Download FASTA + GFF into genome_dir (skips files that already exist)."""
    os.makedirs(genome_dir, exist_ok=True)
    fasta_out = os.path.join(genome_dir, "{}.fa".format(canonical_name))
    gff_out   = os.path.join(genome_dir, "genes.gff")

    if os.path.isfile(fasta_out) and os.path.getsize(fasta_out) > 0:
        print("  FASTA already present, skipping download")
    else:
        _efetch_fasta(accession, fasta_out)

    if os.path.isfile(gff_out) and os.path.getsize(gff_out) > 0:
        print("  GFF already present, skipping download")
    else:
        _efetch_gff(accession, gff_out)

    missing = []
    if not os.path.isfile(fasta_out) or os.path.getsize(fasta_out) == 0:
        missing.append("FASTA ({})".format(fasta_out))
    if not os.path.isfile(gff_out) or os.path.getsize(gff_out) == 0:
        missing.append("GFF ({})".format(gff_out))
    if missing:
        fatal("Failed to obtain: {}".format(", ".join(missing)))

    print("  FASTA → {}".format(fasta_out))
    print("  GFF   → {}".format(gff_out))
    return os.path.abspath(fasta_out), os.path.abspath(gff_out)


# Local file validation

def _fasta_seqs(path):
    """{sequence id: length} from a FASTA. Ids are the first whitespace token."""
    seqs, name, length = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = length
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        seqs[name] = length
    return seqs


def _annotation_features(path):
    """(seqid, end) per feature row of a GFF3/GTF. Both share columns 1 and 5."""
    feats = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            try:
                feats.append((f[0], int(f[4])))
            except ValueError:
                continue
    return feats


def _validate_local(fasta, annotation):
    """Check the FASTA and annotation describe the same genome. snpEff runs with
    -noCheckCds -noCheckProtein, so a mismatched pair would annotate silently wrong."""
    errors = []
    if not os.path.isfile(fasta):
        errors.append("FASTA not found: {}".format(fasta))
    else:
        with open(fasta) as fh:
            if fh.read(1) != ">":
                errors.append("Does not look like a FASTA: {}".format(fasta))
    if not os.path.isfile(annotation):
        errors.append("Annotation not found: {}".format(annotation))
    if errors:
        fatal("\n\t" + "\n\t".join(errors))

    seqs  = _fasta_seqs(fasta)
    feats = _annotation_features(annotation)
    if not seqs:
        fatal("\n\tNo sequences found in {}".format(fasta))
    if not feats:
        fatal(
            "\n\tNo feature rows found in {}\n"
            "\tExpected GFF3 or GTF with tab-separated columns."
            .format(annotation)
        )

    unknown = sorted({sid for sid, _ in feats} - set(seqs))
    if unknown:
        fatal(
            "\n\tAnnotation refers to sequences that are not in the FASTA:\n"
            "\t  annotation: {}\n"
            "\t  FASTA has : {}\n"
            "\tThe two files describe different genomes, or the FASTA headers\n"
            "\thave been renamed. snpEff would annotate against coordinates\n"
            "\tthat do not belong to this sequence and report the result with\n"
            "\tfull confidence."
            .format(", ".join(unknown[:5]), ", ".join(sorted(seqs)[:5]))
        )

    over = [(sid, end, seqs[sid]) for sid, end in feats if end > seqs[sid]]
    if over:
        sid, end, ln = over[0]
        fatal(
            "\n\tAnnotation has {} feature(s) running past the end of their\n"
            "\tsequence, e.g. {} ends at {} but is only {} bp.\n"
            "\tThis is the signature of an annotation built against a\n"
            "\tdifferent version of the assembly."
            .format(len(over), sid, end, ln)
        )

    print("  validated: {} sequence(s), {} feature(s), coordinates consistent"
          .format(len(seqs), len(feats)))


# Tool execution: every tool runs from a pinned image resolved through
# config/containers.json, never from `module load`.

# Directories a build step reads or writes inside an image, set by build() from
# --output so no cluster path is hardcoded here.
_CONTAINER_BINDS = []

# Platform the images resolve against, set by build() from --platform. Module
# state, so every _run_cmd caller does not have to pass it through.
_PLATFORM = containers.DEFAULT_PLATFORM


def _set_platform(platform):
    global _PLATFORM
    _PLATFORM = platform or containers.DEFAULT_PLATFORM


def _load_images():
    """Image map for the platform this build is running on."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        return containers.resolve_images(here, _PLATFORM)
    except Exception as exc:
        fatal("Cannot resolve container images for platform {}: {}"
              .format(_PLATFORM, exc))


def _singularity_prefix(image):
    """argv prefix that runs a command inside `image`."""
    here  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    roots = list(containers.image_roots(here, _PLATFORM).values())
    binds = ",".join(dict.fromkeys(d for d in _CONTAINER_BINDS + roots if os.path.isdir(d)))
    cmd = ["singularity", "exec"]
    if binds:
        cmd += ["--bind", binds]
    return cmd + [image]


def _run_cmd(cmd, label, log_file, image=None, check=True):
    """Run a command, inside `image` when one is given, appending output to log_file.
    With check=False a failure is returned rather than fatal."""
    import shlex
    print("    \u2192 {}".format(label))
    argv = _singularity_prefix(image) + list(cmd) if image else list(cmd)
    cmd_str = " ".join(shlex.quote(str(c)) for c in argv)
    with open(log_file, "a") as lf:
        lf.write("\n$ {}\n".format(cmd_str))
        lf.flush()
        ret = subprocess.call(
            cmd_str, shell=True, executable="/bin/bash",
            stdout=lf, stderr=subprocess.STDOUT,
        )
    if ret != 0 and check:
        if subprocess.call("command -v singularity", shell=True,
                           executable="/bin/bash",
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL) != 0:
            fatal(
                "singularity not found.\n"
                "  'viralrecon build' runs every tool from a container:\n"
                "    module load singularity"
            )
        if image and not os.path.isfile(image):
            fatal("Container image missing: {}".format(image))
        fatal("{} failed (exit {}).  See log: {}".format(label, ret, log_file))
    return ret


def _build_index(canonical_name, genome_dir):
    """Build every index the pipeline needs inside genome_dir (faidx, dict, bowtie2,
    snpEff). Idempotent: each step is skipped if its output already exists."""
    images = _load_images()
    fasta_path = os.path.join(genome_dir, "{}.fa".format(canonical_name))
    fai_path   = fasta_path + ".fai"
    dict_path  = os.path.join(genome_dir, "{}.dict".format(canonical_name))
    bt2_done   = os.path.join(genome_dir, "{}.1.bt2".format(canonical_name))
    snpeff_cfg = os.path.join(genome_dir, "snpEff.config")
    log        = os.path.join(genome_dir, "build_index.log")

    # Detect annotation file
    gff_path = os.path.join(genome_dir, "genes.gff")
    gtf_path = os.path.join(genome_dir, "genes.gtf")
    if os.path.isfile(gff_path):
        ann_src, gtf_ext, gtf_flag = gff_path, ".gff", "-gff3"
    elif os.path.isfile(gtf_path):
        ann_src, gtf_ext, gtf_flag = gtf_path, ".gtf", "-gtf22"
    else:
        fatal("No annotation file found in {}.\n"
              "  Expected: genes.gff or genes.gtf".format(genome_dir))

    print("\n  Indexing '{}' …".format(canonical_name))
    print("  Log → {}".format(log))

    # 1. samtools faidx
    if os.path.isfile(fai_path):
        print("    ✓ .fa.fai exists, skip")
    else:
        _run_cmd(["samtools", "faidx", fasta_path], "samtools faidx", log,
                 image=images["samtools"])

    # 2. samtools dict
    if os.path.isfile(dict_path):
        print("    ✓ .dict exists, skip")
    else:
        _run_cmd(
            ["samtools", "dict", fasta_path, "-o", dict_path],
            "samtools dict", log, image=images["samtools"],
        )

    # 3. bowtie2-build
    if os.path.isfile(bt2_done):
        print("    ✓ .1.bt2 exists, skip")
    else:
        _run_cmd(
            [
                "bowtie2-build",
                "--threads", "4",
                "-f", fasta_path,
                os.path.join(genome_dir, canonical_name),
            ],
            "bowtie2-build",
            log, image=images["bowtie2"],
        )

    # 4. snpEff build
    if os.path.isfile(snpeff_cfg):
        print("    ✓ snpEff.config exists, skip")
    else:
        # snpEff resolves data.dir/<genome_name>/genes.gff, so data.dir must be
        # the parent of genome_dir.
        data_dir = os.path.dirname(genome_dir)
        with open(snpeff_cfg, "w") as fh:
            fh.write("data.dir = {}\n".format(data_dir))
            fh.write("{}.genome: {}\n".format(canonical_name, canonical_name))

        # Stage copies that snpEff expects (sequences.fa + genes.{ext})
        seq_dst = os.path.join(genome_dir, "sequences.fa")
        if not os.path.isfile(seq_dst):
            shutil.copy(fasta_path, seq_dst)
        ann_dst = os.path.join(genome_dir, "genes{}".format(gtf_ext))
        if ann_src != ann_dst and not os.path.isfile(ann_dst):
            shutil.copy(ann_src, ann_dst)

        _run_cmd(
            [
                "snpEff", "build", gtf_flag,
                "-noCheckCds", "-noCheckProtein",
                "-dataDir", data_dir,
                "-config", snpeff_cfg,
                "-v", canonical_name,
            ],
            "snpEff build",
            log, image=images["snpeff"],
        )

    print("  ✓ Indexing complete")


# genome.json management

def _fetch_nextclade_dataset(dataset_name, genome_dir, log):
    """Download a Nextclade dataset into the reference directory, where genome.json
    points to it per target. Needs network, so set http_proxy on compute nodes."""
    images = _load_images()
    out_dir = os.path.join(genome_dir, "nextclade")
    if os.path.isdir(out_dir) and os.path.isfile(os.path.join(out_dir, "pathogen.json")):
        print("    \u2192 nextclade dataset already present")
        return out_dir

    print("  Fetching Nextclade dataset '{}' \u2026".format(dataset_name))
    _run_cmd(
        ["nextclade", "dataset", "get",
         "--name", dataset_name, "--output-dir", out_dir],
        "nextclade dataset get", log, image=images["nextclade"],
    )
    if not os.path.isfile(os.path.join(out_dir, "pathogen.json")):
        fatal(
            "\n\tNextclade dataset '{}' did not download.\n"
            "\tCheck the name against `nextclade dataset list`, and that this node\n"
            "\treaches the internet (Biowulf: export https_proxy=http://dtn20-e0:3128)"
            .format(dataset_name)
        )
    return out_dir


def _match_nextclade_dataset(fasta, genome_dir, log):
    """The one dataset `nextclade sort` matches every FASTA record to, else None.
    No match, a split match across segments, or no network all skip lineage."""
    images = _load_images()
    tsv = os.path.join(genome_dir, "nextclade_sort.tsv")
    print("  Matching a Nextclade dataset \u2026")
    ret = _run_cmd(
        ["nextclade", "sort", "--output-results-tsv", tsv, fasta],
        "nextclade sort", log, image=images["nextclade"], check=False,
    )
    if ret != 0 or not os.path.isfile(tsv):
        print("    nextclade sort failed (no internet?); lineage skipped. See {}"
              .format(log))
        return None
    with open(tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    matched = {r["seqName"]: r["dataset"] for r in rows if r.get("dataset")}
    names = sorted(set(matched.values()))
    if not names:
        print("    no Nextclade dataset matches this reference; lineage skipped")
        return None
    if len(names) > 1 or len(matched) < len({r["seqName"] for r in rows}):
        print("    records do not all match one dataset ({}); lineage skipped"
              .format(", ".join(names)))
        return None
    print("    matched {}".format(names[0]))
    return names[0]


def _choose_nextclade_dataset(no_nextclade, fasta, genome_dir, log):
    """The automatically matched dataset, or None when --no-nextclade opts out."""
    if no_nextclade:
        print("  Nextclade: skipped (--no-nextclade)")
        return None
    return _match_nextclade_dataset(fasta, genome_dir, log)


def _existing_annotation(genome_dir):
    """Whichever annotation form this reference was built with."""
    for ext in (".gff", ".gtf"):
        p = os.path.join(genome_dir, "genes" + ext)
        if os.path.isfile(p):
            return p
    return os.path.join(genome_dir, "genes.gff")


def _reference_artifacts(genome_dir, name):
    """Every file the pipeline links from a built reference. Keep in step with the
    ln -sf list in custom_virmapDB."""
    required = [
        os.path.join(genome_dir, "{}.fa".format(name)),
        os.path.join(genome_dir, "{}.fa.fai".format(name)),
        os.path.join(genome_dir, "{}.dict".format(name)),
        os.path.join(genome_dir, "snpEff.config"),
        os.path.join(genome_dir, "sequences.fa"),
        os.path.join(genome_dir, "snpEffectPredictor.bin"),
    ]
    # bowtie2 writes either a standard or a large index, never both
    bt2 = [os.path.join(genome_dir, "{}.{}.bt2".format(name, p))
           for p in ("1", "2", "3", "4", "rev.1", "rev.2")]
    bt2l = [p + "l" for p in bt2]
    # annotation is .gff or .gtf depending on what the source provided
    ann = [os.path.join(genome_dir, "genes.gff"),
           os.path.join(genome_dir, "genes.gtf")]
    return required, bt2, bt2l, ann


def _missing_artifacts(genome_dir, name):
    """Which required files are absent or empty. Empty list means complete."""
    def bad(p):
        return not os.path.isfile(p) or os.path.getsize(p) == 0

    required, bt2, bt2l, ann = _reference_artifacts(genome_dir, name)
    missing = [p for p in required if bad(p)]
    if any(bad(p) for p in bt2) and any(bad(p) for p in bt2l):
        missing.append(os.path.join(genome_dir, "{}.*.bt2".format(name)))
    if all(bad(p) for p in ann):
        missing.append(os.path.join(genome_dir, "genes.gff|gtf"))
    return missing


def _update_genome_json(genome_json_path, canonical_name, platforms,
                        fasta_path, gff_path, taxid=None, notes=None,
                        nextclade_dataset=None):
    """Add or update one target in genome.json atomically, keeping every other entry.
    An existing note is kept unless --notes is passed (--notes "" clears it)."""
    if os.path.isfile(genome_json_path):
        try:
            with open(genome_json_path) as fh:
                data = json.load(fh)
        except ValueError as exc:
            fatal(
                "\n\tgenome.json is not valid JSON and was left untouched:\n"
                "\t  {}\n"
                "\t  {}\n"
                "\tIt is an accumulating registry of every target built into\n"
                "\tthis directory, so it is not overwritten automatically.\n"
                "\tFix or remove it, then re-run the build."
                .format(genome_json_path, exc)
            )
    else:
        data = {"references": {"contamination": {}, "kraken2": {}, "target": {}}}

    data.setdefault("references", {}).setdefault("target", {})
    existing = set()
    for plat in platforms:
        tgts = data["references"]["target"].setdefault(plat, {})
        if canonical_name in tgts:
            existing.add(plat)
        prior = tgts.get(canonical_name, {})
        entry = {"fasta": fasta_path, "gtf": gff_path}
        # The entry is rebuilt from scratch, so carry over any field not passed in.
        if taxid:
            entry["taxid"] = str(taxid)
        elif prior.get("taxid"):
            entry["taxid"] = prior["taxid"]
        if notes is None:
            if prior.get("notes"):
                entry["notes"] = prior["notes"]
        elif notes != "":
            entry["notes"] = notes
        # Carried forward like notes, so a failed match cannot drop a dataset;
        # "" (from --no-nextclade) clears it.
        if nextclade_dataset is None:
            if prior.get("nextclade_dataset"):
                entry["nextclade_dataset"] = prior["nextclade_dataset"]
        elif nextclade_dataset != "":
            entry["nextclade_dataset"] = nextclade_dataset
        tgts[canonical_name] = entry
    if existing:
        print("  updating existing entry for '{}'".format(canonical_name))

    tmp = genome_json_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=4)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, genome_json_path)

    print("  genome.json → {}".format(genome_json_path))
    print("  Registered as '{}'".format(canonical_name))
    _note = data["references"]["target"][platforms[0]][canonical_name].get("notes")
    if _note:
        print("  note          : {}".format(_note))


# Entry point

def build(sub_args, repo_path):
    """Entry point for `viralrecon build`: fetch FASTA and GFF, build the indices,
    and register the target in --output/genome.json."""
    virus      = sub_args.virus
    accession  = getattr(sub_args, "accession", None)
    outdir     = os.path.abspath(sub_args.output)
    platform   = getattr(sub_args, "platform", None)
    local_fasta      = getattr(sub_args, "fasta",       None)
    local_annotation = getattr(sub_args, "annotation",  None)
    no_index         = getattr(sub_args, "no_index",    False)

    local_pair = bool(local_fasta and local_annotation)

    # ── Source selection: accession XOR local files ─────────────────────────
    # Both say where the reference comes from, so accepting both is ambiguous.
    if accession and local_pair:
        fatal(
            "\n\t--accession and --fasta/--annotation are mutually exclusive.\n"
            "\tEither download a reference from NCBI:\n"
            "\t  viralrecon build --virus SARS --accession NC_045512.2 --output DIR\n"
            "\tor supply your own pair:\n"
            "\t  viralrecon build --name MYVIRUS --fasta my.fa --annotation my.gff --output DIR\n"
            "\tMixing them is refused rather than guessed at because the pair\n"
            "\thas to be coordinate-consistent, and pairing a downloaded FASTA\n"
            "\twith an unrelated annotation cannot be detected after the fact."
        )
    if bool(local_fasta) != bool(local_annotation):
        fatal(
            "\n\t--fasta and --annotation must be given together.\n"
            "\tAn annotation is only meaningful against the sequence it was\n"
            "\tbuilt from; supplying one without the other would silently pair\n"
            "\tit with a different genome."
        )
    if not accession and not local_pair and not virus:
        fatal(
            "\n\tNothing to build. Supply either --accession (optionally with\n"
            "\t--virus), or --fasta with --annotation and a --name."
        )

    # Resolve accession from preset if not provided
    if not accession and not local_pair:
        preset_acc = VIRUS_PRESETS.get(virus.upper())
        if not preset_acc:
            fatal(
                "\n\tNo --accession given and '{}' is not a built-in preset.\n"
                "\tBuilt-in presets: {}\n"
                "\tSupply --accession <NCBI_accession> explicitly."
                .format(virus, ", ".join(sorted(VIRUS_PRESETS)))
            )
        accession = preset_acc
        print("Using preset accession for {}: {}".format(virus.upper(), accession))

    # Local pairs are named by the user: there is no accession to derive a
    # name from, and the target name is what appears in every output path.
    if local_pair:
        canonical_name = (getattr(sub_args, "name", None)
                          or (virus.upper().replace(" ", "_") if virus else None))
        if not canonical_name:
            fatal(
                "\n\tA --name (or --virus) is required with --fasta/--annotation.\n"
                "\tThe target name appears in every output path and in\n"
                "\tgenome.json, so it cannot be inferred from a filename."
            )
    else:
        canonical_name = _canonical_name(virus, accession)
    genome_dir     = os.path.join(outdir, canonical_name)
    # Register only the platform being built on, since genome.json paths are
    # absolute on this filesystem. Defaults to BIOWULF, matching `viralrecon run`.
    platforms      = [platform] if platform else ["BIOWULF"]
    _set_platform(platforms[0])
    _CONTAINER_BINDS[:] = [outdir, os.path.realpath(outdir)]

    genome_json = os.path.join(outdir, "genome.json")
    force = getattr(sub_args, "force", False)

    no_nc = getattr(sub_args, "no_nextclade", False)

    # ── Already built? Verify, then skip ────────────────────────────────────
    # A registry entry is not proof the files exist, so check both before skipping.
    registered = False
    if os.path.isfile(genome_json):
        try:
            with open(genome_json) as fh:
                _reg = json.load(fh).get("references", {}).get("target", {})
            registered = any(canonical_name in _reg.get(pl, {}) for pl in platforms)
        except ValueError:
            pass          # reported properly by _update_genome_json later

    missing = _missing_artifacts(genome_dir, canonical_name)

    if registered and not missing and not force:
        # Complete files do not mean a complete entry: backfill a missing taxid
        # without redoing the file work.
        try:
            with open(genome_json) as fh:
                _reg2 = json.load(fh)["references"]["target"]
            _entry = {}
            for pl in platforms:
                _entry = _reg2.get(pl, {}).get(canonical_name, {}) or _entry
        except Exception:
            _entry = {"taxid": True, "nextclade_dataset": True}

        _taxid_arg = getattr(sub_args, "taxid", None)
        taxid = None
        if not _entry.get("taxid"):
            # A local pair has no accession to look the taxid up from, so it
            # gets one only when --taxid says so.
            taxid = _taxid_arg or (None if local_pair else _efetch_taxid(accession))
            if taxid:
                print("\n  backfilling missing taxid {} for '{}'"
                      .format(taxid, canonical_name))

        # Likewise add a missing Nextclade dataset, which complete files cannot
        # imply and later builds would otherwise never add.
        nc_path = None
        if no_nc and _entry.get("nextclade_dataset"):
            print("\n  --no-nextclade: removing the nextclade dataset from '{}'"
                  .format(canonical_name))
            nc_path = ""
        elif not _entry.get("nextclade_dataset"):
            _log = os.path.join(genome_dir, "build_index.log")
            nc_dataset = _choose_nextclade_dataset(
                no_nc, os.path.join(genome_dir, "{}.fa".format(canonical_name)),
                genome_dir, _log)
            if nc_dataset:
                print("\n  backfilling nextclade dataset '{}' for '{}'"
                      .format(nc_dataset, canonical_name))
                nc_path = _fetch_nextclade_dataset(nc_dataset, genome_dir, _log)

        if taxid or nc_path is not None:
            _update_genome_json(
                genome_json, canonical_name, platforms,
                os.path.join(genome_dir, "{}.fa".format(canonical_name)),
                _existing_annotation(genome_dir),
                taxid=taxid,
                nextclade_dataset=nc_path,
            )
        print("\n✓ Genome '{}' is already built and complete, skipping."
              .format(canonical_name))
        print("  Reference dir : {}".format(genome_dir))
        print("  genome.json   : {}".format(genome_json))
        print("  Re-build anyway with --force.")
        return

    if registered and missing:
        print("\nGenome '{}' is registered but incomplete; rebuilding what is "
              "missing:".format(canonical_name))
        for m in missing[:8]:
            print("    missing: {}".format(os.path.basename(m)))
    elif force and registered:
        print("\n--force given; rebuilding '{}' from scratch."
              .format(canonical_name))

    print("\nBuilding genome '{}' …".format(canonical_name))
    os.makedirs(genome_dir, exist_ok=True)

    # Obtain FASTA + GFF
    if local_pair:
        _validate_local(local_fasta, local_annotation)
        fasta_dst = os.path.join(genome_dir, "{}.fa".format(canonical_name))
        gff_dst   = os.path.join(genome_dir, "genes.gff")
        if not os.path.isfile(fasta_dst):
            shutil.copy(os.path.abspath(local_fasta), fasta_dst)
        if not os.path.isfile(gff_dst):
            shutil.copy(os.path.abspath(local_annotation), gff_dst)
        fasta_path = os.path.abspath(fasta_dst)
        gff_path   = os.path.abspath(gff_dst)
        print("  FASTA → {}".format(fasta_path))
        print("  GFF   → {}".format(gff_path))
    else:
        fasta_path, gff_path = _download_reference(accession, genome_dir, canonical_name)

    # Build indices (idempotent)
    if not no_index:
        _build_index(canonical_name, genome_dir)

    nc_path = "" if no_nc else None
    _log = os.path.join(genome_dir, "build_index.log")
    nc_dataset = _choose_nextclade_dataset(no_nc, fasta_path, genome_dir, _log)
    if nc_dataset:
        nc_path = _fetch_nextclade_dataset(nc_dataset, genome_dir, _log)

    # Write genome.json inside the reference output directory
    taxid = getattr(sub_args, "taxid", None)
    if not taxid and not local_pair:
        taxid = _efetch_taxid(accession)
        if taxid:
            print("  NCBI taxid   : {}".format(taxid))
        else:
            print("  NCBI taxid   : not resolved; the composition profile will "
                  "not break this target out separately")
    _update_genome_json(
        genome_json, canonical_name, platforms,
        os.path.abspath(fasta_path),
        os.path.abspath(gff_path),
        taxid=taxid,
        notes=getattr(sub_args, "notes", None),
        nextclade_dataset=nc_path,
    )

    print("\n✓ Genome '{}' ready.".format(canonical_name))
    print("  Reference dir : {}".format(genome_dir))
    print("  genome.json   : {}".format(genome_json))
    print("\nTo run the pipeline against this (and other built) targets:")
    print("  viralrecon run --input <FASTQ...> --output <DIR> \\")
    print("                 --genome {}".format(genome_json))
