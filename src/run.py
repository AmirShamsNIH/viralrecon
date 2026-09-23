#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Pipeline initialisation, config assembly and launch logic."""

import json
import os
import re
import subprocess
import sys
from shutil import copytree, ignore_patterns, rmtree

try:
    from .utils import Colors, err, exists, fatal, git_commit_hash, join_jsons, unpacked, which
except ImportError:
    from utils import Colors, err, exists, fatal, git_commit_hash, join_jsons, unpacked, which

try:
    from . import containers
except ImportError:
    import containers

try:
    from . import version as __version__
except ImportError:
    import os as _os
    __version__ = open(
        _os.path.join(_os.path.dirname(__file__), '..', 'VERSION')
    ).read().strip()

# Input-file normalisation

ILLUMINA_R1 = '.R1.fastq.gz'
ILLUMINA_R2 = '.R2.fastq.gz'
ILLUMINA_MATE_RE = re.compile(r'[_\.]R[12][_\.]?\d*\.f(?:ast)?q\.gz')

ILLUMINA_RENAME = {
    r'\.R1\.f(ast)?q\.gz$':                  ILLUMINA_R1,
    r'\.R2\.f(ast)?q\.gz$':                  ILLUMINA_R2,
    r'_R1_\d+\.f(ast)?q\.gz$':               ILLUMINA_R1,  # _R1_001.fastq.gz
    r'_R2_\d+\.f(ast)?q\.gz$':               ILLUMINA_R2,  # _R2_001.fastq.gz
    r'\.R1\.(?P<lane>...).f(ast)?q\.gz$':    ILLUMINA_R1,
    r'\.R2\.(?P<lane>...).f(ast)?q\.gz$':    ILLUMINA_R2,
    r'_R1\.f(ast)?q\.gz$':                   ILLUMINA_R1,
    r'_R2\.f(ast)?q\.gz$':                   ILLUMINA_R2,
    r'_1\.f(ast)?q\.gz$':                    ILLUMINA_R1,
    r'_2\.f(ast)?q\.gz$':                    ILLUMINA_R2,
}


def _rename_fastq(filename):
    """Normalise an Illumina FastQ name to canonical .R1/.R2.fastq.gz form."""
    if filename.endswith(ILLUMINA_R1) or filename.endswith(ILLUMINA_R2):
        return filename
    for pattern, canonical in ILLUMINA_RENAME.items():
        if re.search(pattern, filename):
            return re.sub(pattern, canonical, filename)
    raise NameError(
        "\n\tFatal: Cannot normalise FastQ name '{}'.\n"
        "\tExpected extensions: _R1.fastq.gz / _R2.fastq.gz / _1.fastq.gz\n"
        "\tPlease rename the file and retry.".format(filename)
    )


def _strip_sample(filename):
    """Return the sample basename by stripping mate suffix."""
    name = os.path.basename(filename)
    return ILLUMINA_MATE_RE.split(name)[0]


def _verify_mates(ifiles):
    """Ensure every paired-end sample has both R1 and R2."""
    counts = {}
    for f in ifiles:
        name = os.path.basename(f)
        if name.endswith(ILLUMINA_R1) or name.endswith(ILLUMINA_R2):
            sample = _strip_sample(name)
            counts[sample] = counts.get(sample, 0) + 1
    missing = [s for s, n in counts.items() if n == 1]
    if missing:
        fatal(
            "\n\tFatal: Paired-end data detected but a mate is missing for:\n"
            "\t  {}\n"
            "\tEnsure each sample has matching R1 and R2 files.".format(missing)
        )


# Output-directory initialisation

def init(repo_path, output_path, links=None, refresh=False):
    """Create output_path, copy workflow, config and resources into it, and symlink
    inputs into inputs/. Returns the renamed symlink paths."""
    links = links or []
    required = ['workflow', 'resources', 'config']

    if not exists(output_path):
        os.makedirs(output_path)
    elif os.path.isfile(output_path):
        fatal(
            "\n\tFatal: --output '{}' already exists as a file.\n"
            "\tPlease choose a different output path.".format(output_path)
        )

    # Scratch for the run, on the same filesystem as its outputs rather than
    # node-local /lscratch, so nothing a job produced is stranded on a node.
    os.makedirs(os.path.join(output_path, 'tmp'), exist_ok=True)

    _copy_safe(repo_path, output_path, required, refresh)
    return _sym_safe(links, output_path)


def _copy_safe(source, target, resources, refresh=False):
    """Copy resource directories into target. config and resources are copied once to
    keep user edits; workflow is always refreshed so rule fixes reach existing runs."""
    always_refresh = {'workflow'}
    for resource in resources:
        src = os.path.join(source, resource)
        dst = os.path.join(target, resource)
        should_refresh = refresh or (resource in always_refresh)
        if exists(dst) and should_refresh:
            if os.path.isdir(dst):
                rmtree(dst)
            else:
                os.remove(dst)
        if not exists(dst) and exists(src):
            # This copy is what executes and is the run's audit trail, so keep
            # interpreter and editor droppings out of it.
            copytree(src, dst, ignore=ignore_patterns(
                '__pycache__', '*.pyc', '*.pyo', '.DS_Store', '._*'))


def _sym_safe(input_files, target, input_dirname='inputs'):
    """Create renamed symlinks for input FastQ files in target/inputs/."""
    input_dir = os.path.join(target, input_dirname)
    if not exists(input_dir):
        os.makedirs(input_dir)

    renamed = []
    collision = {}
    for f in input_files:
        new_name = _rename_fastq(os.path.basename(f))
        dst      = os.path.join(input_dir, new_name)
        renamed.append(dst)
        collision.setdefault(new_name, []).append(f)

    clashes = {n: s for n, s in collision.items() if len(s) > 1}
    if clashes:
        fatal(
            "\n\tFatal: Input filename collision: two files would produce the "
            "same symlink:\n\t  {}".format(clashes)
        )

    for src, dst in zip(input_files, renamed):
        if not os.path.lexists(dst):
            os.symlink(os.path.abspath(os.path.realpath(src)), dst)

    return renamed


# Config assembly

def setup(sub_args, ifiles, repo_path, output_path):
    """Build the merged pipeline config from template JSON files + user inputs."""
    config_dir = os.path.join(output_path, 'config')
    template_files = [
        os.path.join(config_dir, 'config.json'),
        os.path.join(config_dir, 'containers.json'),
        os.path.join(config_dir, 'genome.json'),
    ]
    config = join_jsons(template_files)
    config['project'] = {}

    # User / environment metadata
    home = os.path.expanduser('~')
    config['project']['userhome'] = home
    config['project']['username'] = os.path.split(home)[-1]
    config['project']['version']  = __version__
    config['project']['git_commit_hash'] = git_commit_hash(repo_path)
    config['project']['pipeline_path']   = repo_path
    config['project']['workpath']        = os.path.abspath(output_path)

    # Input samples
    samples = []
    for f in ifiles:
        s = _strip_sample(f)
        if s not in samples:
            samples.append(s)
    config['samples'] = samples

    # Paired-end check
    has_r2 = any(os.path.basename(f).endswith(ILLUMINA_R2) for f in ifiles)
    if has_r2:
        _verify_mates(ifiles)
    config['project']['paired'] = has_r2

    # Raw-data bind paths for Singularity
    rawdata_dirs = list({
        os.path.dirname(os.path.abspath(os.path.realpath(f)))
        for f in (sub_args.input or [])
    })
    config['project']['datapath'] = ','.join(rawdata_dirs)

    # Target selection: driven by the genome.json from 'viralrecon build'
    platform = getattr(sub_args, 'platform', 'BIOWULF')

    genome_json_path = getattr(sub_args, 'genome', None)
    if genome_json_path:
        genome_json_path = os.path.abspath(genome_json_path)
        if not exists(genome_json_path):
            fatal(
                "\n\tFatal: genome.json not found: {}\n"
                "\tRun 'viralrecon build --virus X --accession Y --output <DIR>' first,\n"
                "\tthen pass the resulting genome.json via --genome."
                .format(genome_json_path)
            )
        with open(genome_json_path) as fh:
            genome_data = json.load(fh)

        build_refs = (
            genome_data
            .get("references", {})
            .get("target", {})
            .get(platform, {})
        )
        if not build_refs:
            fatal(
                "\n\tFatal: No targets for platform '{}' in: {}\n"
                "\tRe-run 'viralrecon build --platform {}' to register targets."
                .format(platform, genome_json_path, platform)
            )

        # Merge build references into config, overriding any stale genome.json copy
        config.setdefault('references', {}).setdefault('target', {})[platform] = build_refs

        # Record the genome directory for Singularity bind-path resolution
        config['project']['genomepath'] = os.path.dirname(genome_json_path)

        # Honour --targets subset if given, otherwise use everything in genome.json
        requested = getattr(sub_args, 'targets', None) or []
        if requested:
            unknown = [t for t in requested if t not in build_refs]
            if unknown:
                fatal(
                    "\n\tFatal: Unknown target(s): {}\n"
                    "\tTargets in {}: {}"
                    .format(unknown, genome_json_path, sorted(build_refs))
                )
            config['targets'] = requested
        else:
            config['targets'] = list(build_refs.keys())

        # A note on a reference exists to be read before a run, not discovered
        # in genome.json afterwards, so echo it for every selected target.
        for _t in config['targets']:
            _n = (build_refs.get(_t) or {}).get('notes')
            if _n:
                print("  note [{}]: {}".format(_t, _n))
    else:
        # Fallback: legacy mode, look up targets in the pipeline's own genome.json
        config['project']['genomepath'] = ''
        all_targets = list(
            config.get('references', {}).get('target', {}).get(platform, {}).keys()
        )
        requested = getattr(sub_args, 'targets', None) or []
        if requested:
            unknown = [t for t in requested if t not in all_targets]
            if unknown:
                fatal(
                    "\n\tFatal: Unknown target(s): {}\n"
                    "\tAvailable for platform {}: {}\n"
                    "\tRegister new ones with: viralrecon build --virus X --accession Y --output <DIR>"
                    .format(unknown, platform, all_targets)
                )
            config['targets'] = requested
        else:
            config['targets'] = all_targets

    # CLI options
    config.setdefault('options', {})
    for opt, val in vars(sub_args).items():
        if opt == 'func':
            continue
        config['options'][opt] = str(val) if not isinstance(val, (list, dict)) else val

    # Resolve image and platform-dependent database paths once, so the Snakemake
    # config carries only real paths and bind paths follow from them.
    _container_data = {k: config[k] for k in ('roots', 'images') if k in config}
    config['images'] = containers.resolve_images(
        repo_path, platform, data=_container_data)
    config.pop('roots', None)
    _apply_platform_paths(config, platform)
    _apply_platform_partition(output_path, platform)

    # Bind paths for Singularity (kept for container-mode runs)
    config['bindpaths'] = _resolve_bind_paths(sub_args, config)

    # ── Slim the config before writing ────────────────────────────────────────
    # Keep only fields the workflow reads, so the file stays auditable by hand.

    # options: only what the Snakefile reads
    _WORKFLOW_OPTS = {"output", "platform"}
    config["options"] = {k: v for k, v in config.get("options", {}).items()
                         if k in _WORKFLOW_OPTS}

    # project: only paired is consumed by the workflow
    config["project"] = {"paired": config.get("project", {}).get("paired", False)}

    # Remove the pipeline-stage metadata block (rule names, execution order)
    config.pop("pipeline", None)

    # Drop references, tools, and paths for the non-active platform
    for section in ("target", "contamination", "kraken2"):
        refs = config.get("references", {}).get(section, {})
        for plat in list(refs.keys()):
            if plat != platform:
                del refs[plat]

    # tools and paths are keyed <name> -> <platform> -> value. Trim only mappings,
    # so commentary keys such as "_comment" lists are left alone.
    for _section in ("tools", "paths"):
        for _name, _map in list(config.get(_section, {}).items()):
            if _name.startswith("_") or not isinstance(_map, dict):
                continue
            for plat in list(_map.keys()):
                if plat != platform and not plat.startswith("_"):
                    del _map[plat]

    return config


def _apply_platform_paths(config, platform):
    """Write the active platform's value of each config['paths'] entry into the
    parameter that reads it. A platform with no entry is fatal, never defaulted."""
    for name, spec in (config.get('paths') or {}).items():
        if name.startswith('_') or not isinstance(spec, dict):
            continue
        stage = spec.get('_parameter_stage')
        if not stage:
            continue
        if platform not in spec:
            fatal(
                "\n\tFatal: config.json paths.{} has no entry for platform {}."
                "\n\tAdd one, or run with a platform that is listed: {}"
                .format(name, platform,
                        sorted(k for k in spec if not k.startswith('_')))
            )
        config.setdefault('parameters', {}).setdefault(stage, {})[name] = spec[platform]


def _apply_platform_partition(output_path, platform):
    """Set __default__.partition in the run directory's cluster.json to this platform's
    queue (norm on Biowulf, all on BigSky and Skyline); the repository file is left untouched."""
    path = os.path.join(output_path, 'config', 'cluster.json')
    if not exists(path):
        return
    with open(path) as fh:
        cluster = json.load(fh)
    queue = (cluster.get('__partition__') or {}).get(platform)
    if not queue:
        return
    cluster.setdefault('__default__', {})['partition'] = queue
    with open(path, 'w') as fh:
        json.dump(cluster, fh, indent=4)


def _resolve_bind_paths(sub_args, config):
    """Collect all filesystem paths that Singularity needs to bind."""
    paths = set()
    workpath = config['project']['workpath']
    paths.add(workpath)

    for raw in config['project']['datapath'].split(','):
        if raw:
            paths.add(os.path.realpath(raw))

    # Bind the pre-built reference directory so Singularity can see it
    genomepath = config['project'].get('genomepath', '')
    if genomepath:
        paths.add(os.path.realpath(genomepath))

    for leaf in unpacked(config):
        if isinstance(leaf, str) and exists(leaf):
            p = os.path.dirname(leaf) if os.path.isfile(leaf) else leaf
            paths.add(p)

    return list(paths - {os.sep})


def build_config(sub_args, pl_home):
    """Top-level: initialise output dir, build config, return it."""
    ifiles = init(
        repo_path=pl_home,
        output_path=sub_args.output,
        links=sub_args.input or [],
        refresh=getattr(sub_args, 'overwrite_pipeline_template', False),
    )
    return setup(sub_args, ifiles, pl_home, sub_args.output)


def save_config(config, output_path):
    """Write the merged config to output_path/config.json."""
    dst = os.path.join(output_path, 'config.json')
    with open(dst, 'w') as fh:
        json.dump(config, fh, indent=4, sort_keys=True)


# Pipeline execution

def dryrun(outdir, config='config.json',
           snakefile=os.path.join('workflow', 'Snakefile')):
    """Dry-run snakemake to surface config/rule errors before a real run."""
    try:
        return subprocess.check_output([
            'snakemake', '-npr',
            '-s', snakefile,
            '--rerun-incomplete',
            '--cores', '256',
            '--configfile={}'.format(config),
        ], cwd=outdir, stderr=subprocess.STDOUT)
    except OSError as e:
        if e.errno == 2 and not which('snakemake'):
            fatal("Fatal: 'snakemake' not found in $PATH.")
        raise
    except subprocess.CalledProcessError as e:
        print(e.output.decode('utf-8'))
        raise


def launch_pipeline(sub_args, bindpaths, pl_home, pl_name):
    """Submit or run the Snakemake master job; returns the Popen object."""
    logdir = os.path.join(sub_args.output, 'logfiles')
    if not exists(logdir):
        os.makedirs(logdir)

    logname = 'snakemake.log' if sub_args.mode == 'local' else 'master.log'
    log = os.path.join(logdir, logname)

    with open(log, 'w') as logfh:
        mjob = _runner(
            mode=sub_args.mode,
            outdir=sub_args.output,
            alt_cache=getattr(sub_args, 'singularity_cache', None),
            threads=int(getattr(sub_args, 'threads', 2)),
            jobname=getattr(sub_args, 'job_name', 'viralrecon'),
            submission_script=os.path.join(pl_home, 'src', 'run.sh'),
            logger=logfh,
            additional_bind_paths=','.join(bindpaths),
            tmp_dir=getattr(sub_args, 'tmp_dir', None) or os.path.join(sub_args.output, 'tmp'),
        )
        if not getattr(sub_args, 'silent', False):
            print("\nRunning {} pipeline in '{}' mode …".format(pl_name, sub_args.mode))
        mjob.wait()

    return mjob


def report_outcome(sub_args, mjob, pl_name):
    """Print success/failure message after the pipeline completes."""
    if sub_args.mode == 'local':
        if mjob.returncode == 0:
            print('{} pipeline completed successfully.'.format(pl_name))
        else:
            fatal('{} pipeline failed. See logfiles/snakemake.log for details.'.format(pl_name))
    elif sub_args.mode == 'slurm':
        jobid_file = os.path.join(sub_args.output, 'logfiles', 'mjobid.log')
        with open(jobid_file) as fh:
            jobid = fh.read().strip()
        if mjob.returncode == 0:
            if not getattr(sub_args, 'silent', False):
                print('Master job submitted: ', end='')
            print(jobid)
        else:
            fatal('Failed to submit the master job.')


def _runner(mode, outdir, alt_cache, logger, additional_bind_paths='',
            threads=2, jobname='viralrecon', submission_script='src/run.sh',
            tmp_dir=None):
    """Internal: launch snakemake locally or via SLURM."""
    outdir = os.path.abspath(outdir)
    env = dict(os.environ)
    cache = os.path.join(outdir, '.singularity')
    env['SINGULARITY_CACHEDIR'] = alt_cache if alt_cache else cache
    if alt_cache:
        cache = alt_cache

    bindpaths = outdir
    tmp = os.path.dirname(tmp_dir.rstrip('/'))
    if tmp == os.sep:
        tmp = tmp_dir.rstrip('/')
    if tmp not in bindpaths.split(','):
        bindpaths = bindpaths + ',' + tmp
    if additional_bind_paths:
        bindpaths = additional_bind_paths + ',' + bindpaths

    if not exists(cache):
        os.makedirs(cache)

    if mode == 'local':
        return subprocess.Popen([
            'snakemake', '-pr', '--rerun-incomplete',
            '--use-singularity',
            '--singularity-args', "'-B {}'".format(bindpaths),
            '--cores', str(threads),
            '--configfile=config.json',
        ], cwd=outdir, stderr=subprocess.STDOUT, stdout=logger, env=env)

    elif mode == 'slurm':
        return subprocess.Popen([
            str(submission_script), mode,
            '-j', jobname,
            '-b', bindpaths,
            '-o', outdir,
            '-c', cache,
            '-t', tmp_dir,
        ], cwd=outdir, stderr=subprocess.STDOUT, stdout=logger, env=env)

    else:
        fatal("Unknown execution mode: '{}'.".format(mode))


# Entry-point functions called by the CLI

def run(sub_args, repo_path):
    """Entry point for `viralrecon run`: build the config, optionally dry-run, then launch."""
    pl_name = 'viralrecon'

    # 1. Initialise output directory and assemble config
    config = build_config(sub_args, repo_path)
    save_config(config, sub_args.output)

    # 2. Dry-run only
    if getattr(sub_args, 'dry_run', False):
        out = dryrun(sub_args.output)
        if out:
            print(out.decode('utf-8'))
        return

    # 3. Real run
    bindpaths = config.get('bindpaths', [])
    mjob = launch_pipeline(sub_args, bindpaths, repo_path, pl_name)
    report_outcome(sub_args, mjob, pl_name)


def unlock(sub_args, repo_path):
    """Entry point for `viralrecon unlock`: remove a stale Snakemake lock."""
    outdir = sub_args.output
    try:
        subprocess.check_call(
            ['snakemake', '--unlock',
             '-s', os.path.join('workflow', 'Snakefile'),
             '--configfile=config.json'],
            cwd=outdir,
        )
        print('Unlocked: {}'.format(outdir))
    except subprocess.CalledProcessError:
        fatal("Could not unlock '{}'. Is snakemake in $PATH?".format(outdir))
