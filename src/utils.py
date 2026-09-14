#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Shared utility functions for the viralrecon CLI."""

import os
import re
import subprocess
import sys


class Colors:
    """ANSI escape sequences for styling terminal output."""
    end        = '\33[0m'
    bold       = '\33[1m'
    italic     = '\33[3m'
    url        = '\33[4m'
    red        = '\33[31m'
    green      = '\33[32m'
    yellow     = '\33[33m'
    cyan       = '\33[96m'
    white      = '\33[37m'
    bg_red     = '\33[41m'
    bg_black   = '\33[40m'


def err(*message, **kwargs):
    """Print to standard error."""
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Print to standard error and exit with code 1."""
    err(*message, **kwargs)
    sys.exit(1)


def exists(path):
    """Return True if path exists on the filesystem."""
    return os.path.exists(path)


def which(cmd, path=None):
    """Return True if executable is in $PATH."""
    if path is None:
        path = os.environ.get('PATH', '').split(os.pathsep)
    for prefix in path:
        candidate = os.path.join(prefix, cmd)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return True
    return False


def require(cmds, suggestions, path=None):
    """Enforce that each executable in cmds is in $PATH."""
    c = Colors()
    error = False
    for i, cmd in enumerate(cmds):
        if not which(cmd, path):
            error = True
            err(
                "\n{0}{1}Fatal: '{2}' is not in $PATH and is required at runtime!{3}"
                "\n  └── Possible fix: module load {4}".format(
                    c.bg_red, c.white, cmd, c.end, suggestions[i]
                )
            )
    if error:
        fatal()


def permissions(parser, path, *args, **kwargs):
    """Return the absolute path if it exists with the requested access, otherwise
    call parser.error()."""
    if not exists(path):
        parser.error("Path '{}' does not exist.".format(path))
    if not os.access(path, *args, **kwargs):
        parser.error("Path '{}' exists but cannot be accessed (permission denied).".format(path))
    return os.path.abspath(path)


def check_cache(parser, cache):
    """Validate a Singularity cache directory, creating it if absent. Errors if it is
    a file or owned by another user."""
    c = Colors()
    if not exists(cache):
        os.makedirs(cache)
    elif os.path.isfile(cache):
        parser.error(
            "\n\t{0}Fatal: --singularity-cache '{1}' already exists as a file.{2}"
            "\n\tPlease re-run with a different path.".format(c.red, cache, c.end)
        )
    elif os.path.isdir(cache):
        inner = os.path.join(cache, 'cache')
        if exists(inner) and os.stat(inner).st_uid != os.getuid():
            parser.error(
                "\n\t{0}Fatal: '{1}' is owned by a different user.{2}"
                "\n\tSingularity caches cannot be shared across users."
                "\n\tPlease re-run with a different --singularity-cache path.".format(
                    c.red, cache, c.end
                )
            )
    return cache


def check_snakemake_version():
    """Verify that snakemake < 8.0.0 is available."""
    try:
        raw = subprocess.check_output(
            ['snakemake', '--version'], stderr=subprocess.STDOUT
        ).strip().decode('utf-8')
    except Exception:
        fatal("Fatal: 'snakemake' is not in $PATH. Please load or install snakemake <= 7.32.4.")

    m = re.search(r'(\d+)\.(\d+)\.(\d+)', raw.split()[-1])
    if m and int(m.group(1)) >= 8:
        fatal(
            "Fatal: Snakemake version '{}' is not supported.\n"
            "Please use snakemake <= 7.32.4.".format(raw)
        )
    return raw


def git_commit_hash(repo_path):
    """Return the HEAD commit hash of the git repo at repo_path."""
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            stderr=subprocess.STDOUT,
            cwd=repo_path,
        ).strip().decode('utf-8')
    except Exception:
        return 'github_release'


def join_jsons(template_files):
    """Merge a list of JSON files into a single dict (last key wins)."""
    import json
    merged = {}
    for f in template_files:
        with open(f, 'r') as fh:
            merged.update(json.load(fh))
    return merged


def unpacked(nested_dict):
    """Recursively yield all leaf values of a nested dict."""
    for v in nested_dict.values():
        if isinstance(v, dict):
            yield from unpacked(v)
        else:
            yield v
