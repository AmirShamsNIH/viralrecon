#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Install sub-command: download reference bundles listed in config/install.json."""

import json
import os

try:
    from .utils import exists, fatal
except ImportError:
    from utils import exists, fatal


def read_install_config(repo_path):
    """Load config/install.json from the pipeline repo."""
    config_file = os.path.join(repo_path, 'config', 'install.json')
    if not exists(config_file):
        return {'install': {}}
    with open(config_file, 'r') as fh:
        return json.load(fh)


def has_downloads(targets):
    """Return True if any download targets are defined."""
    return bool(targets)


def download_target(sub_args, chunks):
    """Download tarball shards for one target (stub, implement as needed)."""
    print("Download target: {}".format(chunks))


def assemble_target(sub_args, chunks):
    """Concatenate shards and extract the tarball (stub, implement as needed)."""
    print("Assemble target: {}".format(chunks))


def install(sub_args, repo_path):
    """Entry point for `viralrecon install`: download data defined in config/install.json."""
    cfg = read_install_config(repo_path)
    targets = cfg.get('install', {})

    if not has_downloads(targets):
        print("Nothing to install: config/install.json has no download targets.")
        return

    outdir = sub_args.output
    if not exists(outdir):
        os.makedirs(outdir)

    for name, chunks in targets.items():
        print("Installing: {}".format(name))
        download_target(sub_args, chunks)
        assemble_target(sub_args, chunks)

    print("Install complete: {}".format(outdir))
