#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Cache sub-command: pull Singularity images from Docker Hub locally."""

import json
import os
import subprocess

try:
    from .utils import err, exists, fatal
except ImportError:
    from utils import err, exists, fatal


def prepare_cache(sif_cache, pipeline_name='viralrecon'):
    """Create the SIF cache directory if it does not exist."""
    if not exists(sif_cache):
        os.makedirs(sif_cache)
    print("Singularity SIF cache: {}".format(sif_cache))


def missing_images(images_config, sif_cache):
    """Return URIs for images not yet present in the local SIF cache."""
    with open(images_config, 'r') as fh:
        data = json.load(fh)
    missing = []
    for name, uri in data.get('images', {}).items():
        sif_name = uri.replace('docker://', '').replace('/', '_').replace(':', '_') + '.sif'
        sif_path = os.path.join(sif_cache, sif_name)
        if not exists(sif_path):
            missing.append((name, uri, sif_path))
    return missing


def pull_images(sif_cache, uris_to_pull):
    """Pull each missing image from its Docker URI and save as a local SIF."""
    for name, uri, sif_path in uris_to_pull:
        print("Pulling '{}'\n  → {}".format(uri, sif_path))
        try:
            subprocess.check_call(['singularity', 'pull', '--force', sif_path, uri])
        except subprocess.CalledProcessError as e:
            err("Warning: failed to pull '{}': {}".format(uri, e))


# Entry-point called by the CLI

def cache(sub_args, repo_path):
    """Entry point for `viralrecon cache`: pull every image in config/containers.json."""
    sif_cache = sub_args.sif_cache
    images_config = os.path.join(repo_path, 'config', 'containers.json')

    prepare_cache(sif_cache)

    if getattr(sub_args, 'dry_run', False):
        missing = missing_images(images_config, sif_cache)
        if missing:
            print("Would pull {} image(s):".format(len(missing)))
            for name, uri, _ in missing:
                print("  {} : {}".format(name, uri))
        else:
            print("All images already cached.")
        return

    missing = missing_images(images_config, sif_cache)
    if not missing:
        print("All images already cached in: {}".format(sif_cache))
        return
    pull_images(sif_cache, missing)
    print("Done. {} image(s) cached.".format(len(missing)))
