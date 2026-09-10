#!/usr/bin/env python3
"""
containers.py
─────────────
Resolve the container image map for the platform the pipeline is running on.

Why this is not just a dict of absolute paths
---------------------------------------------
Every image used to be written out in full in config/containers.json, which
worked while Biowulf was the only platform: the shared library lives at
/data/OpenOmics/SIFs and images we build ourselves live under
/data/RTB_GRS/references/singularity, and both are fixed institutional paths.

Neither path exists on another cluster, and on a cluster where the pipeline is
simply cloned into a group directory there is no institutional reference tree to
put images in at all. So containers.json now names two *roots* per platform and
gives each image as "{shared}/..." or "{ours}/...". Only the roots change
between platforms; the image file names, which are the version pins of record,
stay in one place and cannot drift apart per platform.

A root may contain "{repo_parent}", which expands to the directory holding the
clone. That keeps a platform's built images beside the checkout instead of in a
separate tree far away -- on a new cluster the whole deployment is then one
directory to find, move or hand over.
"""

import json
import os


DEFAULT_PLATFORM = "BIOWULF"


def _repo_root(repo_path):
    return os.path.abspath(repo_path)


def load_container_config(repo_path):
    """Read config/containers.json as-is."""
    cfg = os.path.join(_repo_root(repo_path), "config", "containers.json")
    with open(cfg) as fh:
        return json.load(fh)


def image_roots(repo_path, platform=DEFAULT_PLATFORM, data=None):
    """
    The {shared} and {ours} directories for `platform`, with {repo_parent}
    expanded. Falls back to the default platform's roots when a platform has
    no entry, so an unknown --platform fails on a missing image with a real
    path in the message rather than on a KeyError here.
    """
    data = data if data is not None else load_container_config(repo_path)
    roots = data.get("roots", {})
    entry = roots.get(platform) or roots.get(DEFAULT_PLATFORM) or {}
    repo_parent = os.path.dirname(_repo_root(repo_path))
    return {
        key: os.path.normpath(val.format(repo_parent=repo_parent))
        for key, val in entry.items()
    }


def resolve_images(repo_path, platform=DEFAULT_PLATFORM, data=None):
    """
    tool name -> absolute .sif path for `platform`.

    Keys beginning with "_" are commentary and are dropped, so callers get a
    map they can index by tool name without filtering.
    """
    data  = data if data is not None else load_container_config(repo_path)
    roots = image_roots(repo_path, platform, data)
    out   = {}
    for tool, path in data.get("images", {}).items():
        if tool.startswith("_"):
            continue
        out[tool] = os.path.normpath(path.format(**roots))
    return out
