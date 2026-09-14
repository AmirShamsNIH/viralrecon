#!/usr/bin/env python3
"""Resolve container image paths for a platform. containers.json names per-platform
{shared} and {ours} roots, so the image file names (the version pins) are kept once."""

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
    """The {shared} and {ours} directories for platform, with {repo_parent} expanded.
    An unknown platform falls back to the default so the error names a real path."""
    data = data if data is not None else load_container_config(repo_path)
    roots = data.get("roots", {})
    entry = roots.get(platform) or roots.get(DEFAULT_PLATFORM) or {}
    repo_parent = os.path.dirname(_repo_root(repo_path))
    return {
        key: os.path.normpath(val.format(repo_parent=repo_parent))
        for key, val in entry.items()
    }


def resolve_images(repo_path, platform=DEFAULT_PLATFORM, data=None):
    """tool name -> absolute .sif path for platform, without the "_" commentary keys."""
    data  = data if data is not None else load_container_config(repo_path)
    roots = image_roots(repo_path, platform, data)
    out   = {}
    for tool, path in data.get("images", {}).items():
        if tool.startswith("_"):
            continue
        out[tool] = os.path.normpath(path.format(**roots))
    return out
