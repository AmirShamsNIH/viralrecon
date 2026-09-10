#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Common helper functions shared across the entire viralrecon workflow.
Mirrors the OpenOmics/baseline workflow/scripts/common.py pattern.
"""


def allocated(resource, rule, lookup, default='__default__'):
    """Return a resource value for a rule from cluster.json.
    Falls back to __default__ if the rule has no explicit entry.

    @param resource <str>: 'threads', 'mem', 'time', 'partition', 'gres', etc.
    @param rule     <str>: Snakemake rule name (must match cluster.json key)
    @param lookup  <dict>: parsed cluster.json dict
    @param default  <str>: fallback key [default: '__default__']
    @return <str>
    """
    try:
        return lookup[rule][resource]
    except KeyError:
        return lookup[default][resource]


def provided(sample_list, condition):
    """Return sample_list when condition is True, empty list otherwise.
    Used in rule all to conditionally include optional output targets.
    """
    return sample_list if condition else []


def ignore(sample_list, condition):
    """Return empty list when condition is True (inverse of provided())."""
    return [] if condition else sample_list


def references(config, ref_keys):
    """Return True only if every key in ref_keys exists and is non-empty in
    config['references']. Used to guard rules that need optional references.
    """
    for key in ref_keys:
        val = config.get('references', {}).get(key, '')
        if not val:
            return False
    return True


def str_bool(s):
    """Safely cast a string to bool.
    Accepts: 'true'/'1'/'y'/'yes' → True; 'false'/'0'/'n'/'no'/'' → False.
    Raises TypeError for anything else.
    """
    v = str(s).lower()
    if v in ('true', '1', 'y', 'yes'):
        return True
    if v in ('false', '0', 'n', 'no', ''):
        return False
    raise TypeError("Cannot convert '{}' to bool.".format(s))
