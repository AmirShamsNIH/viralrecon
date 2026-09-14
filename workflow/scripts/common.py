#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Helper functions shared across the viralrecon workflow."""


def allocated(resource, rule, lookup, default='__default__'):
    """Return a resource value ('threads', 'mem', 'time', ...) for rule from the
    parsed cluster.json, falling back to the __default__ entry."""
    try:
        return lookup[rule][resource]
    except KeyError:
        return lookup[default][resource]


def provided(sample_list, condition):
    """Return sample_list when condition is True, otherwise an empty list."""
    return sample_list if condition else []


def ignore(sample_list, condition):
    """Return empty list when condition is True (inverse of provided())."""
    return [] if condition else sample_list


def references(config, ref_keys):
    """True only if every key in ref_keys is present and non-empty in
    config['references']."""
    for key in ref_keys:
        val = config.get('references', {}).get(key, '')
        if not val:
            return False
    return True


def str_bool(s):
    """Cast a string to bool: true/1/y/yes or false/0/n/no/empty. Raises TypeError
    for anything else."""
    v = str(s).lower()
    if v in ('true', '1', 'y', 'yes'):
        return True
    if v in ('false', '0', 'n', 'no', ''):
        return False
    raise TypeError("Cannot convert '{}' to bool.".format(s))
