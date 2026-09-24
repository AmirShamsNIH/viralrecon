# common.smk: shared helpers and wildcard constraints

from scripts.common import allocated

# Names must not contain path separators.
wildcard_constraints:
    sample = r"[^/]+",
    target = r"[^/]+",
    # Constrained so {vset} cannot absorb part of a target or file name.
    vset = r"raw|filtered",
