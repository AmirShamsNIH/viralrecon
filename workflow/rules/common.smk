# common.smk: shared helpers and wildcard constraints

from scripts.common import allocated

# Wildcard constraints: sample names must not contain path separators.
wildcard_constraints:
    sample = r"[^/]+",
    target = r"[^/]+",
    # the variant sets carried in parallel by make_variants_matrix; constrained
    # so {vset} cannot absorb part of a target or file name
    vset   = r"raw|filtered",
