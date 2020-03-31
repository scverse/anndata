"""
Code for merging/ concatenating AnnData objects.
"""
from collections.abc import Mapping
from functools import singledispatch, reduce
from typing_extensions import Literal

import numpy as np
from scipy import sparse

#################
# Utilities
#################


class MissingVal:
    """Represents a missing value."""

    pass


def is_missing(v):
    return v is MissingVal


def not_missing(v):
    return v is not MissingVal


@singledispatch
def equal(a, b):
    return np.array_equal(a, b)


@equal.register(sparse.spmatrix)
def equal_sparse(a, b):
    if isinstance(b, sparse.spmatrix):
        return len((a != b).data) == 0
    else:
        return False


def union_keys(ds):
    return set().union(*(d.keys() for d in ds))


def intersect_keys(ds):
    return reduce(lambda x, y: x.intersection(y), map(set, (d.keys() for d in ds)))


###################
# Per element logic
###################


def unique_value(vals):
    """Given a collection vals, returns the unique value (if one exists), otherwise returns MissingValue."""
    unique_val = vals[0]
    for v in vals[1:]:
        if not equal(v, unique_val):
            return MissingVal
    return unique_val


def first(vals):
    for val in vals:
        if not_missing(val):
            return val
    return MissingVal


###################
# Merging
###################

def merge_nested(ds, keys_join, value_join):
    out = {}
    for k in keys_join(ds):
        v = _merge_nested(ds, k, keys_join, value_join)
        if not_missing(v):
            out[k] = v
    return out


def _merge_nested(ds, k, keys_join, value_join):
    vals = [d[k] for d in ds if k in d]
    if len(vals) == 0:
        return MissingVal
    elif all(isinstance(v, Mapping) for v in vals):
        new_map = merge_nested(vals, keys_join, value_join)
        if len(new_map) == 0:
            return MissingVal
        else:
            return new_map
    else:
        return value_join(vals)


def merge_unique(ds):
    return merge_nested(ds, union_keys, unique_value)


def merge_common(ds):
    return merge_nested(ds, intersect_keys, unique_value)


def merge_first(ds):
    return merge_nested(ds, union_keys, first)

######################
# Interface
######################


UNS_STRATEGIES_TYPE = Literal[None, "common", "unique", "first"]
UNS_STRATEGIES = {
    None: lambda x: {},
    "common": merge_common,
    "unique": merge_unique,
    "first": merge_first,
}

# TODO: I should be making copies of all sub-elements
# TODO: Should I throw a warning about sparse arrays in uns?
def merge_uns(unss, strategy: UNS_STRATEGIES_TYPE):
    return UNS_STRATEGIES[strategy](unss)
