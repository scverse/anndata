"""
Code for merging/ concatenating AnnData objects.
"""
from collections.abc import Mapping
from functools import singledispatch, reduce
from typing import Callable, Collection, TypeVar, Union

import numpy as np
from scipy import sparse

T = TypeVar("T")

###################
# Utilities
###################


class MissingVal:
    """Represents a missing value."""

    pass


def is_missing(v) -> bool:
    return v is MissingVal


def not_missing(v) -> bool:
    return v is not MissingVal


# Since it's difficult to check equality of sparse arrays
@singledispatch
def equal(a, b) -> bool:
    return np.array_equal(a, b)


@equal.register(sparse.spmatrix)
def equal_sparse(a, b) -> bool:
    if isinstance(b, sparse.spmatrix):
        return len((a != b).data) == 0
    else:
        return False


def union_keys(ds: Collection[Mapping]) -> set:
    return set().union(*(d.keys() for d in ds))


def intersect_keys(ds: Collection[Mapping]) -> set:
    return reduce(lambda x, y: x.intersection(y), map(set, (d.keys() for d in ds)))


###################
# Per element logic
###################


def unique_value(vals: Collection[T]) -> Union[T, MissingVal]:
    """
    Given a collection vals, returns the unique value (if one exists), otherwise
    returns MissingValue.
    """
    unique_val = vals[0]
    for v in vals[1:]:
        if not equal(v, unique_val):
            return MissingVal
    return unique_val


def first(vals: Collection[T]) -> Union[T, MissingVal]:
    """
    Given a collection of vals, return the first non-missing one.If they're all missing,
    return MissingVal.
    """
    for val in vals:
        if not_missing(val):
            return val
    return MissingVal


###################
# Merging
###################


def merge_nested(ds: Collection[Mapping], keys_join: Callable, value_join: Callable):
    out = {}
    for k in keys_join(ds):
        v = _merge_nested(ds, k, keys_join, value_join)
        if not_missing(v):
            out[k] = v
    return out


def _merge_nested(
    ds: Collection[Mapping], k, keys_join: Callable, value_join: Callable
):
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


def merge_unique(ds: Collection[Mapping]) -> Mapping:
    return merge_nested(ds, union_keys, unique_value)


def merge_same(ds: Collection[Mapping]) -> Mapping:
    return merge_nested(ds, intersect_keys, unique_value)


def merge_first(ds: Collection[Mapping]) -> Mapping:
    return merge_nested(ds, union_keys, first)


###################
# Interface
###################

# Leaving out for now, it's ugly in the rendered docs and would be adding a dependency.
# from typing_extensions import Literal
# UNS_STRATEGIES_TYPE = Literal[None, "same", "unique", "first"]
UNS_STRATEGIES = {
    None: lambda x: {},
    "same": merge_same,
    "unique": merge_unique,
    "first": merge_first,
}

# TODO: I should be making copies of all sub-elements
# TODO: Should I throw a warning about sparse arrays in uns?
def merge_uns(unss, strategy):
    return UNS_STRATEGIES[strategy](unss)
