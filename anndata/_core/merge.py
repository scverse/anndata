"""
Code for merging/ concatenating AnnData objects.
"""
from collections import OrderedDict
from collections.abc import Mapping, MutableSet
from copy import deepcopy
from functools import partial, reduce, singledispatch
from itertools import repeat
from operator import and_, or_, sub
from typing import Callable, Collection, Iterable, TypeVar, Union
from warnings import warn

import numpy as np
import pandas as pd
from scipy import sparse

from .anndata import AnnData

T = TypeVar("T")

###################
# Utilities
###################


class OrderedSet(MutableSet):
    def __init__(self, vals=()):
        self.dict = OrderedDict(zip(vals, repeat(None)))

    def __contains__(self, val):
        return val in self.dict

    def __iter__(self):
        return iter(self.dict)

    def __len__(self):
        return len(self.dict)

    def __repr__(self):
        return "OrderedSet: {" + ", ".join(map(str, self)) + "}"

    def copy(self):
        return OrderedSet(self.dict.copy())

    def add(self, val):
        self.dict[val] = None

    def union(self, *vals):
        return reduce(or_, (self, *vals))

    def discard(self, val):
        if val in self:
            del self.dict[val]

    def difference(self, *vals):
        return reduce(sub, (self, *vals))


class MissingVal:
    """Represents a missing value."""


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
    # It's a weird api, don't blame me
    if isinstance(b, sparse.spmatrix):
        comp = a != b
        if isinstance(comp, bool):
            return not comp
        else:
            return len((a != b).data) == 0
    else:
        return False


def union_keys(ds: Collection[Mapping]) -> set:
    return reduce(or_, (d.keys() for d in ds), OrderedSet())


def intersect_keys(ds: Collection[Mapping]) -> set:
    return reduce(and_, map(OrderedSet, (d.keys() for d in ds)))


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


def only(vals: Collection[T]) -> Union[T, MissingVal]:
    """Return the only value in the collection, otherwise MissingVal."""
    if len(vals) == 1:
        return vals[0]
    else:
        return MissingVal


###################
# Merging
###################


def merge_nested(ds: Collection[Mapping], keys_join: Callable, value_join: Callable):
    out = {}
    for k in keys_join(ds):
        v = _merge_nested(ds, k, keys_join, value_join)
        if not_missing(v):
            out[k] = deepcopy(v)
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


def merge_only(ds: Collection[Mapping]) -> Mapping:
    return merge_nested(ds, union_keys, only)


###################
# Interface
###################

# Leaving out for now, it's ugly in the rendered docs and would be adding a dependency.
# from typing_extensions import Literal
# UNS_STRATEGIES_TYPE = Literal[None, "same", "unique", "first", "only"]
UNS_STRATEGIES = {
    None: lambda x: {},
    "same": merge_same,
    "unique": merge_unique,
    "first": merge_first,
    "only": merge_only,
}

# TODO: I should be making copies of all sub-elements
# TODO: Should I throw a warning about sparse arrays in uns?
def merge_uns(unss, strategy):
    return UNS_STRATEGIES[strategy](unss)


def as_sparse(x):
    if not isinstance(x, sparse.spmatrix):
        return sparse.csr_matrix(x)
    else:
        return x


def resolve_index(inds: Iterable[pd.Index], join):
    if join == "inner":
        return reduce(lambda x, y: x.intersection(y), inds)
    elif join == "outer":
        return reduce(lambda x, y: x.union(y), inds)
    else:
        raise ValueError()


def gen_reindexer(new_var: pd.Index, cur_var: pd.Index):
    """
    Given a new set of var_names, and a current set, generates a function which will reindex
    a matrix to be aligned with the new set.

    Usage
    -----

    >>> reindexer = gen_reindexer(a.var_names, b.var_names)
    >>> sparse.vstack([a.X, reindexer(b.X)])
    """
    new_size = len(new_var)
    old_size = len(cur_var)
    new_pts = new_var.get_indexer(cur_var)
    cur_pts = np.arange(len(new_pts))

    mask = new_pts != -1

    new_pts = new_pts[mask]
    cur_pts = cur_pts[mask]

    def reindexer(X):
        idxmtx = sparse.coo_matrix(
            (np.ones(len(new_pts), dtype=int), (cur_pts, new_pts)),
            shape=(old_size, new_size),
        )
        return X @ idxmtx

    return reindexer


def concat_arrays(arrays, reindexers, index=None):
    arrays = list(arrays)
    if any(isinstance(a, pd.DataFrame) for a in arrays):
        if not all(isinstance(a, pd.DataFrame) for a in arrays):
            return MissingVal
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat([f(x) for f, x in zip(reindexers, arrays)], ignore_index=True)
        df.index = index
        return df
    elif any(isinstance(a, sparse.spmatrix) for a in arrays):
        return sparse.vstack([f(as_sparse(a)) for f, a in zip(reindexers, arrays)])
    else:
        return np.vstack([f(x) for f, x in zip(reindexers, arrays)])


def concat_aligned_mapping(mappings, reindexers, index=None):
    return {
        k: concat_arrays([m[k] for m in mappings], reindexers, index=index)
        for k in intersect_keys(mappings)
    }


def merge_dataframes(dfs, new_index, merge_strategy=merge_unique):
    dfs = [df.reindex(index=new_index) for df in dfs]
    # New dataframe with all shared data
    new_df = pd.DataFrame(merge_strategy(dfs), index=new_index)
    return new_df


def merge_outer(mappings, batch_keys, *, join_index="-", merge=merge_unique):
    """
    Concate elements of two mappings, such that non-overlapping entries are added with their batch-key appended.

    Note: this currently does NOT work for nested mappings. Additionally, values are not promised to be unique, and may be overwritten.
    """
    all_keys = union_keys(mappings)
    out = merge(mappings)
    for key in all_keys.difference(out.keys()):
        for b, m in zip(batch_keys, mappings):
            val = m.get(key, None)
            if val is not None:
                out[f"{key}{join_index}{b}"] = val
    return out


def concat(
    adatas,
    *,
    join="inner",
    batch_key="batch",
    batch_categories=None,
    uns_merge=None,
    index_unique="-",
):
    """Re-implementation of `AnnData.concatenate`, but better."""
    # Argument normalization
    adatas = list(adatas)
    if batch_categories is None:
        batch_categories = np.arange(len(adatas)).astype(str)

    # Combining indexes
    obs_names = pd.Index(
        np.concatenate(
            [
                pd.Series(a.obs_names) + f"{index_unique}{batch}"
                for batch, a in zip(batch_categories, adatas)
            ]
        )
    )
    var_names = resolve_index([a.var_names for a in adatas], join=join)
    reindexers = [gen_reindexer(var_names, a.var_names) for a in adatas]

    # Obs
    batch = pd.Series(
        np.repeat(np.arange(len(adatas)), [a.n_obs for a in adatas]), dtype="category"
    ).map(dict(zip(np.arange(len(adatas)), batch_categories)))
    obs = pd.concat([a.obs for a in adatas], ignore_index=True)
    obs.index = obs_names
    obs[batch_key] = batch.values

    # Var
    var = merge_dataframes(
        [a.var for a in adatas],
        var_names,
        # TODO: Allow use of other strategies, like passing merge_unique here.
        # Current behaviour is mostly for backwards compat
        partial(merge_outer, batch_keys=batch_categories, merge=merge_same),
    )
    # The current behaviour for making column names unique is like make_index_unique, but unfortunatley does it's own thing.

    # Everything else
    X = concat_arrays([a.X for a in adatas], reindexers)
    layers = concat_aligned_mapping([a.layers for a in adatas], reindexers)
    # TODO: allow computing how to merge for each array seperately
    obsm = concat_aligned_mapping(
        [a.obsm for a in adatas], repeat(lambda x: x), index=obs_names
    )
    uns = merge_uns([a.uns for a in adatas], strategy=uns_merge)

    raw = None
    has_raw = [a.raw is not None for a in adatas]
    if all(has_raw):
        raw = concat(
            [
                AnnData(
                    X=a.raw.X,
                    obs=pd.DataFrame(index=a.obs_names),
                    var=a.raw.var,
                    varm=a.raw.varm,
                )
                for a in adatas
            ],
            join=join,
            batch_key=batch_key,
            batch_categories=batch_categories,
            index_unique=index_unique,
        )
    elif any(has_raw):
        warn(
            "Only some adata objects have `.raw` attribute, "
            "not concatenating `.raw` attributes.",
            UserWarning,
        )

    return AnnData(X=X, layers=layers, obsm=obsm, obs=obs, var=var, uns=uns, raw=raw)
