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
import warnings
from warnings import warn

import numpy as np
import pandas as pd
from scipy import sparse

from .anndata import AnnData

T = TypeVar("T")

###################
# Utilities
###################


# Pretty much just for maintaining order of keys
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
        return reduce(or_, vals, self)

    def discard(self, val):
        if val in self:
            del self.dict[val]

    def difference(self, *vals):
        return reduce(sub, vals, self)


def union_keys(ds: Collection) -> OrderedSet:
    return reduce(or_, ds, OrderedSet())


def intersect_keys(ds: Collection) -> OrderedSet:
    return reduce(and_, map(OrderedSet, ds))


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


def as_sparse(x):
    if not isinstance(x, sparse.spmatrix):
        return sparse.csr_matrix(x)
    else:
        return x


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


def resolve_index(inds: Iterable[pd.Index], join):
    if join == "inner":
        return reduce(lambda x, y: x.intersection(y), inds)
    elif join == "outer":
        return reduce(lambda x, y: x.union(y), inds)
    else:
        raise ValueError()


def default_fill_value(els):
    """Given some arrays, returns what the default fill value should be.

    This is largely due to backwards compat, and might not be the ideal solution.
    """
    if any(isinstance(el, sparse.spmatrix) for el in els):
        return 0
    else:
        return np.nan


def gen_reindexer(new_var: pd.Index, cur_var: pd.Index, *, fill_value=0):
    """
    Given a new set of var_names, and a current set, generates a function which will reindex
    a matrix to be aligned with the new set.

    Usage
    -----

    >>> a = AnnData(sparse.eye(3), var=pd.DataFrame(index=list("abc")))
    >>> b = AnnData(sparse.eye(2), var=pd.DataFrame(index=list("ba")))
    >>> reindexer = gen_reindexer(a.var_names, b.var_names)
    >>> sparse.vstack([a.X, reindexer(b.X)]).toarray()
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.],
           [0., 1., 0.],
           [1., 0., 0.]], dtype=float32)
    >>> reindexer_nan = gen_reindexer(a.var_names, b.var_names, fill_value=np.nan)
    >>> sparse.vstack([a.X, reindexer_nan(b.X)]).toarray()
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.],
           [ 0.,  1., nan],
           [ 1.,  0., nan]], dtype=float32)
    """
    new_size = len(new_var)
    old_size = len(cur_var)
    new_pts = new_var.get_indexer(cur_var)
    cur_pts = np.arange(len(new_pts))

    mask = new_pts != -1

    new_pts = new_pts[mask]
    cur_pts = cur_pts[mask]

    def reindexer(X, fill_value=fill_value):
        if not np.can_cast(fill_value, X.dtype):
            out_dtype = np.promote_types(np.array(fill_value).dtype, X.dtype)
        else:
            out_dtype = X.dtype

        # Case where layer was empty
        # TODO: This could be done more elegantly
        if X.shape[1] != 0:
            idxmtx = sparse.coo_matrix(
                (np.ones(len(new_pts), dtype=int), (cur_pts, new_pts)),
                shape=(old_size, new_size),
                dtype=out_dtype,
            )
        else:
            idxmtx = sparse.coo_matrix((0, new_size))
        out = X @ idxmtx

        if fill_value != 0:
            if X.shape[1] != 0:
                to_fill = new_var.get_indexer(new_var.difference(cur_var))
            else:
                to_fill = range(0, len(new_var))
            if len(to_fill) > 0:
                # More efficient to set columns on csc
                if sparse.issparse(out):
                    out = sparse.csc_matrix(out)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
                    out[:, to_fill] = fill_value

        return out

    return reindexer


def concat_arrays(arrays, reindexers, index=None, fill_value=None):
    arrays = list(arrays)
    if fill_value is None:
        fill_value = default_fill_value(arrays)

    if any(isinstance(a, pd.DataFrame) for a in arrays):
        if not all(isinstance(a, pd.DataFrame) for a in arrays):
            raise NotImplementedError(
                "Cannot concatenate a dataframe with other array types."
            )
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat([f(x) for f, x in zip(reindexers, arrays)], ignore_index=True)
        df.index = index
        return df
    elif any(isinstance(a, sparse.spmatrix) for a in arrays):
        return sparse.vstack(
            [
                f(as_sparse(a), fill_value=fill_value)
                for f, a in zip(reindexers, arrays)
            ],
            format="csr",
        )
    else:
        return np.vstack(
            [f(x, fill_value=fill_value) for f, x in zip(reindexers, arrays)]
        )


def inner_concat_aligned_mapping(mappings, reindexers=None, index=None):
    result = {}

    for k in intersect_keys(mappings):
        els = [m[k] for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_inner_reindexers(els, new_index=index)
        else:
            cur_reindexers = reindexers

        result[k] = concat_arrays(els, cur_reindexers, index=index)
    return result


def gen_inner_reindexers(els, new_index):
    if all(isinstance(el, pd.DataFrame) for el in els if not_missing(el)):
        common_cols = reduce(lambda x, y: x.intersection(y), (el.columns for el in els))
        reindexers = repeat(lambda x: x[common_cols])
    else:
        min_col = min(el.shape[1] for el in els)
        reindexers = [
            gen_reindexer(pd.RangeIndex(min_col), pd.RangeIndex(el.shape[1]))
            for el in els
        ]
    return reindexers


def gen_outer_reindexers(els, shapes, new_index: pd.Index, fill_value=None):
    if all(isinstance(el, pd.DataFrame) for el in els if not_missing(el)):
        reindexers = [
            lambda x: x if not_missing(el) else pd.DataFrame(index=range(shape))
            for el, shape in zip(els, shapes)
        ]
    else:
        if fill_value is None:
            fill_value = default_fill_value(els)

        max_col = max(el.shape[1] for el in els if not_missing(el))
        orig_cols = [el.shape[1] if not_missing(el) else 0 for el in els]
        reindexers = [
            gen_reindexer(
                pd.RangeIndex(max_col), pd.RangeIndex(n), fill_value=fill_value,
            )
            for n in orig_cols
        ]
    return reindexers


def outer_concat_aligned_mapping(
    mappings, reindexers=None, index=None, fill_value=None
):
    result = {}
    ns = [m.parent.n_obs for m in mappings]

    for k in union_keys(mappings):
        els = [m.get(k, MissingVal) for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_outer_reindexers(
                els, ns, new_index=index, fill_value=fill_value
            )
        else:
            cur_reindexers = reindexers

        result[k] = concat_arrays(
            [
                el if not_missing(el) else np.zeros((n, 0), dtype=bool)
                for el, n in zip(els, ns)
            ],
            cur_reindexers,
            index=index,
            fill_value=fill_value,
        )
    return result


def merge_dataframes(dfs, new_index, merge_strategy=merge_unique):
    dfs = [df.reindex(index=new_index) for df in dfs]
    # New dataframe with all shared data
    new_df = pd.DataFrame(merge_strategy(dfs), index=new_index)
    return new_df


def merge_outer(mappings, batch_keys, *, join_index="-", merge=merge_unique):
    """
    Combine elements of two mappings, such that non-overlapping entries are added with their batch-key appended.

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
    fill_value=None,
):
    """Re-implementation of `AnnData.concatenate`, but better."""
    # Argument normalization
    adatas = list(adatas)
    if batch_categories is None:
        batch_categories = np.arange(len(adatas)).astype(str)

    # Batch column
    batch = pd.Categorical.from_codes(
        np.repeat(np.arange(len(adatas)), [a.n_obs for a in adatas]),
        categories=batch_categories,
    )

    # Combining indexes
    obs_names = pd.concat([pd.Series(a.obs_names) for a in adatas], ignore_index=True)
    if index_unique is not None:
        obs_names = obs_names.str.cat(batch.map(str), sep=index_unique)
    obs_names = pd.Index(obs_names)

    var_names = resolve_index([a.var_names for a in adatas], join=join)
    reindexers = [
        gen_reindexer(var_names, a.var_names, fill_value=fill_value) for a in adatas
    ]

    # Obs
    obs = pd.concat([a.obs for a in adatas], ignore_index=True)
    obs.index = obs_names
    obs[batch_key] = batch

    # Var
    var = merge_dataframes(
        [a.var for a in adatas],
        var_names,
        # TODO: Allow use of other strategies, like passing merge_unique here.
        # Current behaviour is mostly for backwards compat. It's like make_names_unique, but
        # unfortunately the behaviour is different.
        partial(merge_outer, batch_keys=batch_categories, merge=merge_same),
    )

    X = concat_arrays([a.X for a in adatas], reindexers, fill_value=fill_value)

    if join == "inner":
        layers = inner_concat_aligned_mapping([a.layers for a in adatas], reindexers)
        obsm = inner_concat_aligned_mapping([a.obsm for a in adatas], index=obs_names)
    elif join == "outer":
        layers = outer_concat_aligned_mapping(
            [a.layers for a in adatas], reindexers, fill_value=fill_value
        )
        obsm = outer_concat_aligned_mapping(
            [a.obsm for a in adatas], index=obs_names, fill_value=fill_value
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
            fill_value=fill_value,
        )
    elif any(has_raw):
        warn(
            "Only some AnnData objects have `.raw` attribute, "
            "not concatenating `.raw` attributes.",
            UserWarning,
        )

    return AnnData(X=X, layers=layers, obs=obs, var=var, obsm=obsm, uns=uns, raw=raw)
