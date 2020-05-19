"""
Code for merging/ concatenating AnnData objects.
"""
from collections import OrderedDict
from collections.abc import Mapping, MutableSet
from copy import deepcopy
from functools import partial, reduce, singledispatch
from itertools import repeat
from operator import and_, or_, sub
from typing import Callable, Collection, Iterable, Tuple, TypeVar, Union
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


#####################
# Concatenation
#####################


class Reindexer(object):
    """
    Indexing to be applied to axis of 2d array orthogonal to the axis being concatenated.

    Attrs
    -----
    old_idx
        Original index
    new_idx
        Target index
    old_pos
        Indices of original index which will be kept
    new_pos
        Indices of new index which data from old_pos will be placed in.
        Together with `old_pos` this forms a mapping.
    """

    def __init__(self, old_idx, new_idx):
        self.old_idx = old_idx
        self.new_idx = new_idx

        new_pos = new_idx.get_indexer(old_idx)
        old_pos = np.arange(len(new_pos))

        mask = new_pos != -1

        self.new_pos = new_pos[mask]
        self.old_pos = old_pos[mask]

    def __call__(self, el, *, axis=1, fill_value=None):
        return self.apply(el, axis=axis, fill_value=fill_value)

    def apply(self, el, *, axis, fill_value=None):
        if isinstance(el, pd.DataFrame):
            return self._apply_to_df(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, sparse.spmatrix):
            return self._apply_to_sparse(el, axis=axis, fill_value=fill_value)
        else:
            return self._apply_to_array(el, axis=axis, fill_value=fill_value)

    def _apply_to_df(self, el, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = np.NaN
        return el.reindex(self.new_idx, axis=axis, fill_value=fill_value)

    def _apply_to_array(self, el, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = default_fill_value([el])
        if 0 in el.shape:
            return np.broadcast_to(fill_value, (el.shape[0], len(self.new_idx)))

        indexer = self.old_idx.get_indexer(self.new_idx)

        # Indexes real fast, and does outer indexing
        return pd.api.extensions.take(
            el, indexer, axis=axis, allow_fill=True, fill_value=fill_value,
        )

    def _apply_to_sparse(self, el, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = default_fill_value([el])
        if fill_value != 0:
            to_fill = self.new_idx.get_indexer(self.new_idx.difference(self.old_idx))
        else:
            to_fill = np.array([])

        # Fixing outer indexing for missing values
        if el.shape[1] == 0:
            to_fill = range(0, len(self.new_idx))
            el = sparse.coo_matrix((el.shape[0], len(self.old_pos)))

        fill_idxer = None

        if len(to_fill) > 0:
            idxmtx_dtype = np.promote_types(el.dtype, np.array(fill_value).dtype)
        else:
            idxmtx_dtype = bool

        if axis == 1:
            idxmtx = sparse.coo_matrix(
                (np.ones(len(self.new_pos), dtype=bool), (self.old_pos, self.new_pos)),
                shape=(len(self.old_idx), len(self.new_idx)),
                dtype=idxmtx_dtype,
            )
            out = el @ idxmtx

            if len(to_fill) > 0:
                out = out.tocsc()
                fill_idxer = (slice(None), to_fill)
        elif axis == 0:
            idxmtx = sparse.coo_matrix(
                (np.ones(len(self.new_pos), dtype=bool), (self.new_pos, self.old_pos)),
                shape=(len(self.new_idx), len(self.old_idx)),
                dtype=idxmtx_dtype,
            )
            out = idxmtx @ el

            if len(to_fill) > 0:
                out = out.tocsr()
                fill_idxer = (to_fill, slice(None))

        if fill_idxer is not None:
            out[fill_idxer] = fill_value

        return out


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


def gen_reindexer(new_var: pd.Index, cur_var: pd.Index):
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
    """
    return Reindexer(cur_var, new_var)


def concat_arrays(arrays, reindexers, axis=0, index=None, fill_value=None):
    arrays = list(arrays)
    if fill_value is None:
        fill_value = default_fill_value(arrays)

    if any(isinstance(a, pd.DataFrame) for a in arrays):
        if not all(isinstance(a, pd.DataFrame) for a in arrays):
            raise NotImplementedError(
                "Cannot concatenate a dataframe with other array types."
            )
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat(
            [f(x) for f, x in zip(reindexers, arrays)], ignore_index=True, axis=axis
        )
        df.index = index
        return df
    elif any(isinstance(a, sparse.spmatrix) for a in arrays):
        sparse_stack = (sparse.vstack, sparse.hstack)[axis]
        return sparse_stack(
            [
                f(as_sparse(a), axis=1 - axis, fill_value=fill_value)
                for f, a in zip(reindexers, arrays)
            ],
            format="csr",
        )
    else:
        return np.concatenate(
            [
                f(x, fill_value=fill_value, axis=1 - axis)
                for f, x in zip(reindexers, arrays)
            ],
            axis=axis,
        )


def inner_concat_aligned_mapping(mappings, reindexers=None, index=None, axis=0):
    result = {}

    for k in intersect_keys(mappings):
        els = [m[k] for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_inner_reindexers(els, new_index=index, axis=axis)
        else:
            cur_reindexers = reindexers

        result[k] = concat_arrays(els, cur_reindexers, index=index, axis=axis)
    return result


def gen_inner_reindexers(els, new_index, axis=0):
    alt_axis = 1 - axis
    if axis == 0:
        df_indices = lambda x: x.columns
    elif axis == 1:
        df_indices = lambda x: x.indices

    if all(isinstance(el, pd.DataFrame) for el in els if not_missing(el)):
        common_ind = reduce(
            lambda x, y: x.intersection(y), (df_indices(el) for el in els)
        )
        reindexers = [Reindexer(df_indices(el), common_ind) for el in els]
    else:
        min_ind = min(el.shape[alt_axis] for el in els)
        reindexers = [
            gen_reindexer(pd.RangeIndex(min_ind), pd.RangeIndex(el.shape[alt_axis]))
            for el in els
        ]
    return reindexers


def gen_outer_reindexers(els, shapes, new_index: pd.Index, *, axis=0):
    if all(isinstance(el, pd.DataFrame) for el in els if not_missing(el)):
        reindexers = [
            lambda x: x if not_missing(el) else pd.DataFrame(index=range(shape))
            for el, shape in zip(els, shapes)
        ]
    else:
        # if fill_value is None:
        # fill_value = default_fill_value(els)

        max_col = max(el.shape[1] for el in els if not_missing(el))
        orig_cols = [el.shape[1] if not_missing(el) else 0 for el in els]
        reindexers = [
            gen_reindexer(pd.RangeIndex(max_col), pd.RangeIndex(n),) for n in orig_cols
        ]
    return reindexers


def outer_concat_aligned_mapping(
    mappings, reindexers=None, index=None, fill_value=None, axis=0
):
    result = {}
    ns = [m.parent.shape[axis] for m in mappings]

    for k in union_keys(mappings):
        els = [m.get(k, MissingVal) for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_outer_reindexers(els, ns, new_index=index, axis=axis,)
        else:
            cur_reindexers = reindexers

        result[k] = concat_arrays(
            [
                el if not_missing(el) else np.zeros((n, 0), dtype=bool)
                for el, n in zip(els, ns)
            ],
            cur_reindexers,
            axis=axis,
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


def _resolve_dim(*, dim: str = None, axis: int = None) -> Tuple[int, str]:
    _dims = ("obs", "var")
    if (dim is None and axis is None) or (dim is not None and axis is not None):
        raise ValueError(
            f"Must pass exactly one of `dim` or `axis`. Got: dim={dim}, axis={axis}."
        )
    elif dim is not None and dim not in _dims:
        raise ValueError(f"`dim` must be one of ('obs', 'var'), was {dim}")
    elif axis is not None and axis not in (0, 1):
        raise ValueError(f"`axis` must be either 0 or 1, was {axis}")
    if dim is not None:
        return _dims.index(dim), dim
    else:
        return axis, _dims[axis]


def dim_indices(adata, *, axis=None, dim=None):
    _, dim = _resolve_dim(axis=axis, dim=dim)
    return getattr(adata, f"{dim}_names")


def dim_size(adata, *, axis=None, dim=None):
    ax, _ = _resolve_dim(axis, dim)
    return adata.shape[ax]


def concat(
    adatas,
    *,
    axis=0,
    join="inner",
    batch_key="batch",
    batch_categories=None,
    uns_merge=None,
    merge=None,
    index_unique="-",
    fill_value=None,
):
    """Re-implementation of `AnnData.concatenate`, but better."""
    # Argument normalization
    adatas = list(adatas)
    if batch_categories is None:
        batch_categories = np.arange(len(adatas)).astype(str)
    if axis == 0:
        dim = "obs"
    elif axis == 1:
        dim = "var"

    alt_axis, alt_dim = _resolve_dim(axis=1 - axis)

    # Batch column
    batch = pd.Categorical.from_codes(
        np.repeat(np.arange(len(adatas)), [a.shape[axis] for a in adatas]),
        categories=batch_categories,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(dim_indices(a, axis=axis)) for a in adatas], ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(batch.map(str), sep=index_unique)
    concat_indices = pd.Index(concat_indices)

    alt_indices = resolve_index(
        [dim_indices(a, axis=1 - axis) for a in adatas], join=join
    )
    reindexers = [
        gen_reindexer(alt_indices, dim_indices(a, axis=1 - axis)) for a in adatas
    ]

    # Annotation for concatenation axis
    concat_annot = pd.concat([getattr(a, dim) for a in adatas], ignore_index=True)
    concat_annot.index = concat_indices
    concat_annot[batch_key] = batch

    # Annotation for other axis
    alt_annot = merge_dataframes(
        [getattr(a, alt_dim) for a in adatas],
        alt_indices,
        # TODO: Allow use of other strategies, like passing merge_unique here.
        # Current behaviour is mostly for backwards compat. It's like make_names_unique, but
        # unfortunately the behaviour is different.
        partial(merge_outer, batch_keys=batch_categories, merge=merge_same),
    )

    X = concat_arrays(
        [a.X for a in adatas], reindexers, axis=axis, fill_value=fill_value
    )

    if join == "inner":
        layers = inner_concat_aligned_mapping(
            [a.layers for a in adatas], axis=axis, reindexers=reindexers
        )
        concat_mapping = inner_concat_aligned_mapping(
            [getattr(a, f"{dim}m") for a in adatas], index=concat_indices
        )
        # obsm = inner_concat_aligned_mapping([a.obsm for a in adatas], index=obs_names)
    elif join == "outer":
        layers = outer_concat_aligned_mapping(
            [a.layers for a in adatas], reindexers, axis=axis, fill_value=fill_value
        )
        concat_mapping = outer_concat_aligned_mapping(
            [getattr(a, f"{dim}m") for a in adatas],
            index=concat_indices,
            fill_value=fill_value,
        )

    # Need to do reindexing along other axis, make reindexer work for dataframes
    alt_mapping = merge_uns(
        [
            {k: r(v, axis=0) for k, v in getattr(a, f"{alt_dim}m").items()}
            for r, a in zip(reindexers, adatas)
        ],
        strategy=merge,
    )
    # varm = merge_uns([a.varm for a in adatas], strategy=merge)
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
            axis=axis,
        )
    elif any(has_raw):
        warn(
            "Only some AnnData objects have `.raw` attribute, "
            "not concatenating `.raw` attributes.",
            UserWarning,
        )
    return AnnData(
        **{
            "X": X,
            "layers": layers,
            dim: concat_annot,
            alt_dim: alt_annot,
            f"{dim}m": concat_mapping,
            f"{alt_dim}m": alt_mapping,
            "uns": uns,
            "raw": raw,
        }
    )
