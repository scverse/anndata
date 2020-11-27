"""
Code for merging/ concatenating AnnData objects.
"""
from collections import OrderedDict
from collections.abc import Mapping, MutableSet
from functools import reduce, singledispatch
from itertools import repeat
from operator import and_, or_, sub
from typing import Any, Callable, Collection, Iterable, Optional, Tuple, TypeVar, Union
import typing
from warnings import warn

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.base import spmatrix

from .anndata import AnnData
from ..compat import Literal
from ..utils import asarray

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

    def union(self, *vals) -> "OrderedSet":
        return reduce(or_, vals, self)

    def discard(self, val):
        if val in self:
            del self.dict[val]

    def difference(self, *vals) -> "OrderedSet":
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


# We need to be able to check for equality of arrays to know which are the same.
# Unfortunatley equality of arrays is poorly defined.
# * `np.array_equal` does not work for sparse arrays
# * `np.array_equal(..., equal_nan=True)` does not work for null values at the moment
#   (see https://github.com/numpy/numpy/issues/16377)
# So we have to define it ourselves with these two issues in mind.
# TODO: Hopefully this will stop being an issue in the future and this code can be removed.
@singledispatch
def equal(a, b) -> bool:
    return np.array_equal(a, asarray(b))


@equal.register(pd.DataFrame)
def equal_dataframe(a, b) -> bool:
    return a.equals(b)


@equal.register(np.ndarray)
def equal_array(a, b) -> bool:
    return equal(pd.DataFrame(a), pd.DataFrame(asarray(b)))


@equal.register(sparse.spmatrix)
def equal_sparse(a, b) -> bool:
    # It's a weird api, don't blame me
    if isinstance(b, sparse.spmatrix):
        comp = a != b
        if isinstance(comp, bool):
            return not comp
        # fmt: off
        return (
            (len(comp.data) == 0)
            or (
                np.isnan(a[comp]).all()
                and np.isnan(b[comp]).all()
            )
        )
        # fmt: on
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


def merge_only(ds: Collection[Mapping]) -> Mapping:
    return merge_nested(ds, union_keys, only)


###################
# Interface
###################

# Leaving out for now, it's ugly in the rendered docs and would be adding a dependency.
# from typing_extensions import Literal
# UNS_STRATEGIES_TYPE = Literal[None, "same", "unique", "first", "only"]
MERGE_STRATEGIES = {
    None: lambda x: {},
    "same": merge_same,
    "unique": merge_unique,
    "first": merge_first,
    "only": merge_only,
}

StrategiesLiteral = Literal["same", "unique", "first", "only"]


def resolve_merge_strategy(
    strategy: Union[str, Callable, None]
) -> Callable[[Collection[Mapping]], Mapping]:
    if not isinstance(strategy, Callable):
        strategy = MERGE_STRATEGIES[strategy]
    return strategy


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

        self.no_change = new_idx.equals(old_idx)

        new_pos = new_idx.get_indexer(old_idx)
        old_pos = np.arange(len(new_pos))

        mask = new_pos != -1

        self.new_pos = new_pos[mask]
        self.old_pos = old_pos[mask]

    def __call__(self, el, *, axis=1, fill_value=None):
        return self.apply(el, axis=axis, fill_value=fill_value)

    def apply(self, el, *, axis, fill_value=None):
        if self.no_change and (el.shape[axis] == len(self.old_idx)):
            return el
        if isinstance(el, pd.DataFrame):
            return self._apply_to_df(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, sparse.spmatrix):
            return self._apply_to_sparse(el, axis=axis, fill_value=fill_value)
        else:
            return self._apply_to_array(el, axis=axis, fill_value=fill_value)

    def _apply_to_df(self, el: pd.DataFrame, *, axis, fill_value=None):
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
            el, indexer, axis=axis, allow_fill=True, fill_value=fill_value
        )

    def _apply_to_sparse(self, el: spmatrix, *, axis, fill_value=None) -> spmatrix:
        if fill_value is None:
            fill_value = default_fill_value([el])
        if fill_value != 0:
            to_fill = self.new_idx.get_indexer(self.new_idx.difference(self.old_idx))
        else:
            to_fill = np.array([])

        # Fixing outer indexing for missing values
        if el.shape[1] == 0:
            if fill_value == 0:
                return sparse.coo_matrix((el.shape[0], len(self.new_idx)))
            else:
                return np.broadcast_to(fill_value, (el.shape[0], len(self.new_idx)))

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


def merge_indices(
    inds: Iterable[pd.Index], join: Literal["inner", "outer"]
) -> pd.Index:
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


def gen_inner_reindexers(els, new_index, axis: Literal[0, 1] = 0):
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
            gen_reindexer(pd.RangeIndex(max_col), pd.RangeIndex(n)) for n in orig_cols
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
            cur_reindexers = gen_outer_reindexers(els, ns, new_index=index, axis=axis)
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


def concat_pairwise_mapping(
    mappings: Collection[Mapping], shapes: Collection[int], join_keys=intersect_keys
):
    result = {}
    for k in join_keys(mappings):
        els = [
            m.get(k, sparse.csr_matrix((s, s), dtype=bool))
            for m, s in zip(mappings, shapes)
        ]
        result[k] = sparse.block_diag(els, format="csr")
    return result


def merge_dataframes(
    dfs: Iterable[pd.DataFrame], new_index, merge_strategy=merge_unique
) -> pd.DataFrame:
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


def dim_indices(adata, *, axis=None, dim=None) -> pd.Index:
    """Helper function to get adata.{dim}_names."""
    _, dim = _resolve_dim(axis=axis, dim=dim)
    return getattr(adata, f"{dim}_names")


def dim_size(adata, *, axis=None, dim=None) -> int:
    """Helper function to get adata.shape[dim]."""
    ax, _ = _resolve_dim(axis, dim)
    return adata.shape[ax]


def concat(
    adatas: Union[Collection[AnnData], "typing.Mapping[str, AnnData]"],
    *,
    axis: Literal[0, 1] = 0,
    join: Literal["inner", "outer"] = "inner",
    merge: Union[StrategiesLiteral, Callable, None] = None,
    uns_merge: Union[StrategiesLiteral, Callable, None] = None,
    label: Optional[str] = None,
    keys: Optional[Collection] = None,
    index_unique: Optional[str] = None,
    fill_value: Optional[Any] = None,
    pairwise: bool = False,
) -> AnnData:
    """Concatenates AnnData objects along an axis.

    See the :doc:`concatenation` section in the docs for a more in-depth description.

    .. warning::

        This function is marked as experimental for the `0.7` release series, and will
        supercede the :meth:`AnnData.concatenate() <anndata.AnnData.concatenate>` method
        in future releases.

    Params
    ------
    adatas
        The objects to be concatenated. If a Mapping is passed, keys are used for the `keys`
        argument and values are concatenated.
    axis
        Which axis to concatenate along.
    join
        How to align values when concatenating. If "outer", the union of the other axis
        is taken. If "inner", the intersection. See :doc:`concatenation` for more.
    merge
        How elements not aligned to the axis being concatenated along are selected.
        Currently implemented strategies include:

        * `None`: No elements are kept.
        * `"same"`: Elements that are the same in each of the objects.
        * `"unique"`: Elements for which there is only one possible value.
        * `"first"`: The first element seen at each from each position.
        * `"only"`: Elements that show up in only one of the objects.
    uns_merge
        How the elements of `.uns` are selected. Uses the same set of strategies as
        the `merge` argument, except applied recursively.
    label
        Column in axis annotation (i.e. `.obs` or `.var`) to place batch information in.
        If it's None, no column is added.
    keys
        Names for each object being added. These values are used for column values for
        `label` or appended to the index if `index_unique` is not `None`. Defaults to
        incrementing integer labels.
    index_unique
        Whether to make the index unique by using the keys. If provided, this
        is the delimeter between "{orig_idx}{index_unique}{key}". When `None`,
        the original indices are kept.
    fill_value
        When `join="outer"`, this is the value that will be used to fill the introduced
        indices. By default, sparse arrays are padded with zeros, while dense arrays and
        DataFrames are padded with missing values.
    pairwise
        Whether pairwise elements along the concatenated dimension should be included.
        This is False by default, since the resulting arrays are often not meaningful.

    Notes
    -----

    .. warning::

        If you use `join='outer'` this fills 0s for sparse data when
        variables are absent in a batch. Use this with care. Dense data is
        filled with `NaN`.

    Examples
    --------

    Preparing example objects

    >>> import anndata as ad, pandas as pd, numpy as np
    >>> from scipy import sparse
    >>> a = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[0, 1], [2, 3]])),
    ...     obs=pd.DataFrame({"group": ["a", "b"]}, index=["s1", "s2"]),
    ...     var=pd.DataFrame(index=["var1", "var2"]),
    ...     varm={"ones": np.ones((2, 5)), "rand": np.random.randn(2, 3), "zeros": np.zeros((2, 5))},
    ...     uns={"a": 1, "b": 2, "c": {"c.a": 3, "c.b": 4}},
    ... )
    >>> b = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[4, 5, 6], [7, 8, 9]])),
    ...     obs=pd.DataFrame({"group": ["b", "c"], "measure": [1.2, 4.3]}, index=["s3", "s4"]),
    ...     var=pd.DataFrame(index=["var1", "var2", "var3"]),
    ...     varm={"ones": np.ones((3, 5)), "rand": np.random.randn(3, 5)},
    ...     uns={"a": 1, "b": 3, "c": {"c.b": 4}},
    ... )
    >>> c = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[10, 11], [12, 13]])),
    ...     obs=pd.DataFrame({"group": ["a", "b"]}, index=["s1", "s2"]),
    ...     var=pd.DataFrame(index=["var3", "var4"]),
    ...     uns={"a": 1, "b": 4, "c": {"c.a": 3, "c.b": 4, "c.c": 5}},
    ... )

    Concatenating along different axes

    >>> ad.concat([a, b]).to_df()
        var1  var2
    s1   0.0   1.0
    s2   2.0   3.0
    s3   4.0   5.0
    s4   7.0   8.0
    >>> ad.concat([a, c], axis=1).to_df()
        var1  var2  var3  var4
    s1   0.0   1.0  10.0  11.0
    s2   2.0   3.0  12.0  13.0

    Inner and outer joins

    >>> inner = ad.concat([a, b])  # Joining on intersection of variables
    >>> inner
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
    >>> (inner.obs_names, inner.var_names)  # doctest: +NORMALIZE_WHITESPACE
    (Index(['s1', 's2', 's3', 's4'], dtype='object'),
    Index(['var1', 'var2'], dtype='object'))
    >>> outer = ad.concat([a, b], join="outer") # Joining on union of variables
    >>> outer
    AnnData object with n_obs × n_vars = 4 × 3
        obs: 'group', 'measure'
    >>> outer.var_names
    Index(['var1', 'var2', 'var3'], dtype='object')
    >>> outer.to_df()  # Sparse arrays are padded with zeroes by default
        var1  var2  var3
    s1   0.0   1.0   0.0
    s2   2.0   3.0   0.0
    s3   4.0   5.0   6.0
    s4   7.0   8.0   9.0

    Keeping track of source objects

    >>> ad.concat({"a": a, "b": b}, label="batch").obs
       group batch
    s1     a     a
    s2     b     a
    s3     b     b
    s4     c     b
    >>> ad.concat([a, b], label="batch", keys=["a", "b"]).obs  # Equivalent to previous
       group batch
    s1     a     a
    s2     b     a
    s3     b     b
    s4     c     b
    >>> ad.concat({"a": a, "b": b}, index_unique="-").obs
         group
    s1-a     a
    s2-a     b
    s3-b     b
    s4-b     c

    Combining values not aligned to axis of concatenation

    >>> ad.concat([a, b], merge="same")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones'
    >>> ad.concat([a, b], merge="unique")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones', 'zeros'
    >>> ad.concat([a, b], merge="first")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'ones', 'rand', 'zeros'
    >>> ad.concat([a, b], merge="only")
    AnnData object with n_obs × n_vars = 4 × 2
        obs: 'group'
        varm: 'zeros'

    The same merge strategies can be used for elements in `.uns`

    >>> dict(ad.concat([a, b, c], uns_merge="same").uns)
    {'a': 1, 'c': {'c.b': 4}}
    >>> dict(ad.concat([a, b, c], uns_merge="unique").uns)
    {'a': 1, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}
    >>> dict(ad.concat([a, b, c], uns_merge="only").uns)
    {'c': {'c.c': 5}}
    >>> dict(ad.concat([a, b, c], uns_merge="first").uns)
    {'a': 1, 'b': 2, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}
    """
    # Argument normalization
    merge = resolve_merge_strategy(merge)
    uns_merge = resolve_merge_strategy(uns_merge)

    if isinstance(adatas, Mapping):
        if keys is not None:
            raise TypeError(
                "Cannot specify categories in both mapping keys and using `keys`. "
                "Only specify this once."
            )
        keys, adatas = list(adatas.keys()), list(adatas.values())
    else:
        adatas = list(adatas)

    if keys is None:
        keys = np.arange(len(adatas)).astype(str)

    axis, dim = _resolve_dim(axis=axis)
    alt_axis, alt_dim = _resolve_dim(axis=1 - axis)

    # Label column
    label_col = pd.Categorical.from_codes(
        np.repeat(np.arange(len(adatas)), [a.shape[axis] for a in adatas]),
        categories=keys,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(dim_indices(a, axis=axis)) for a in adatas], ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(label_col.map(str), sep=index_unique)
    concat_indices = pd.Index(concat_indices)

    alt_indices = merge_indices(
        [dim_indices(a, axis=alt_axis) for a in adatas], join=join
    )
    reindexers = [
        gen_reindexer(alt_indices, dim_indices(a, axis=alt_axis)) for a in adatas
    ]

    # Annotation for concatenation axis
    concat_annot = pd.concat(
        [getattr(a, dim) for a in adatas], join=join, ignore_index=True
    )
    concat_annot.index = concat_indices
    if label is not None:
        concat_annot[label] = label_col

    # Annotation for other axis
    alt_annot = merge_dataframes(
        [getattr(a, alt_dim) for a in adatas], alt_indices, merge
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
        if pairwise:
            concat_pairwise = concat_pairwise_mapping(
                mappings=[getattr(a, f"{dim}p") for a in adatas],
                shapes=[a.shape[axis] for a in adatas],
                join_keys=intersect_keys,
            )
        else:
            concat_pairwise = {}
    elif join == "outer":
        layers = outer_concat_aligned_mapping(
            [a.layers for a in adatas], reindexers, axis=axis, fill_value=fill_value
        )
        concat_mapping = outer_concat_aligned_mapping(
            [getattr(a, f"{dim}m") for a in adatas],
            index=concat_indices,
            fill_value=fill_value,
        )
        if pairwise:
            concat_pairwise = concat_pairwise_mapping(
                mappings=[getattr(a, f"{dim}p") for a in adatas],
                shapes=[a.shape[axis] for a in adatas],
                join_keys=union_keys,
            )
        else:
            concat_pairwise = {}

    # TODO: Reindex lazily, so we don't have to make those copies until we're sure we need the element
    alt_mapping = merge(
        [
            {k: r(v, axis=0) for k, v in getattr(a, f"{alt_dim}m").items()}
            for r, a in zip(reindexers, adatas)
        ],
    )
    alt_pairwise = merge(
        [
            {k: r(r(v, axis=0), axis=1) for k, v in getattr(a, f"{alt_dim}p").items()}
            for r, a in zip(reindexers, adatas)
        ]
    )
    uns = uns_merge([a.uns for a in adatas])

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
            label=label,
            keys=keys,
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
            f"{dim}p": concat_pairwise,
            f"{alt_dim}p": alt_pairwise,
            "uns": uns,
            "raw": raw,
        }
    )
