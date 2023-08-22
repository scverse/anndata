"""
Code for merging/ concatenating AnnData objects.
"""
from __future__ import annotations

from collections import OrderedDict
from collections.abc import (
    Callable,
    Collection,
    Mapping,
    MutableSet,
    Iterable,
    Sequence,
)
from functools import reduce, singledispatch
from itertools import repeat
from operator import and_, or_, sub
from typing import Any, Optional, TypeVar, Union, Literal
import typing
from warnings import warn, filterwarnings

from natsort import natsorted
import numpy as np
import pandas as pd
from pandas.api.extensions import ExtensionDtype
from scipy import sparse
from scipy.sparse import spmatrix

from .anndata import AnnData
from ..compat import AwkArray, DaskArray, CupySparseMatrix, CupyArray, CupyCSRMatrix
from ..utils import asarray, dim_len
from .index import _subset, make_slice
from anndata._warnings import ExperimentalFeatureWarning

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


@equal.register(DaskArray)
def equal_dask_array(a, b) -> bool:
    import dask.array as da
    from dask.base import tokenize

    if a is b:
        return True
    if a.shape != b.shape:
        return False
    if isinstance(b, DaskArray):
        if tokenize(a) == tokenize(b):
            return True
    return da.equal(a, b, where=~(da.isnan(a) == da.isnan(b))).all()


@equal.register(np.ndarray)
def equal_array(a, b) -> bool:
    return equal(pd.DataFrame(a), pd.DataFrame(asarray(b)))


@equal.register(CupyArray)
def equal_cupyarray(a, b) -> bool:
    import cupy as cp

    return bool(cp.array_equal(a, b, equal_nan=True))


@equal.register(pd.Series)
def equal_series(a, b) -> bool:
    return a.equals(b)


@equal.register(sparse.spmatrix)
@equal.register(CupySparseMatrix)
def equal_sparse(a, b) -> bool:
    # It's a weird api, don't blame me
    import array_api_compat

    xp = array_api_compat.array_namespace(a.data)

    if isinstance(b, (CupySparseMatrix, sparse.spmatrix)):
        if isinstance(a, CupySparseMatrix):
            # Comparison broken for CSC matrices
            # https://github.com/cupy/cupy/issues/7757
            a, b = CupyCSRMatrix(a), CupyCSRMatrix(b)
        comp = a != b
        if isinstance(comp, bool):
            return not comp
        if isinstance(comp, CupySparseMatrix):
            # https://github.com/cupy/cupy/issues/7751
            comp = comp.get()
        # fmt: off
        return (
            (len(comp.data) == 0)
            or (
                xp.isnan(a[comp]).all()
                and xp.isnan(b[comp]).all()
            )
        )
        # fmt: on
    else:
        return False


@equal.register(AwkArray)
def equal_awkward(a, b) -> bool:
    from ..compat import awkward as ak

    return ak.almost_equal(a, b)


def as_sparse(x):
    if not isinstance(x, sparse.spmatrix):
        return sparse.csr_matrix(x)
    else:
        return x


def as_cp_sparse(x) -> CupySparseMatrix:
    import cupyx.scipy.sparse as cpsparse

    if isinstance(x, cpsparse.spmatrix):
        return x
    elif isinstance(x, np.ndarray):
        return cpsparse.csr_matrix(as_sparse(x))
    else:
        return cpsparse.csr_matrix(x)


def unify_dtypes(dfs: Iterable[pd.DataFrame]) -> list[pd.DataFrame]:
    """
    Attempts to unify datatypes from multiple dataframes.

    For catching cases where pandas would convert to object dtype.
    """
    # Get shared categorical columns
    df_dtypes = [dict(df.dtypes) for df in dfs]
    columns = reduce(lambda x, y: x.union(y), [df.columns for df in dfs])

    dtypes: dict[str, list[np.dtype | ExtensionDtype]] = {col: [] for col in columns}
    for col in columns:
        for df in df_dtypes:
            dtypes[col].append(df.get(col, None))

    if len(dtypes) == 0:
        return dfs
    else:
        dfs = [df.copy(deep=False) for df in dfs]

    new_dtypes = {}
    for col in dtypes.keys():
        target_dtype = try_unifying_dtype(dtypes[col])
        if target_dtype is not None:
            new_dtypes[col] = target_dtype

    for df in dfs:
        for col, dtype in new_dtypes.items():
            if col in df:
                df[col] = df[col].astype(dtype)

    return dfs


def try_unifying_dtype(
    col: Sequence[np.dtype | ExtensionDtype],
) -> pd.core.dtypes.base.ExtensionDtype | None:
    """
    If dtypes can be unified, returns the dtype they would be unified to.

    Returns None if they can't be unified, or if we can expect pandas to unify them for
    us.

    Params
    ------
    col:
        A list of dtypes to unify. Can be numpy/ pandas dtypes, or None (which denotes
        a missing value)
    """
    dtypes: set[pd.CategoricalDtype] = set()
    # Categorical
    if any(isinstance(dtype, pd.CategoricalDtype) for dtype in col):
        ordered = False
        for dtype in col:
            if isinstance(dtype, pd.CategoricalDtype):
                dtypes.add(dtype)
                ordered = ordered | dtype.ordered
            elif not pd.isnull(dtype):
                return False
        if len(dtypes) > 0 and not ordered:
            categories = reduce(
                lambda x, y: x.union(y),
                [dtype.categories for dtype in dtypes if not pd.isnull(dtype)],
            )

            return pd.CategoricalDtype(natsorted(categories), ordered=False)
    # Boolean
    elif all(pd.api.types.is_bool_dtype(dtype) or dtype is None for dtype in col):
        if any(dtype is None for dtype in col):
            return pd.BooleanDtype()
        else:
            return None
    else:
        return None


def check_combinable_cols(cols: list[pd.Index], join: Literal["inner", "outer"]):
    """Given columns for a set of dataframes, checks if the can be combined.

    Looks for if there are duplicated column names that would show up in the result.
    """
    repeated_cols = reduce(lambda x, y: x.union(y[y.duplicated()]), cols, set())
    if join == "inner":
        intersecting_cols = intersect_keys(cols)
        problem_cols = repeated_cols.intersection(intersecting_cols)
    elif join == "outer":
        problem_cols = repeated_cols
    else:
        raise ValueError()

    if len(problem_cols) > 0:
        problem_cols = list(problem_cols)
        raise pd.errors.InvalidIndexError(
            f"Cannot combine dataframes as some contained duplicated column names - "
            "causing ambiguity.\n\n"
            f"The problem columns are: {problem_cols}"
        )


# TODO: open PR or feature request to cupy
def _cpblock_diag(mats, format=None, dtype=None):
    """
    Modified version of scipy.sparse.block_diag for cupy sparse.
    """
    import cupy as cp
    from cupyx.scipy import sparse as cpsparse

    row = []
    col = []
    data = []
    r_idx = 0
    c_idx = 0
    for a in mats:
        # if isinstance(a, (list, numbers.Number)):
        #     a = cpsparse.coo_matrix(a)
        nrows, ncols = a.shape
        if cpsparse.issparse(a):
            a = a.tocoo()
            row.append(a.row + r_idx)
            col.append(a.col + c_idx)
            data.append(a.data)
        else:
            a_row, a_col = cp.divmod(cp.arange(nrows * ncols), ncols)
            row.append(a_row + r_idx)
            col.append(a_col + c_idx)
            data.append(a.reshape(-1))
        r_idx += nrows
        c_idx += ncols
    row = cp.concatenate(row)
    col = cp.concatenate(col)
    data = cp.concatenate(data)
    return cpsparse.coo_matrix(
        (data, (row, col)), shape=(r_idx, c_idx), dtype=dtype
    ).asformat(format)


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


class Reindexer:
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
        """
        Reindex element so el[axis] is aligned to self.new_idx.

        Missing values are to be replaced with `fill_value`.
        """
        if self.no_change and (dim_len(el, axis) == len(self.old_idx)):
            return el
        if isinstance(el, pd.DataFrame):
            return self._apply_to_df(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, sparse.spmatrix):
            return self._apply_to_sparse(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, AwkArray):
            return self._apply_to_awkward(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, DaskArray):
            return self._apply_to_dask_array(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, CupyArray):
            return self._apply_to_cupy_array(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, CupySparseMatrix):
            return self._apply_to_sparse(el, axis=axis, fill_value=fill_value)
        else:
            return self._apply_to_array(el, axis=axis, fill_value=fill_value)

    def _apply_to_df(self, el: pd.DataFrame, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = np.NaN
        return el.reindex(self.new_idx, axis=axis, fill_value=fill_value)

    def _apply_to_dask_array(self, el: DaskArray, *, axis, fill_value=None):
        import dask.array as da

        if fill_value is None:
            fill_value = default_fill_value([el])
        shape = list(el.shape)
        if el.shape[axis] == 0:
            # Presumably faster since it won't allocate the full array
            shape[axis] = len(self.new_idx)
            return da.broadcast_to(fill_value, tuple(shape))

        indexer = self.old_idx.get_indexer(self.new_idx)
        sub_el = _subset(el, make_slice(indexer, axis, len(shape)))

        if any(indexer == -1):
            sub_el[make_slice(indexer == -1, axis, len(shape))] = fill_value

        return sub_el

    def _apply_to_cupy_array(self, el, *, axis, fill_value=None):
        import cupy as cp

        if fill_value is None:
            fill_value = default_fill_value([el])
        if el.shape[axis] == 0:
            # Presumably faster since it won't allocate the full array
            shape = list(el.shape)
            shape[axis] = len(self.new_idx)
            return cp.broadcast_to(cp.asarray(fill_value), tuple(shape))

        old_idx_tuple = [slice(None)] * len(el.shape)
        old_idx_tuple[axis] = self.old_pos
        old_idx_tuple = tuple(old_idx_tuple)
        new_idx_tuple = [slice(None)] * len(el.shape)
        new_idx_tuple[axis] = self.new_pos
        new_idx_tuple = tuple(new_idx_tuple)

        out_shape = list(el.shape)
        out_shape[axis] = len(self.new_idx)

        out = cp.full(tuple(out_shape), fill_value)
        out[new_idx_tuple] = el[old_idx_tuple]

        return out

    def _apply_to_array(self, el, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = default_fill_value([el])
        if el.shape[axis] == 0:
            # Presumably faster since it won't allocate the full array
            shape = list(el.shape)
            shape[axis] = len(self.new_idx)
            return np.broadcast_to(fill_value, tuple(shape))

        indexer = self.old_idx.get_indexer(self.new_idx)

        # Indexes real fast, and does outer indexing
        return pd.api.extensions.take(
            el, indexer, axis=axis, allow_fill=True, fill_value=fill_value
        )

    def _apply_to_sparse(self, el: spmatrix, *, axis, fill_value=None) -> spmatrix:
        if isinstance(el, CupySparseMatrix):
            from cupyx.scipy import sparse
        else:
            from scipy import sparse
        import array_api_compat

        xp = array_api_compat.array_namespace(el.data)

        if fill_value is None:
            fill_value = default_fill_value([el])
        if fill_value != 0:
            to_fill = self.new_idx.get_indexer(self.new_idx.difference(self.old_idx))
        else:
            to_fill = xp.array([])

        # Fixing outer indexing for missing values
        if el.shape[axis] == 0:
            shape = list(el.shape)
            shape[axis] = len(self.new_idx)
            shape = tuple(shape)
            if fill_value == 0:
                return sparse.csr_matrix(shape)
            else:
                return type(el)(xp.broadcast_to(xp.asarray(fill_value), shape))

        fill_idxer = None

        if len(to_fill) > 0 or isinstance(el, CupySparseMatrix):
            idxmtx_dtype = xp.promote_types(el.dtype, xp.array(fill_value).dtype)
        else:
            idxmtx_dtype = bool

        if axis == 1:
            idxmtx = sparse.coo_matrix(
                (
                    xp.ones(len(self.new_pos), dtype=idxmtx_dtype),
                    (xp.asarray(self.old_pos), xp.asarray(self.new_pos)),
                ),
                shape=(len(self.old_idx), len(self.new_idx)),
                dtype=idxmtx_dtype,
            )
            out = el @ idxmtx

            if len(to_fill) > 0:
                out = out.tocsc()
                fill_idxer = (slice(None), to_fill)
        elif axis == 0:
            idxmtx = sparse.coo_matrix(
                (
                    xp.ones(len(self.new_pos), dtype=idxmtx_dtype),
                    (xp.asarray(self.new_pos), xp.asarray(self.old_pos)),
                ),
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

    def _apply_to_awkward(self, el: AwkArray, *, axis, fill_value=None):
        import awkward as ak

        if self.no_change:
            return el
        elif axis == 1:  # Indexing by field
            if self.new_idx.isin(self.old_idx).all():  # inner join
                return el[self.new_idx]
            else:  # outer join
                # TODO: this code isn't actually hit, we should refactor
                raise Exception("This should be unreachable, please open an issue.")
        else:
            if len(self.new_idx) > len(self.old_idx):
                el = ak.pad_none(el, 1, axis=axis)  # axis == 0
            return el[self.old_idx.get_indexer(self.new_idx)]


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
           [1., 0., 0.]])
    """
    return Reindexer(cur_var, new_var)


def np_bool_to_pd_bool_array(df: pd.DataFrame):
    for col_name, col_type in dict(df.dtypes).items():
        if col_type is np.dtype(bool):
            df[col_name] = pd.array(df[col_name].values)
    return df


def concat_arrays(arrays, reindexers, axis=0, index=None, fill_value=None):
    arrays = list(arrays)
    if fill_value is None:
        fill_value = default_fill_value(arrays)

    if any(isinstance(a, pd.DataFrame) for a in arrays):
        # TODO: This is hacky, 0 is a sentinel for outer_concat_aligned_mapping
        if not all(
            isinstance(a, pd.DataFrame) or a is MissingVal or 0 in a.shape
            for a in arrays
        ):
            raise NotImplementedError(
                "Cannot concatenate a dataframe with other array types."
            )
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat(
            unify_dtypes([f(x) for f, x in zip(reindexers, arrays)]),
            ignore_index=True,
            axis=axis,
        )
        df.index = index
        return df
    elif any(isinstance(a, AwkArray) for a in arrays):
        from ..compat import awkward as ak

        if not all(
            isinstance(a, AwkArray) or a is MissingVal or 0 in a.shape for a in arrays
        ):
            raise NotImplementedError(
                "Cannot concatenate an AwkwardArray with other array types."
            )

        return ak.concatenate([f(a) for f, a in zip(reindexers, arrays)], axis=axis)
    elif any(isinstance(a, CupySparseMatrix) for a in arrays):
        import cupyx.scipy.sparse as cpsparse

        sparse_stack = (cpsparse.vstack, cpsparse.hstack)[axis]
        return sparse_stack(
            [
                f(as_cp_sparse(a), axis=1 - axis, fill_value=fill_value)
                for f, a in zip(reindexers, arrays)
            ],
            format="csr",
        )
    elif any(isinstance(a, CupyArray) for a in arrays):
        import cupy as cp

        if not all(isinstance(a, CupyArray) or 0 in a.shape for a in arrays):
            raise NotImplementedError(
                "Cannot concatenate a cupy array with other array types."
            )
        return cp.concatenate(
            [
                f(cp.asarray(x), fill_value=fill_value, axis=1 - axis)
                for f, x in zip(reindexers, arrays)
            ],
            axis=axis,
        )
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
    elif any(isinstance(el, AwkArray) for el in els if not_missing(el)):
        if not all(isinstance(el, AwkArray) for el in els if not_missing(el)):
            raise NotImplementedError(
                "Cannot concatenate an AwkwardArray with other array types."
            )
        common_keys = intersect_keys(el.fields for el in els)
        reindexers = [
            Reindexer(pd.Index(el.fields), pd.Index(list(common_keys))) for el in els
        ]
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
            (lambda x: x)
            if not_missing(el)
            else (lambda _, shape=shape: pd.DataFrame(index=range(shape)))
            for el, shape in zip(els, shapes)
        ]
    elif any(isinstance(el, AwkArray) for el in els if not_missing(el)):
        import awkward as ak

        if not all(isinstance(el, AwkArray) for el in els if not_missing(el)):
            raise NotImplementedError(
                "Cannot concatenate an AwkwardArray with other array types."
            )
        warn(
            "Outer joins on awkward.Arrays will have different return values in the future."
            "For details, and to offer input, please see:\n\n\t"
            "https://github.com/scverse/anndata/issues/898",
            ExperimentalFeatureWarning,
        )
        filterwarnings(
            "ignore",
            category=ExperimentalFeatureWarning,
            message=r"Outer joins on awkward.Arrays will have different return values.*",
        )
        # all_keys = union_keys(el.fields for el in els if not_missing(el))
        reindexers = []
        for el in els:
            if not_missing(el):
                reindexers.append(lambda x: x)
            else:
                reindexers.append(
                    lambda x: ak.pad_none(
                        ak.Array([]),
                        len(x),
                        0,
                    )
                )
    else:
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

        # Handling of missing values here is hacky for dataframes
        # We should probably just handle missing elements for all types
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
        if all(isinstance(el, (CupySparseMatrix, CupyArray)) for el in els):
            result[k] = _cpblock_diag(els, format="csr")
        else:
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


def _resolve_dim(*, dim: str = None, axis: int = None) -> tuple[int, str]:
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


# TODO: Resolve https://github.com/scverse/anndata/issues/678 and remove this function
def concat_Xs(adatas, reindexers, axis, fill_value):
    """
    Shimy until support for some missing X's is implemented.

    Basically just checks if it's one of the two supported cases, or throws an error.

    This is not done inline in `concat` because we don't want to maintain references
    to the values of a.X.
    """
    Xs = [a.X for a in adatas]
    if all(X is None for X in Xs):
        return None
    elif any(X is None for X in Xs):
        raise NotImplementedError(
            "Some (but not all) of the AnnData's to be concatenated had no .X value. "
            "Concatenation is currently only implemented for cases where all or none of"
            " the AnnData's have .X assigned."
        )
    else:
        return concat_arrays(Xs, reindexers, axis=axis, fill_value=fill_value)


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

    See the :doc:`concatenation <../concatenation>` section in the docs for a more in-depth description.

    Params
    ------
    adatas
        The objects to be concatenated. If a Mapping is passed, keys are used for the `keys`
        argument and values are concatenated.
    axis
        Which axis to concatenate along.
    join
        How to align values when concatenating. If "outer", the union of the other axis
        is taken. If "inner", the intersection. See :doc:`concatenation <../concatenation>`
        for more.
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
    s1     0     1
    s2     2     3
    s3     4     5
    s4     7     8
    >>> ad.concat([a, c], axis=1).to_df()
        var1  var2  var3  var4
    s1     0     1    10    11
    s2     2     3    12    13

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
    s1     0     1     0
    s2     2     3     0
    s3     4     5     6
    s4     7     8     9

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
    check_combinable_cols([getattr(a, dim).columns for a in adatas], join=join)
    concat_annot = pd.concat(
        unify_dtypes([getattr(a, dim) for a in adatas]),
        join=join,
        ignore_index=True,
    )
    concat_annot.index = concat_indices
    if label is not None:
        concat_annot[label] = label_col

    # Annotation for other axis
    alt_annot = merge_dataframes(
        [getattr(a, alt_dim) for a in adatas], alt_indices, merge
    )

    X = concat_Xs(adatas, reindexers, axis=axis, fill_value=fill_value)

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
