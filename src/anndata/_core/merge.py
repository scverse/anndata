"""
Code for merging/ concatenating AnnData objects.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable, Mapping, MutableSet
from functools import partial, reduce, singledispatch
from itertools import repeat
from operator import and_, or_, sub
from typing import TYPE_CHECKING, Literal, TypeVar
from warnings import warn

import numpy as np
import pandas as pd
import scipy
from natsort import natsorted
from packaging.version import Version
from scipy import sparse

from anndata._core.file_backing import to_memory
from anndata._warnings import ExperimentalFeatureWarning

from ..compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    CupyArray,
    CupyCSRMatrix,
    CupySparseMatrix,
    DaskArray,
    _map_cat_to_str,
)
from ..utils import asarray, axis_len, warn_once
from .anndata import AnnData
from .index import _subset, make_slice
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Collection, Generator, Iterable, Sequence
    from typing import Any

    from pandas.api.extensions import ExtensionDtype

    from anndata._types import Join_T

    from ..compat import XDataArray, XDataset

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

    def union(self, *vals) -> OrderedSet:
        return reduce(or_, vals, self)

    def discard(self, val):
        if val in self:
            del self.dict[val]

    def difference(self, *vals) -> OrderedSet:
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
# Unfortunately equality of arrays is poorly defined.
# * `np.array_equal` does not work for sparse arrays
# * `np.array_equal(..., equal_nan=True)` does not work for null values at the moment
#   (see https://github.com/numpy/numpy/issues/16377)
# So we have to define it ourselves with these two issues in mind.
# TODO: Hopefully this will stop being an issue in the future and this code can be removed.
@singledispatch
def equal(a, b) -> bool:
    a = asarray(a)
    b = asarray(b)
    if a.ndim == b.ndim == 0:
        return bool(a == b)
    a_na = (
        pd.isna(a) if a.dtype.names is None else np.False_
    )  # pd.isna doesn't work for record arrays
    b_na = pd.isna(b) if b.dtype.names is None else np.False_
    return np.array_equal(a_na, b_na) and np.array_equal(a[~a_na], b[~b_na])


@equal.register(pd.DataFrame)
@equal.register(Dataset2D)
@equal.register(pd.Series)
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
    if isinstance(b, DaskArray) and tokenize(a) == tokenize(b):
        return True
    if isinstance(a._meta, CSMatrix):
        # TODO: Maybe also do this in the other case?
        return da.map_blocks(equal, a, b, drop_axis=(0, 1)).all()
    else:
        return da.equal(a, b, where=~(da.isnan(a) == da.isnan(b))).all()


@equal.register(np.ndarray)
def equal_array(a, b) -> bool:
    # Reshaping allows us to compare inputs with >2 dimensions
    # We cast to pandas since it will still work with non-numeric types
    b = asarray(b)
    if a.shape != b.shape:
        return False

    return equal(pd.DataFrame(a.reshape(-1)), pd.DataFrame(b.reshape(-1)))


@equal.register(CupyArray)
def equal_cupyarray(a, b) -> bool:
    import cupy as cp

    return bool(cp.array_equal(a, b, equal_nan=True))


@equal.register(CSMatrix)
@equal.register(CSArray)
@equal.register(CupySparseMatrix)
def equal_sparse(a, b) -> bool:
    # It's a weird api, don't blame me
    import array_api_compat

    xp = array_api_compat.array_namespace(a.data)

    if isinstance(b, CupySparseMatrix | CSMatrix | CSArray):
        if isinstance(a, CupySparseMatrix):
            # Comparison broken for CSC matrices
            # https://github.com/cupy/cupy/issues/7757
            a, b = CupyCSRMatrix(a), CupyCSRMatrix(b)
        if Version(scipy.__version__) >= Version("1.16.0rc1"):
            # TODO: https://github.com/scipy/scipy/issues/23068
            return bool(
                a.format == b.format
                and (a.shape == b.shape)
                and np.all(a.indptr == b.indptr)
                and np.all(a.indices == b.indices)
                and np.all((a.data == b.data) | (np.isnan(a.data) & np.isnan(b.data)))
            )
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


def as_sparse(x, *, use_sparse_array: bool = False) -> CSMatrix | CSArray:
    if not isinstance(x, CSMatrix | CSArray):
        in_memory_array_class = (
            sparse.csr_array if use_sparse_array else sparse.csr_matrix
        )
        if isinstance(x, DaskArray):
            x = x.map_blocks(
                sparse.csr_matrix,
                meta=sparse.csr_matrix(x._meta),
                dtype=x.dtype,
            ).compute()
        return in_memory_array_class(x)
    return x


def as_cp_sparse(x) -> CupySparseMatrix:
    import cupyx.scipy.sparse as cpsparse

    if isinstance(x, cpsparse.spmatrix):
        return x
    elif isinstance(x, np.ndarray):
        return cpsparse.csr_matrix(as_sparse(x))
    else:
        return cpsparse.csr_matrix(x)


def unify_dtypes(
    dfs: Iterable[pd.DataFrame | Dataset2D],
) -> list[pd.DataFrame | Dataset2D]:
    """
    Attempts to unify datatypes from multiple dataframes.

    For catching cases where pandas would convert to object dtype.
    """
    dfs = list(dfs)
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

    new_dtypes = {
        col: target_dtype
        for col, dtype in dtypes.items()
        if (target_dtype := try_unifying_dtype(dtype)) is not None
    }

    for df in dfs:
        for col, dtype in new_dtypes.items():
            if col in df:
                df[col] = df[col].astype(dtype)

    return dfs


def try_unifying_dtype(  # noqa PLR0911, PLR0912
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
                return None
        if len(dtypes) > 0:
            categories = reduce(
                lambda x, y: x.union(y),
                (dtype.categories for dtype in dtypes if not pd.isnull(dtype)),
            )

            if not ordered:
                return pd.CategoricalDtype(natsorted(categories), ordered=False)
            else:  # for xarray Datasets, see https://github.com/pydata/xarray/issues/10247
                categories_intersection = reduce(
                    lambda x, y: x.intersection(y),
                    (
                        dtype.categories
                        for dtype in dtypes
                        if not pd.isnull(dtype) and len(dtype.categories) > 0
                    ),
                )
                if len(categories_intersection) < len(categories):
                    return object
                else:
                    same_orders = all(
                        dtype.ordered
                        for dtype in dtypes
                        if not pd.isnull(dtype) and len(dtype.categories) > 0
                    )
                    same_orders &= all(
                        np.all(categories == dtype.categories)
                        for dtype in dtypes
                        if not pd.isnull(dtype) and len(dtype.categories) > 0
                    )
                    if same_orders:
                        return next(iter(dtypes))
                    return object
    # Boolean
    elif all(pd.api.types.is_bool_dtype(dtype) or dtype is None for dtype in col):
        if any(dtype is None for dtype in col):
            return pd.BooleanDtype()
        else:
            return None
    else:
        return None


def check_combinable_cols(cols: list[pd.Index], join: Join_T):
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
        msg = (
            f"Cannot combine dataframes as some contained duplicated column names - "
            "causing ambiguity.\n\n"
            f"The problem columns are: {problem_cols}"
        )
        raise pd.errors.InvalidIndexError(msg)


# TODO: open PR or feature request to cupy
def _cp_block_diag(mats, format=None, dtype=None):
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


def _dask_block_diag(mats):
    from itertools import permutations

    import dask.array as da

    blocks = np.zeros((len(mats), len(mats)), dtype=object)
    for i, j in permutations(range(len(mats)), 2):
        blocks[i, j] = da.from_array(
            sparse.csr_matrix((mats[i].shape[0], mats[j].shape[1]))
        )
    for i, x in enumerate(mats):
        if not isinstance(x._meta, sparse.csr_matrix):
            x = x.map_blocks(sparse.csr_matrix)
        blocks[i, i] = x

    return da.block(blocks.tolist())


###################
# Per element logic
###################


def unique_value(vals: Collection[T]) -> T | MissingVal:
    """
    Given a collection vals, returns the unique value (if one exists), otherwise
    returns MissingValue.
    """
    unique_val = vals[0]
    for v in vals[1:]:
        if not equal(v, unique_val):
            return MissingVal
    return unique_val


def first(vals: Collection[T]) -> T | MissingVal:
    """
    Given a collection of vals, return the first non-missing one.If they're all missing,
    return MissingVal.
    """
    for val in vals:
        if not_missing(val):
            return val
    return MissingVal


def only(vals: Collection[T]) -> T | MissingVal:
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
    elif all(isinstance(v, Mapping) and not isinstance(v, Dataset2D) for v in vals):
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
    strategy: str | Callable | None,
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

    def apply(self, el, *, axis, fill_value=None):  # noqa: PLR0911
        """
        Reindex element so el[axis] is aligned to self.new_idx.

        Missing values are to be replaced with `fill_value`.
        """
        if self.no_change and (axis_len(el, axis) == len(self.old_idx)):
            return el
        if isinstance(el, pd.DataFrame | Dataset2D):
            return self._apply_to_df_like(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, CSMatrix | CSArray | CupySparseMatrix):
            return self._apply_to_sparse(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, AwkArray):
            return self._apply_to_awkward(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, DaskArray):
            return self._apply_to_dask_array(el, axis=axis, fill_value=fill_value)
        elif isinstance(el, CupyArray):
            return self._apply_to_cupy_array(el, axis=axis, fill_value=fill_value)
        else:
            return self._apply_to_array(el, axis=axis, fill_value=fill_value)

    def _apply_to_df_like(self, el: pd.DataFrame | Dataset2D, *, axis, fill_value=None):
        if fill_value is None:
            fill_value = np.nan
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

        indexer = self.idx
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

        indexer = self.idx

        # Indexes real fast, and does outer indexing
        return pd.api.extensions.take(
            el, indexer, axis=axis, allow_fill=True, fill_value=fill_value
        )

    def _apply_to_sparse(  # noqa: PLR0912
        self, el: CSMatrix | CSArray, *, axis, fill_value=None
    ) -> CSMatrix:
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
                if isinstance(el, CSArray):
                    memory_class = sparse.csr_array
                else:
                    memory_class = sparse.csr_matrix
                return memory_class(shape)
            else:
                return type(el)(xp.broadcast_to(xp.asarray(fill_value), shape))

        fill_idxer = None

        if len(to_fill) > 0 or isinstance(el, CupySparseMatrix):
            idxmtx_dtype = xp.promote_types(el.dtype, xp.array(fill_value).dtype)
        else:
            idxmtx_dtype = bool
        if isinstance(el, CSArray):
            memory_class = sparse.coo_array
        else:
            memory_class = sparse.coo_matrix
        if axis == 1:
            idxmtx = memory_class(
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
            idxmtx = memory_class(
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
                msg = "This should be unreachable, please open an issue."
                raise Exception(msg)
        else:
            if len(self.new_idx) > len(self.old_idx):
                el = ak.pad_none(el, 1, axis=axis)  # axis == 0
            return el[self.idx]

    @property
    def idx(self):
        return self.old_idx.get_indexer(self.new_idx)


def merge_indices(inds: Iterable[pd.Index], join: Join_T) -> pd.Index:
    if join == "inner":
        return reduce(lambda x, y: x.intersection(y), inds)
    elif join == "outer":
        return reduce(lambda x, y: x.union(y), inds)
    else:
        msg = f"`join` must be one of 'inner' or 'outer', got {join!r}"
        raise ValueError(msg)


def default_fill_value(els):
    """Given some arrays, returns what the default fill value should be.

    This is largely due to backwards compat, and might not be the ideal solution.
    """
    if any(
        isinstance(el, CSMatrix | CSArray)
        or (isinstance(el, DaskArray) and isinstance(el._meta, CSMatrix | CSArray))
        for el in els
    ):
        return 0
    else:
        return np.nan


def gen_reindexer(new_var: pd.Index, cur_var: pd.Index):
    """
    Given a new set of var_names, and a current set, generates a function which will reindex
    a matrix to be aligned with the new set.

    Usage
    -----

    >>> a = AnnData(sparse.eye(3, format="csr"), var=pd.DataFrame(index=list("abc")))
    >>> b = AnnData(sparse.eye(2, format="csr"), var=pd.DataFrame(index=list("ba")))
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


def concat_arrays(  # noqa: PLR0911, PLR0912
    arrays, reindexers, axis=0, index=None, fill_value=None, *, force_lazy: bool = False
):
    from anndata.experimental.backed._compat import Dataset2D

    arrays = list(arrays)
    if fill_value is None:
        fill_value = default_fill_value(arrays)

    if any(isinstance(a, Dataset2D) for a in arrays):
        if any(isinstance(a, pd.DataFrame) for a in arrays):
            arrays = [to_memory(a) if isinstance(a, Dataset2D) else a for a in arrays]
        elif not all(isinstance(a, Dataset2D) for a in arrays):
            msg = f"Cannot concatenate a Dataset2D with other array types {[type(a) for a in arrays if not isinstance(a, Dataset2D)]}."
            raise ValueError(msg)
        else:
            return concat_dataset2d_on_annot_axis(
                arrays, join="outer", force_lazy=force_lazy
            )
    if any(isinstance(a, pd.DataFrame) for a in arrays):
        # TODO: This is hacky, 0 is a sentinel for outer_concat_aligned_mapping
        if not all(
            isinstance(a, pd.DataFrame) or a is MissingVal or 0 in a.shape
            for a in arrays
        ):
            msg = "Cannot concatenate a dataframe with other array types."
            raise NotImplementedError(msg)
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat(
            unify_dtypes(f(x) for f, x in zip(reindexers, arrays, strict=True)),
            axis=axis,
            ignore_index=True,
        )
        df.index = index
        return df
    elif any(isinstance(a, AwkArray) for a in arrays):
        from ..compat import awkward as ak

        if not all(
            isinstance(a, AwkArray) or a is MissingVal or 0 in a.shape for a in arrays
        ):
            msg = "Cannot concatenate an AwkwardArray with other array types."
            raise NotImplementedError(msg)

        return ak.concatenate(
            [f(a) for f, a in zip(reindexers, arrays, strict=True)], axis=axis
        )
    elif any(isinstance(a, CupySparseMatrix) for a in arrays):
        import cupyx.scipy.sparse as cpsparse

        if not all(
            isinstance(a, CupySparseMatrix | CupyArray) or 0 in a.shape for a in arrays
        ):
            msg = "Cannot concatenate a cupy array with other array types."
            raise NotImplementedError(msg)
        sparse_stack = (cpsparse.vstack, cpsparse.hstack)[axis]
        return sparse_stack(
            [
                f(as_cp_sparse(a), axis=1 - axis, fill_value=fill_value)
                for f, a in zip(reindexers, arrays, strict=True)
            ],
            format="csr",
        )
    elif any(isinstance(a, CupyArray) for a in arrays):
        import cupy as cp

        if not all(isinstance(a, CupyArray) or 0 in a.shape for a in arrays):
            msg = "Cannot concatenate a cupy array with other array types."
            raise NotImplementedError(msg)
        return cp.concatenate(
            [
                f(cp.asarray(x), fill_value=fill_value, axis=1 - axis)
                for f, x in zip(reindexers, arrays, strict=True)
            ],
            axis=axis,
        )
    elif any(isinstance(a, CSMatrix | CSArray) for a in arrays):
        sparse_stack = (sparse.vstack, sparse.hstack)[axis]
        use_sparse_array = any(issubclass(type(a), CSArray) for a in arrays)
        mat = sparse_stack(
            [
                f(
                    as_sparse(a, use_sparse_array=use_sparse_array),
                    axis=1 - axis,
                    fill_value=fill_value,
                )
                for f, a in zip(reindexers, arrays, strict=True)
            ],
            format="csr",
        )
        return mat
    else:
        return np.concatenate(
            [
                f(x, fill_value=fill_value, axis=1 - axis)
                for f, x in zip(reindexers, arrays, strict=True)
            ],
            axis=axis,
        )


def inner_concat_aligned_mapping(
    mappings,
    *,
    reindexers=None,
    index=None,
    axis=0,
    concat_axis=None,
    force_lazy: bool = False,
):
    if concat_axis is None:
        concat_axis = axis
    result = {}

    for k in intersect_keys(mappings):
        els = [m[k] for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_inner_reindexers(
                els, new_index=index, axis=concat_axis
            )
        else:
            cur_reindexers = reindexers

        result[k] = concat_arrays(
            els, cur_reindexers, index=index, axis=concat_axis, force_lazy=force_lazy
        )
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
            msg = "Cannot concatenate an AwkwardArray with other array types."
            raise NotImplementedError(msg)
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
            for el, shape in zip(els, shapes, strict=True)
        ]
    elif any(isinstance(el, AwkArray) for el in els if not_missing(el)):
        import awkward as ak

        if not all(isinstance(el, AwkArray) for el in els if not_missing(el)):
            msg = "Cannot concatenate an AwkwardArray with other array types."
            raise NotImplementedError(msg)
        warn_once(
            "Outer joins on awkward.Arrays will have different return values in the future. "
            "For details, and to offer input, please see:\n\n\t"
            "https://github.com/scverse/anndata/issues/898",
            ExperimentalFeatureWarning,
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


def missing_element(
    n: int,
    els: list[CSArray | CSMatrix | np.ndarray | DaskArray],
    axis: Literal[0, 1] = 0,
    fill_value: Any | None = None,
    off_axis_size: int = 0,
) -> np.ndarray | DaskArray:
    """Generates value to use when there is a missing element."""
    should_return_dask = any(isinstance(el, DaskArray) for el in els)
    # 0 sized array for in-memory prevents allocating unnecessary memory while preserving broadcasting.
    shape = (n, off_axis_size) if axis == 0 else (off_axis_size, n)
    if should_return_dask:
        import dask.array as da

        return da.full(
            shape, default_fill_value(els) if fill_value is None else fill_value
        )
    return np.zeros(shape, dtype=bool)


def outer_concat_aligned_mapping(
    mappings,
    *,
    reindexers=None,
    index=None,
    axis=0,
    concat_axis=None,
    fill_value=None,
    force_lazy: bool = False,
):
    if concat_axis is None:
        concat_axis = axis
    result = {}
    ns = [m.parent.shape[axis] for m in mappings]

    for k in union_keys(mappings):
        els = [m.get(k, MissingVal) for m in mappings]
        if reindexers is None:
            cur_reindexers = gen_outer_reindexers(
                els, ns, new_index=index, axis=concat_axis
            )
        else:
            cur_reindexers = reindexers

        # Dask needs to create a full array and can't do the size-0 trick
        off_axis_size = 0
        if any(isinstance(e, DaskArray) for e in els):
            if not isinstance(cur_reindexers[0], Reindexer):  # pragma: no cover
                msg = "Cannot re-index a dask array without a Reindexer"
                raise ValueError(msg)
            off_axis_size = cur_reindexers[0].idx.shape[0]
        # Handling of missing values here is hacky for dataframes
        # We should probably just handle missing elements for all types
        result[k] = concat_arrays(
            [
                el
                if not_missing(el)
                else missing_element(
                    n,
                    axis=concat_axis,
                    els=els,
                    fill_value=fill_value,
                    off_axis_size=off_axis_size,
                )
                for el, n in zip(els, ns, strict=True)
            ],
            cur_reindexers,
            axis=concat_axis,
            index=index,
            fill_value=fill_value,
            force_lazy=force_lazy,
        )
    return result


def concat_pairwise_mapping(
    mappings: Collection[Mapping], shapes: Collection[int], join_keys=intersect_keys
):
    result = {}
    if any(any(isinstance(v, CSArray) for v in m.values()) for m in mappings):
        sparse_class = sparse.csr_array
    else:
        sparse_class = sparse.csr_matrix

    for k in join_keys(mappings):
        els = [
            m.get(k, sparse_class((s, s), dtype=bool))
            for m, s in zip(mappings, shapes, strict=True)
        ]
        if all(isinstance(el, CupySparseMatrix | CupyArray) for el in els):
            result[k] = _cp_block_diag(els, format="csr")
        elif all(isinstance(el, DaskArray) for el in els):
            result[k] = _dask_block_diag(els)
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
        for b, m in zip(batch_keys, mappings, strict=True):
            val = m.get(key, None)
            if val is not None:
                out[f"{key}{join_index}{b}"] = val
    return out


def _resolve_axis(
    axis: Literal["obs", 0, "var", 1],
) -> tuple[Literal[0], Literal["obs"]] | tuple[Literal[1], Literal["var"]]:
    if axis in {0, "obs"}:
        return (0, "obs")
    if axis in {1, "var"}:
        return (1, "var")
    msg = f"`axis` must be either 0, 1, 'obs', or 'var', was {axis}"
    raise ValueError(msg)


def axis_indices(adata: AnnData, axis: Literal["obs", 0, "var", 1]) -> pd.Index:
    """Helper function to get adata.{dim}_names."""
    _, axis_name = _resolve_axis(axis)
    attr = getattr(adata, axis_name)
    if isinstance(attr, Dataset2D):
        return attr.true_index
    else:
        return attr.index


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
        msg = (
            "Some (but not all) of the AnnData's to be concatenated had no .X value. "
            "Concatenation is currently only implemented for cases where all or none of"
            " the AnnData's have .X assigned."
        )
        raise NotImplementedError(msg)
    else:
        return concat_arrays(Xs, reindexers, axis=axis, fill_value=fill_value)


def make_dask_col_from_extension_dtype(
    col: XDataArray, *, use_only_object_dtype: bool = False
) -> DaskArray:
    """
    Creates dask arrays from :class:`pandas.api.extensions.ExtensionArray` dtype :class:`xarray.DataArray`s.

    Parameters
    ----------
    col
        The columns to be converted
    use_only_object_dtype
        Whether or not to cast all :class:`pandas.api.extensions.ExtensionArray` dtypes to `object` type, by default False

    Returns
    -------
    A :class:`dask.Array`: representation of the column.
    """
    import dask.array as da

    from anndata._io.specs.lazy_methods import (
        compute_chunk_layout_for_axis_size,
        get_chunksize,
        maybe_open_h5,
    )
    from anndata.compat import XDataArray
    from anndata.compat import xarray as xr
    from anndata.experimental import read_elem_lazy

    base_path_or_zarr_group = col.attrs.get("base_path_or_zarr_group")
    elem_name = col.attrs.get("elem_name")
    if (
        base_path_or_zarr_group is not None and elem_name is not None
    ):  # lazy, backed by store
        dims = col.dims
        coords = col.coords.copy()
        with maybe_open_h5(base_path_or_zarr_group, elem_name) as f:
            maybe_chunk_size = get_chunksize(read_elem_lazy(f))
            chunk_size = (
                compute_chunk_layout_for_axis_size(
                    1000 if maybe_chunk_size is None else maybe_chunk_size[0],
                    col.shape[0],
                ),
            )

        def get_chunk(block_info=None):
            # reopening is important to get around h5py's unserializable lock in processes
            with maybe_open_h5(base_path_or_zarr_group, elem_name) as f:
                v = read_elem_lazy(f)
                variable = xr.Variable(
                    data=xr.core.indexing.LazilyIndexedArray(v), dims=dims
                )
                data_array = XDataArray(
                    variable,
                    coords=coords,
                    dims=dims,
                )
                idx = tuple(
                    slice(start, stop)
                    for start, stop in block_info[None]["array-location"]
                )
                chunk = np.array(data_array.data[idx])
            return chunk

        if col.dtype == "category" or col.dtype == "string" or use_only_object_dtype:  # noqa PLR1714
            dtype = "object"
        else:
            dtype = col.dtype.numpy_dtype
        return da.map_blocks(
            get_chunk,
            chunks=chunk_size,
            meta=np.array([], dtype=dtype),
            dtype=dtype,
        )

    return da.from_array(col.values, chunks=-1)  # in-memory


def make_xarray_extension_dtypes_dask(
    annotations: Iterable[Dataset2D], *, use_only_object_dtype: bool = False
) -> Generator[XDataset, None, None]:
    """
    Creates a generator of Dataset2D objects with dask arrays in place of :class:`pandas.api.extensions.ExtensionArray` dtype columns.

    Parameters
    ----------
    annotations
        The datasets to be altered
    use_only_object_dtype
        Whether or not to cast all :class:`pandas.api.extensions.ExtensionArray` dtypes to `object` type, by default False

    Yields
    ------
    An altered dataset.
    """
    for a in annotations:
        extension_cols = {
            col for col in a.columns if pd.api.types.is_extension_array_dtype(a[col])
        }

        yield a.copy(
            data={
                name: (
                    make_dask_col_from_extension_dtype(
                        col, use_only_object_dtype=use_only_object_dtype
                    )
                    if name in extension_cols
                    else col
                )
                for name, col in a._items()
            }
        )


DS_CONCAT_DUMMY_INDEX_NAME = "concat_index"


def concat_dataset2d_on_annot_axis(
    annotations: Iterable[Dataset2D],
    join: Join_T,
    *,
    force_lazy: bool,
    concat_indices: pd.Index | None = None,
) -> Dataset2D:
    """Create a concatenate dataset from a list of :class:`~anndata.experimental.backed.Dataset2D` objects.
    The goal of this function is to mimic `pd.concat(..., ignore_index=True)` so has some complicated logic
    for handling the "index" to ensure (a) nothing is loaded into memory and (b) the true index is always tracked.

    Parameters
    ----------
    annotations
        The :class:`~anndata.experimental.backed.Dataset2D` objects to be concatenated.
    join
        Type of join operation
    force_lazy
        Whether to lazily concatenate elements using dask even when eager concatenation is possible.
    concat_indices
        Already calculated indices to be used as the index on the concatenated object.

    Returns
    -------
    Concatenated :class:`~anndata.experimental.backed.Dataset2D`
    """
    from anndata._core.xarray import Dataset2D
    from anndata._io.specs.lazy_methods import DUMMY_RANGE_INDEX_KEY
    from anndata.compat import xarray as xr

    annotations_re_indexed = []
    have_backed = any(a.is_backed for a in annotations)
    if have_backed or force_lazy:
        annotations = make_xarray_extension_dtypes_dask(annotations)
    else:
        annotations = unify_dtypes(annotations)
    for a in annotations:
        old_key = a.index_dim
        is_fake_index = old_key != a.true_index_dim
        # First create a dummy index
        a.ds.coords[DS_CONCAT_DUMMY_INDEX_NAME] = (
            old_key,
            pd.RangeIndex(a.shape[0]),
        )
        # Set all the dimensions to this new dummy index
        ds_swapped = a.ds.swap_dims({old_key: DS_CONCAT_DUMMY_INDEX_NAME})
        # Move the old coordinate into a variable
        old_coord = ds_swapped.coords[old_key]
        del ds_swapped.coords[old_key]
        ds_swapped[old_key] = old_coord
        a = Dataset2D(ds_swapped)
        if not is_fake_index:
            a.true_index_dim = old_key
        annotations_re_indexed.append(a)
    # Concat along the dummy index
    ds_concat = xr.concat(
        [a.ds for a in annotations_re_indexed],
        join=join,
        dim=DS_CONCAT_DUMMY_INDEX_NAME,
    )
    ds_concat.attrs.pop("indexing_key", None)
    # Wrapping allows us to use the Dataset2D methods
    # directly for setting certain attrs/coords without duplicating here.
    ds_concat_2d = Dataset2D(ds_concat)
    ds_concat_2d.is_backed = have_backed
    if concat_indices is not None:
        concat_indices.name = DS_CONCAT_DUMMY_INDEX_NAME
        ds_concat_2d.index = concat_indices
        ds_concat = ds_concat_2d.ds
    else:
        ds_concat.coords[DS_CONCAT_DUMMY_INDEX_NAME] = pd.RangeIndex(
            ds_concat.coords[DS_CONCAT_DUMMY_INDEX_NAME].shape[0]
        )
    # Drop any lingering dimensions (swap doesn't delete)
    ds_concat = ds_concat.drop_dims(
        d for d in ds_concat.dims if d != DS_CONCAT_DUMMY_INDEX_NAME
    )
    # Create a new true index and then delete the columns resulting from the concatenation for each index.
    # This includes the dummy column (which is neither a dimension nor a true indexing column)
    if concat_indices is None:
        index = xr.concat(
            [a.true_xr_index for a in annotations_re_indexed],
            dim=DS_CONCAT_DUMMY_INDEX_NAME,
        )
        # prevent duplicate values
        index.coords[DS_CONCAT_DUMMY_INDEX_NAME] = ds_concat.coords[
            DS_CONCAT_DUMMY_INDEX_NAME
        ]
        ds_concat.coords[DS_CONCAT_DUMMY_INDEX_NAME] = index
    for key in {
        true_index
        for a in annotations_re_indexed
        if (true_index := a.true_index_dim) != a.index_dim
    }:
        del ds_concat[key]
    if DUMMY_RANGE_INDEX_KEY in ds_concat:
        del ds_concat[DUMMY_RANGE_INDEX_KEY]
    ds_concat_2d = Dataset2D(ds_concat)
    return ds_concat_2d


def concat(  # noqa: PLR0912, PLR0913, PLR0915
    adatas: Collection[AnnData] | Mapping[str, AnnData],
    *,
    axis: Literal["obs", 0, "var", 1] = "obs",
    join: Join_T = "inner",
    merge: StrategiesLiteral | Callable | None = None,
    uns_merge: StrategiesLiteral | Callable | None = None,
    label: str | None = None,
    keys: Collection | None = None,
    index_unique: str | None = None,
    fill_value: Any | None = None,
    pairwise: bool = False,
    force_lazy: bool = False,
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

        For :class:`xarray.Dataset` objects, we use their :func:`xarray.merge` with `override` to stay lazy.
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
        is the delimiter between "{orig_idx}{index_unique}{key}". When `None`,
        the original indices are kept.
    fill_value
        When `join="outer"`, this is the value that will be used to fill the introduced
        indices. By default, sparse arrays are padded with zeros, while dense arrays and
        DataFrames are padded with missing values.
    pairwise
        Whether pairwise elements along the concatenated dimension should be included.
        This is False by default, since the resulting arrays are often not meaningful.
    force_lazy
        Whether to lazily concatenate elements using dask even when eager concatenation is possible.
        At the moment, this only affects obs/var and elements of obsm/varm that are xarray Datasets.

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
    ...     varm={
    ...         "ones": np.ones((2, 5)),
    ...         "rand": np.random.randn(2, 3),
    ...         "zeros": np.zeros((2, 5)),
    ...     },
    ...     uns={"a": 1, "b": 2, "c": {"c.a": 3, "c.b": 4}},
    ... )
    >>> b = ad.AnnData(
    ...     X=sparse.csr_matrix(np.array([[4, 5, 6], [7, 8, 9]])),
    ...     obs=pd.DataFrame(
    ...         {"group": ["b", "c"], "measure": [1.2, 4.3]}, index=["s3", "s4"]
    ...     ),
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
    >>> ad.concat([a, c], axis="var").to_df()
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
    >>> outer = ad.concat([a, b], join="outer")  # Joining on union of variables
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

    Using the axis’ index instead of its name

    >>> ad.concat([a, b], axis=0).to_df()  # Equivalent to axis="obs"
        var1  var2
    s1     0     1
    s2     2     3
    s3     4     5
    s4     7     8
    >>> ad.concat([a, c], axis=1).to_df()  # Equivalent to axis="var"
        var1  var2  var3  var4
    s1     0     1    10    11
    s2     2     3    12    13

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

    from anndata._core.xarray import Dataset2D
    from anndata.compat import xarray as xr

    # Argument normalization
    merge = resolve_merge_strategy(merge)
    uns_merge = resolve_merge_strategy(uns_merge)

    if isinstance(adatas, Mapping):
        if keys is not None:
            msg = (
                "Cannot specify categories in both mapping keys and using `keys`. "
                "Only specify this once."
            )
            raise TypeError(msg)
        keys, adatas = list(adatas.keys()), list(adatas.values())
    else:
        adatas = list(adatas)

    if keys is None:
        keys = np.arange(len(adatas)).astype(str)

    axis, axis_name = _resolve_axis(axis)
    alt_axis, alt_axis_name = _resolve_axis(axis=1 - axis)

    # Label column
    label_col = pd.Categorical.from_codes(
        np.repeat(np.arange(len(adatas)), [a.shape[axis] for a in adatas]),
        categories=keys,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(axis_indices(a, axis=axis)) for a in adatas], ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(
            _map_cat_to_str(label_col), sep=index_unique
        )
    concat_indices = pd.Index(concat_indices)

    alt_indices = merge_indices(
        [axis_indices(a, axis=alt_axis) for a in adatas], join=join
    )
    reindexers = [
        gen_reindexer(alt_indices, axis_indices(a, axis=alt_axis)) for a in adatas
    ]

    # Annotation for concatenation axis
    check_combinable_cols([getattr(a, axis_name).columns for a in adatas], join=join)
    annotations = [getattr(a, axis_name) for a in adatas]
    are_any_annotations_dataframes = any(
        isinstance(a, pd.DataFrame) for a in annotations
    )
    if are_any_annotations_dataframes:
        annotations_in_memory = (
            to_memory(a) if isinstance(a, Dataset2D) else a for a in annotations
        )
        concat_annot = pd.concat(
            unify_dtypes(annotations_in_memory),
            join=join,
            ignore_index=True,
        )
        concat_annot.index = concat_indices
    else:
        concat_annot = concat_dataset2d_on_annot_axis(
            annotations,
            join,
            force_lazy=force_lazy,
            concat_indices=concat_indices,
        )
    if label is not None:
        concat_annot[label] = label_col

    # Annotation for other axis
    alt_annotations = [getattr(a, alt_axis_name) for a in adatas]
    are_any_alt_annotations_dataframes = any(
        isinstance(a, pd.DataFrame) for a in alt_annotations
    )
    if are_any_alt_annotations_dataframes:
        alt_annotations_in_memory = [
            to_memory(a) if isinstance(a, Dataset2D) else a for a in alt_annotations
        ]
        alt_annot = merge_dataframes(alt_annotations_in_memory, alt_indices, merge)
    else:
        # TODO: figure out mapping of our merge to theirs instead of just taking first, although this appears to be
        # the only "lazy" setting so I'm not sure we really want that.
        # Because of xarray's merge upcasting, it's safest to simply assume that all dtypes are objects.
        annotations_with_only_dask = list(
            make_xarray_extension_dtypes_dask(
                alt_annotations, use_only_object_dtype=True
            )
        )
        annotations_with_only_dask = [
            a.ds.rename({a.true_index_dim: "merge_index"})
            for a in annotations_with_only_dask
        ]
        alt_annot = Dataset2D(
            xr.merge(annotations_with_only_dask, join=join, compat="override")
        )
        alt_annot.true_index_dim = "merge_index"

    X = concat_Xs(adatas, reindexers, axis=axis, fill_value=fill_value)

    if join == "inner":
        concat_aligned_mapping = inner_concat_aligned_mapping
        join_keys = intersect_keys
    elif join == "outer":
        concat_aligned_mapping = partial(
            outer_concat_aligned_mapping, fill_value=fill_value
        )
        join_keys = union_keys
    else:
        msg = f"{join=} should have been validated above by pd.concat"
        raise AssertionError(msg)

    layers = concat_aligned_mapping(
        [a.layers for a in adatas], axis=axis, reindexers=reindexers
    )
    concat_mapping = concat_aligned_mapping(
        [getattr(a, f"{axis_name}m") for a in adatas],
        axis=axis,
        concat_axis=0,
        index=concat_indices,
        force_lazy=force_lazy,
    )
    if pairwise:
        concat_pairwise = concat_pairwise_mapping(
            mappings=[getattr(a, f"{axis_name}p") for a in adatas],
            shapes=[a.shape[axis] for a in adatas],
            join_keys=join_keys,
        )
    else:
        concat_pairwise = {}

    # TODO: Reindex lazily, so we don't have to make those copies until we're sure we need the element
    alt_mapping = merge(
        [
            {k: r(v, axis=0) for k, v in getattr(a, f"{alt_axis_name}m").items()}
            for r, a in zip(reindexers, adatas, strict=True)
        ],
    )
    alt_pairwise = merge(
        [
            {
                k: r(r(v, axis=0), axis=1)
                for k, v in getattr(a, f"{alt_axis_name}p").items()
            }
            for r, a in zip(reindexers, adatas, strict=True)
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
        msg = (
            "Only some AnnData objects have `.raw` attribute, "
            "not concatenating `.raw` attributes."
        )
        warn(msg, UserWarning, stacklevel=2)
    return AnnData(
        **{
            "X": X,
            "layers": layers,
            axis_name: concat_annot,
            alt_axis_name: alt_annot,
            f"{axis_name}m": concat_mapping,
            f"{alt_axis_name}m": alt_mapping,
            f"{axis_name}p": concat_pairwise,
            f"{alt_axis_name}p": alt_pairwise,
            "uns": uns,
            "raw": raw,
        }
    )
