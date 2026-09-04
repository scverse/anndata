from __future__ import annotations

from collections.abc import Iterable, Sequence
from enum import Enum
from functools import singledispatch, wraps
from itertools import repeat
from types import EllipsisType
from typing import TYPE_CHECKING, Literal, NamedTuple, cast, overload

import h5py
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse

from .._settings import settings
from ..compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    DaskArray,
    IndexManager,
    XDataArray,
    has_xp_base,
)
from ..types import SupportsArrayApiBase
from ._dataframe_backend import DataFrameLike, frame_annotation_columns, subset_frame
from .xarray import Dataset2D

if TYPE_CHECKING:
    import sys
    from collections.abc import Callable
    from typing import TypeAlias

    if sys.version_info >= (3, 13):
        from typing import TypeIs
    else:
        from typing_extensions import TypeIs

    from ..acc import Array
    from ..typing import AlignedArray, Index, Index1D, _Index1DNorm
    from .anndata import AnnData
    from .raw import Raw


__all__ = [
    "_ensure_numpy_idx",
    "_fix_slice_bounds",
    "_get_vector_ambiguous",
    "_normalize_indices",
    "_safe_fancy_index_h5py",
    "_subset",
    "array_api_ix",
    "make_slice",
    "unpack_index",
]

ArrayApiDtype = Literal["real floating", "signed integer", "unsigned integer", "bool"]
DtypeKind = Literal["f", "i", "u", "b"]


def _normalize_indices(
    index: Index[IndexManager], names0: pd.Index, names1: pd.Index
) -> tuple[_Index1DNorm | int | np.integer, _Index1DNorm | int | np.integer]:
    # deal with tuples of length 1
    if isinstance(index, tuple) and len(index) == 1:
        index = index[0]
    ax0, ax1 = unpack_index(index)
    return _normalize_index(ax0, names0), _normalize_index(ax1, names1)


class _XrDtV(NamedTuple):
    type: type[np.generic]
    array_api_str: ArrayApiDtype
    kind: DtypeKind

    def __call__(
        self,
        indexer: SupportsArrayApiBase
        | pd.api.extensions.ExtensionArray
        | np.matrix
        | pd.MultiIndex,
    ) -> bool:
        return (
            has_xp_base(indexer)
            and indexer.__array_namespace__().isdtype(indexer.dtype, self.array_api_str)
        ) or (
            isinstance(indexer, pd.api.extensions.ExtensionArray)
            and (
                issubclass(indexer.dtype.type, self.type)
                or indexer.dtype.kind == self.kind
            )
        )


class XArrayDtype(_XrDtV, Enum):
    Float = (np.floating, "real floating", "f")
    SignedInt = (np.signedinteger, "signed integer", "i")
    UnsignedInt = (np.unsignedinteger, "unsigned integer", "u")
    Bool = (np.bool, "bool", "b")


@singledispatch
def _gen_anndata_index(
    indexer: Index1D, index: pd.Index
) -> _Index1DNorm | int | np.integer:
    msg = f"Unknown indexer {indexer!r} of type {type(indexer)}"
    raise IndexError(msg)


@_gen_anndata_index.register(slice)
def _from_slice(indexer: slice, index: pd.Index) -> slice:
    def name_idx(i):
        if isinstance(i, str):
            i = index.get_loc(i)
        return i

    start = name_idx(indexer.start)
    stop = name_idx(indexer.stop)
    # string slices can only be inclusive, so +1 in that case
    if isinstance(indexer.stop, str):
        stop = None if stop is None else stop + 1
    step = indexer.step
    return slice(start, stop, step)


@_gen_anndata_index.register(np.integer | int)
def _from_int(indexer: np.integer | int, index: pd.Index) -> np.integer | int:
    return indexer


@_gen_anndata_index.register(str)
def _from_str(indexer: str, index: pd.Index) -> int | slice | NDArray[np.bool_]:
    # non-unique names resolve to a slice or a mask rather than a position
    return index.get_loc(indexer)


@_gen_anndata_index.register(XDataArray)
def _from_xarray(indexer: XDataArray, index: pd.Index) -> np.ndarray:
    if isinstance(indexer.data, DaskArray):
        return indexer.data.compute()
    return indexer.data


@_gen_anndata_index.register(CSMatrix | CSArray)
def _from_sparse(
    indexer: CSMatrix | CSArray,
    index: pd.Index,
) -> _Index1DNorm | int | np.integer:
    return _gen_anndata_index(indexer.toarray(), index)


@_gen_anndata_index.register(Sequence)
def _from_sequence(
    indexer: Sequence,
    index: pd.Index,
) -> _Index1DNorm | int | np.integer:
    arr = np.array(indexer)
    if len(arr) == 0:
        arr = arr.astype(int)
    return _gen_anndata_index(arr, index)


def is_pandas_idx(
    indexer: object,
) -> TypeIs[pd.api.extensions.ExtensionArray | pd.MultiIndex]:
    return isinstance(indexer, pd.api.extensions.ExtensionArray | pd.MultiIndex)


@_gen_anndata_index.register(
    SupportsArrayApiBase | pd.api.extensions.ExtensionArray | np.matrix | pd.MultiIndex
)
def _from_array(
    indexer: SupportsArrayApiBase
    | pd.api.extensions.ExtensionArray
    | np.matrix
    | pd.MultiIndex,
    index: pd.Index,
) -> SupportsArrayApiBase:
    # convert to the 1D if it's accidentally 2D column/row vector
    # convert sparse into dense arrays if needed
    xp = indexer.__array_namespace__() if has_xp_base(indexer) else np
    if indexer.shape == (index.shape[0], 1) or indexer.shape == (1, index.shape[0]):
        # the array API has no `ravel`, and `np.matrix` is always 2D
        indexer = (
            np.asarray(indexer).reshape(-1)
            if isinstance(indexer, np.matrix)
            else xp.reshape(indexer, (-1,))
        )
    # https://github.com/numpy/numpy/issues/27545
    is_numpy_string = indexer.dtype == np.dtypes.StringDType()
    # MultiIndex objects are not turned into arrays in _normalize_index, so handle them explicitly
    if not is_numpy_string or is_pandas_idx(indexer):
        # if it is a float array or something along those lines, convert it to integers
        if XArrayDtype.Float(indexer):
            indexer_int = xp.astype(indexer, xp.int64)
            if xp.all((indexer - indexer_int) != 0):
                msg = f"Indexer {indexer!r} has floating point values."
                raise IndexError(msg)
            return indexer_int
        if XArrayDtype.SignedInt(indexer) or XArrayDtype.UnsignedInt(indexer):
            # Might not work for range indexes
            return np.asarray(indexer) if is_pandas_idx(indexer) else indexer
        if XArrayDtype.Bool(indexer):
            if indexer.shape != index.shape:
                msg = (
                    f"Boolean index does not match AnnData’s shape along this "
                    f"dimension. Boolean index has shape {indexer.shape} while "
                    f"AnnData index has shape {index.shape}."
                )
                raise IndexError(msg)
            return np.asarray(indexer) if is_pandas_idx(indexer) else indexer
        if not isinstance(indexer, np.ndarray) and has_xp_base(indexer):
            msg = f"indexer is array-api compatible but has unsupported dtype: {indexer.dtype}"
            raise ValueError(msg)
    # indexer is a string array here; `get_indexer` needs a collection
    names = (
        indexer
        if isinstance(indexer, pd.Index)
        else pd.Index(
            indexer
            if isinstance(indexer, pd.api.extensions.ExtensionArray)
            else np.asarray(indexer)
        )
    )
    positions = index.get_indexer(names)
    if xp.any(positions < 0):
        not_found = names[positions < 0]
        msg = (
            f"Values {list(not_found)}, from {list(names)}, "
            "are not valid obs/ var names or indices."
        )
        raise KeyError(msg)
    return positions  # np.ndarray[int]


def _normalize_index(
    indexer: Index1D[IndexManager], index: pd.Index
) -> _Index1DNorm | int | np.integer:

    if isinstance(indexer, pd.Index | pd.Series) and (
        not isinstance(indexer, pd.MultiIndex) or settings.restrict_index_types
    ):
        indexer = indexer.array
    if isinstance(indexer, IndexManager):
        indexer = indexer.get_default()

    return _gen_anndata_index(indexer, index)


def _fix_slice_bounds(
    s: slice[int | None, int | None, int | None], length: int
) -> slice[int, int, int]:
    """The slice will be clipped to length, and the step won't be None.

    E.g. infer None valued attributes.
    """
    step = s.step if s.step is not None else 1

    # slice constructor would have errored if step was 0
    if step > 0:
        start = s.start if s.start is not None else 0
        stop = s.stop if s.stop is not None else length
    elif step < 0:
        # Reverse
        start = s.start if s.start is not None else length
        stop = s.stop if s.stop is not None else 0
    else:
        msg = "step must be non-zero"
        raise AssertionError(msg)

    return slice(start, stop, step)


def unpack_index[M: IndexManager](index: Index[M]) -> tuple[Index1D[M], Index1D[M]]:
    if not isinstance(index, tuple):
        if index is Ellipsis:
            index = slice(None)
        return index, slice(None)
    num_ellipsis = sum(i is Ellipsis for i in index)
    if num_ellipsis > 1:
        msg = "an index can only have a single ellipsis ('...')"
        raise IndexError(msg)
    # If index has Ellipsis, filter it out (and if not, error)
    if len(index) > 2:
        if not num_ellipsis:
            msg = "Received a length 3 index without an ellipsis"
            raise IndexError(msg)
        ax0, ax1 = (i for i in index if not isinstance(i, EllipsisType))
        return ax0, ax1
    # If index has Ellipsis, replace it with slice
    if len(index) == 2:
        ax0, ax1 = (slice(None) if isinstance(i, EllipsisType) else i for i in index)
        return ax0, ax1
    if len(index) == 1:
        index = index[0]
        if index is Ellipsis:
            index = slice(None)
        return index, slice(None)
    msg = "invalid number of indices"
    raise IndexError(msg)


Idx1D: TypeAlias = (  # noqa: UP040
    slice
    | NDArray[np.bool_]
    | NDArray[np.integer]
    | SupportsArrayApiBase
    | IndexManager
)
"""`_Index1DNorm[IndexManager]`, spelled out so it resolves at runtime."""

SubsetIdx: TypeAlias = tuple[Idx1D] | tuple[Idx1D, Idx1D]  # noqa: UP040
"""A one- or two-dimensional index as stored on a view."""

NumpyIdx1D: TypeAlias = slice | NDArray[np.bool_] | NDArray[np.integer]  # noqa: UP040
NumpySubsetIdx: TypeAlias = tuple[NumpyIdx1D] | tuple[NumpyIdx1D, NumpyIdx1D]  # noqa: UP040
"""A :data:`SubsetIdx` whose array parts are numpy arrays."""


def _index_manager_to_numpy_idx(idx: _Index1DNorm[IndexManager]) -> _Index1DNorm:
    """Unwrap an `IndexManager`; every other index is passed through untouched."""
    return np.asarray(idx) if isinstance(idx, IndexManager) else idx


def _as_numpy_idx(idx: _Index1DNorm[IndexManager]) -> NumpyIdx1D:
    """Materialise an index for consumers that only take numpy, e.g. pandas and h5py."""
    return idx if isinstance(idx, slice | np.ndarray) else np.asarray(idx)


def _index_manager_to_numpy_idx_in_tuple(subset_idx: SubsetIdx) -> SubsetIdx:
    """Unwrap any `IndexManager` in a tuple of indices."""
    if len(subset_idx) == 1:
        return (_index_manager_to_numpy_idx(subset_idx[0]),)
    return (
        _index_manager_to_numpy_idx(subset_idx[0]),
        _index_manager_to_numpy_idx(subset_idx[1]),
    )


def _as_numpy_subset_idx(subset_idx: SubsetIdx) -> NumpySubsetIdx:
    """Materialise a tuple of indices for the numpy-only implementations."""
    if len(subset_idx) == 1:
        return (_as_numpy_idx(subset_idx[0]),)
    return (_as_numpy_idx(subset_idx[0]), _as_numpy_idx(subset_idx[1]))


def _ensure_numpy_idx[T, R](
    func: Callable[[T, NumpySubsetIdx], R],
) -> Callable[[T, SubsetIdx], R]:
    """Materialise non-numpy indices for the wrapped implementation."""

    @wraps(func)
    def _ensure(a: T, subset_idx: SubsetIdx) -> R:
        return func(a, _as_numpy_subset_idx(subset_idx))

    return _ensure


def array_api_ix(*args: SupportsArrayApiBase) -> tuple[SupportsArrayApiBase, ...]:
    """Construct an open mesh from multiple sequences.

    Vendored version of `numpy.ix_` for the array-api.

    For each sequence `args[i]`, it returns an array with `.ndim == len(args)`,
    `.size == .shape[i] == len(args[i])` (i.e. `shape[…] == 1` for each other dimension)
    """
    out = []
    n_dims = len(args)
    for k, new in enumerate(args):
        xp = new.__array_namespace__()
        if new.ndim != 1:  # pragma: no cover
            msg = "Cross index must be 1 dimensional"
            raise ValueError(msg)
        if xp.isdtype(new.dtype, "bool"):
            (new,) = xp.nonzero(new)
        new = xp.reshape(new, (1,) * k + (new.size,) + (1,) * (n_dims - k - 1))
        out.append(new)
    return tuple(out)


def _prepare_array_api_idx(
    a: SupportsArrayApiBase, subset_idx: SubsetIdx
) -> tuple[slice | SupportsArrayApiBase, ...]:
    """Prepare indices for array-api compatible subsetting."""
    xp = a.__array_namespace__()

    def get_idx(idx):
        if isinstance(idx, IndexManager):
            return idx.get_for_array(a)
        elif isinstance(idx, slice) or has_xp_base(idx):
            return idx
        else:  # pragma: no cover
            # Convert numpy/list to array-api array on the target device
            # In theory should be unreachable so this is a last resort since xp.asarray is pretty undefined.
            return xp.asarray(idx, device=a.device)

    maybe_array_api_idxs = tuple(get_idx(idx) for idx in subset_idx)
    if all(isinstance(x, type(a)) for x in maybe_array_api_idxs):
        return array_api_ix(*maybe_array_api_idxs)

    return maybe_array_api_idxs


@singledispatch
def _subset_dispatch(
    a: AlignedArray | AnnData, subset_idx: SubsetIdx
) -> AlignedArray | AnnData:
    """Subset a numpy or array-API array; everything else is registered below.

    For numpy arrays with array indices (not slices), this uses np.ix_ for
    coordinate-based indexing. For array-api arrays, it uses device-aware
    indexing with IndexManager support.
    """
    if not isinstance(a, np.ndarray):
        if has_xp_base(a):
            return a[_prepare_array_api_idx(a, subset_idx)]
        # zarr and cupy arrays are indexed just like numpy ones, but aren’t typed as such
        a = cast("np.ndarray", a)

    numpy_idx = _as_numpy_subset_idx(subset_idx)

    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, Iterable) for x in numpy_idx):
        return a[np.ix_(*numpy_idx)]
    return a[numpy_idx]


def _subset[T: AlignedArray](
    a: T, subset_idx: tuple[_Index1DNorm[IndexManager], ...]
) -> T | np.ndarray:
    """Select a subset of `a` using the given indices.

    Dispatch cannot express “returns what it was given”, so the generic
    signature lives here and the implementations register on
    `_subset_dispatch`.
    """
    match subset_idx:
        case (i,):
            idx: SubsetIdx = (i,)
        case (i, j):
            idx = (i, j)
        case _:
            msg = f"Can only subset along one or two axes, got {len(subset_idx)}."
            raise IndexError(msg)
    return cast("T | np.ndarray", _subset_dispatch(a, idx))


@_subset_dispatch.register(DaskArray)
@_ensure_numpy_idx
def _subset_dask(a: DaskArray, subset_idx: NumpySubsetIdx) -> DaskArray:
    if len(subset_idx) > 1 and all(isinstance(x, Iterable) for x in subset_idx):
        if isinstance(a._meta, sparse.csc_matrix | sparse.csc_array):
            return a[:, subset_idx[1]][subset_idx[0], :]
        return a[subset_idx[0], :][:, subset_idx[1]]
    return a[subset_idx]


@_subset_dispatch.register(CSMatrix)
@_subset_dispatch.register(CSArray)
@_ensure_numpy_idx
def _subset_sparse(
    a: CSMatrix | CSArray, subset_idx: NumpySubsetIdx
) -> CSMatrix | CSArray:
    if len(subset_idx) == 1:
        return _subset_sparse_axis(a, subset_idx[0], axis=0)
    rows, cols = subset_idx
    # a pair of arrays would index coordinate-wise, so subset one axis at a time
    if _is_full_slice(rows):
        return _subset_sparse_axis(a, cols, axis=1)
    if _is_full_slice(cols):
        return _subset_sparse_axis(a, rows, axis=0)
    return _subset_sparse_axis(_subset_sparse_axis(a, rows, axis=0), cols, axis=1)


def _is_full_slice(idx: NumpyIdx1D) -> bool:
    return isinstance(idx, slice) and idx == slice(None)


def _subset_sparse_axis(
    a: CSMatrix | CSArray, idx: NumpyIdx1D, *, axis: Literal[0, 1]
) -> CSMatrix | CSArray:
    # slice and array indices are served by separate `__getitem__` overloads
    if isinstance(idx, slice):
        return a[idx, :] if axis == 0 else a[:, idx]
    return a[idx, :] if axis == 0 else a[:, idx]


@_subset_dispatch.register(pd.DataFrame)
@_ensure_numpy_idx
def _subset_df(df: pd.DataFrame, subset_idx: NumpySubsetIdx) -> pd.DataFrame:
    if len(subset_idx) == 1:
        return df.iloc[subset_idx[0]]
    return df.iloc[subset_idx[0], subset_idx[1]]


@_subset_dispatch.register(DataFrameLike)
@_ensure_numpy_idx
def _subset_frame[T: DataFrameLike](df: T, subset_idx: NumpySubsetIdx) -> T:
    return cast("T", subset_frame(df, subset_idx))


@_subset_dispatch.register(Dataset2D)
@_ensure_numpy_idx
def _subset_dataset2d(ds: Dataset2D, subset_idx: NumpySubsetIdx) -> Dataset2D:
    return ds.iloc[subset_idx]


@_subset_dispatch.register(AwkArray)
@_ensure_numpy_idx
def _subset_awkarray(a: AwkArray, subset_idx: NumpySubsetIdx) -> AwkArray:
    if all(isinstance(x, Iterable) for x in subset_idx):
        return a[np.ix_(*subset_idx)]
    return a[subset_idx]


# Registration for SparseDataset occurs in sparse_dataset.py
@_subset_dispatch.register(h5py.Dataset)
@_ensure_numpy_idx
def _subset_dataset(d: h5py.Dataset, subset_idx: NumpySubsetIdx) -> np.ndarray:
    order: tuple[NDArray[np.integer] | slice, ...]
    inv_order: tuple[NDArray[np.integer] | slice, ...]
    order, inv_order = zip(
        *(_index_order_and_inverse(idx) for idx in subset_idx), strict=True
    )
    # check for duplicates or multi-dimensional fancy indexing
    array_dims = [i for i in order if isinstance(i, np.ndarray)]
    has_duplicates = any(len(np.unique(i)) != len(i) for i in array_dims)
    # Use safe indexing if there are duplicates OR multiple array dimensions
    # (h5py doesn't support multi-dimensional fancy indexing natively)
    if has_duplicates or len(array_dims) > 1:
        # For multi-dimensional indexing, bypass the sorting logic and use original indices
        return _safe_fancy_index_h5py(d, subset_idx)
    # from hdf5, then to real order
    return d[order][inv_order]


@overload
def _index_order_and_inverse(
    axis_idx: NDArray[np.integer] | NDArray[np.bool_],
) -> tuple[NDArray[np.integer], NDArray[np.integer]]: ...
@overload
def _index_order_and_inverse(
    axis_idx: _Index1DNorm,
) -> tuple[NDArray[np.integer] | slice, NDArray[np.integer] | slice]: ...
def _index_order_and_inverse(
    axis_idx: _Index1DNorm,
) -> tuple[NDArray[np.integer] | slice, NDArray[np.integer] | slice]:
    """Order and get inverse index array."""
    if isinstance(axis_idx, slice):
        return axis_idx, slice(None)
    # h5py only understands numpy indices
    idx = np.asarray(axis_idx)
    if idx.dtype == bool:
        idx = np.flatnonzero(idx)
    order = np.argsort(idx)
    return idx[order], np.argsort(order)


def _process_index_for_h5py(
    idx: _Index1DNorm,
) -> tuple[NDArray[np.integer] | slice, NDArray[np.integer] | None]:
    """Process a single index for h5py compatibility, handling sorting and duplicates."""
    if isinstance(idx, slice):
        # Not an array - no special processing needed
        return idx, None

    # h5py only understands numpy indices
    idx = np.asarray(idx)
    if idx.dtype == bool:
        idx = np.flatnonzero(idx)

    # For h5py fancy indexing, we need sorted indices
    # But we also need to track how to reverse the sorting
    unique, inverse = np.unique(idx, return_inverse=True)
    return (
        # Has duplicates - use unique + inverse mapping approach
        (unique, inverse)
        if len(unique) != len(idx)
        # No duplicates - just sort and track reverse mapping
        else _index_order_and_inverse(idx)
    )


def _safe_fancy_index_h5py(
    dataset: h5py.Dataset, subset_idx: NumpySubsetIdx
) -> np.ndarray:
    # Handle multi-dimensional indexing of h5py dataset
    # This avoids h5py's limitation with multi-dimensional fancy indexing
    # without loading the entire dataset into memory

    # Convert boolean arrays to integer arrays and handle sorting for h5py
    processed_indices: tuple[NDArray[np.integer] | slice, ...]
    reverse_indices: tuple[NDArray[np.integer] | None, ...]
    processed_indices, reverse_indices = zip(
        *(_process_index_for_h5py(idx) for idx in subset_idx), strict=True
    )

    # First find the index that reduces the size of the dataset the most
    i_min = np.argmin([
        _get_index_size(inds, dataset.shape[i]) / dataset.shape[i]
        for i, inds in enumerate(processed_indices)
    ])

    # Apply the most selective index first to h5py dataset
    first_index: list[NDArray[np.integer] | slice] = [slice(None)] * len(
        processed_indices
    )
    first_index[i_min] = processed_indices[i_min]
    in_memory_array = cast("np.ndarray", dataset[tuple(first_index)])

    # Apply remaining indices to the numpy array
    remaining_indices = list(processed_indices)
    remaining_indices[i_min] = slice(None)  # Already applied
    result = in_memory_array[tuple(remaining_indices)]

    # Now apply reverse mappings to get the original order
    for dim, reverse_map in enumerate(reverse_indices):
        if reverse_map is not None:
            result = result.take(reverse_map, axis=dim)

    return result


def _get_index_size(idx: _Index1DNorm, dim_size: int) -> int:
    """Get size for any index type."""
    if isinstance(idx, slice):
        return len(range(*idx.indices(dim_size)))
    elif isinstance(idx, int):
        return 1
    else:  # For other types, try to get length
        return idx.shape[0]


def make_slice(idx, dimidx: int, n: int = 2) -> tuple[slice, ...]:
    mut = list(repeat(slice(None), n))
    mut[dimidx] = idx
    return tuple(mut)


def _get_vector_ambiguous(
    adata: AnnData | Raw, k: str, dim: Literal["obs", "var"], layer: str | None = None
) -> Array:
    from ..acc import A

    idxdim = "var" if dim == "obs" else "obs"
    in_annotation = k in frame_annotation_columns(getattr(adata, dim), dim=dim)
    match (in_annotation, k in getattr(adata, f"{idxdim}_names")):
        case True, True:
            msg = f"Key {k} could be found in both .{idxdim}_names and .{dim}.columns"
            raise ValueError(msg)
        case False, False:
            msg = f"Could not find key {k} in .{idxdim}_names or .{dim}.columns."
            raise KeyError(msg)
        case True, False:
            ref = (A.obs if dim == "obs" else A.var)[k]
        case False, True:
            acc = A.layers[layer]
            ref = acc[k, :] if idxdim == "obs" else acc[:, k]
    return adata[ref]
