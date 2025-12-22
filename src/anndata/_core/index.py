from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import singledispatch
from itertools import repeat
from typing import TYPE_CHECKING, cast, overload

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from ..compat import AwkArray, CSArray, CSMatrix, DaskArray, XDataArray
from .xarray import Dataset2D

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from ..compat import Index, Index1D, Index1DNorm


def _normalize_indices(
    index: Index | None, names0: pd.Index, names1: pd.Index
) -> tuple[Index1DNorm | int | np.integer, Index1DNorm | int | np.integer]:
    # deal with tuples of length 1
    if isinstance(index, tuple) and len(index) == 1:
        index = index[0]
    ax0, ax1 = unpack_index(index)
    ax0 = _normalize_index(ax0, names0)
    ax1 = _normalize_index(ax1, names1)
    return ax0, ax1


def _normalize_index(  # noqa: PLR0911, PLR0912
    indexer: Index1D, index: pd.Index
) -> Index1DNorm | int | np.integer:
    # TODO: why is this here? All tests pass without it and it seems at the minimum not strict enough.
    if not isinstance(index, pd.RangeIndex) and index.dtype in (np.float64, np.int64):
        msg = f"Don’t call _normalize_index with non-categorical/string names and non-range index {index}"
        raise TypeError(msg)

    if isinstance(indexer, pd.Index | pd.Series):
        indexer = indexer.array

    # the following is insanely slow for sequences,
    # we replaced it using pandas below
    def name_idx(i):
        if isinstance(i, str):
            i = index.get_loc(i)
        return i

    if isinstance(indexer, slice):
        start = name_idx(indexer.start)
        stop = name_idx(indexer.stop)
        # string slices can only be inclusive, so +1 in that case
        if isinstance(indexer.stop, str):
            stop = None if stop is None else stop + 1
        step = indexer.step
        return slice(start, stop, step)
    elif isinstance(indexer, np.integer | int):
        return indexer
    elif isinstance(indexer, str):
        return index.get_loc(indexer)  # int
    elif isinstance(
        indexer,
        Sequence
        | np.ndarray
        | pd.api.extensions.ExtensionArray
        | CSMatrix
        | np.matrix
        | CSArray,
    ):
        if (shape := getattr(indexer, "shape", None)) is not None and (
            shape == (index.shape[0], 1) or shape == (1, index.shape[0])
        ):
            if isinstance(indexer, CSMatrix | CSArray):
                indexer = indexer.toarray()
            indexer = np.ravel(indexer)
        if not isinstance(indexer, np.ndarray):
            indexer = np.array(indexer)
            if len(indexer) == 0:
                indexer = indexer.astype(int)
        if isinstance(indexer, np.ndarray) and np.issubdtype(
            indexer.dtype, np.floating
        ):
            indexer_int = indexer.astype(int)
            if np.all((indexer - indexer_int) != 0):
                msg = f"Indexer {indexer!r} has floating point values."
                raise IndexError(msg)
        if issubclass(indexer.dtype.type, np.integer | np.floating):
            return indexer  # Might not work for range indexes
        elif issubclass(indexer.dtype.type, np.bool_):
            if indexer.shape != index.shape:
                msg = (
                    f"Boolean index does not match AnnData’s shape along this "
                    f"dimension. Boolean index has shape {indexer.shape} while "
                    f"AnnData index has shape {index.shape}."
                )
                raise IndexError(msg)
            return indexer
        else:  # indexer should be string array
            positions = index.get_indexer(indexer)
            if np.any(positions < 0):
                not_found = indexer[positions < 0]
                msg = (
                    f"Values {list(not_found)}, from {list(indexer)}, "
                    "are not valid obs/ var names or indices."
                )
                raise KeyError(msg)
            return positions  # np.ndarray[int]
    elif isinstance(indexer, XDataArray):
        if isinstance(indexer.data, DaskArray):
            return indexer.data.compute()
        return indexer.data
    msg = f"Unknown indexer {indexer!r} of type {type(indexer)}"
    raise IndexError(msg)


def _fix_slice_bounds(s: slice, length: int) -> slice:
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

    return slice(start, stop, step)


def unpack_index(index: Index) -> tuple[Index1D, Index1D]:
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
        index = tuple(i for i in index if i is not Ellipsis)
        return index
    # If index has Ellipsis, replace it with slice
    if len(index) == 2:
        index = tuple(slice(None) if i is Ellipsis else i for i in index)
        return index
    if len(index) == 1:
        index = index[0]
        if index is Ellipsis:
            index = slice(None)
        return index, slice(None)
    msg = "invalid number of indices"
    raise IndexError(msg)


@singledispatch
def _subset(
    a: np.ndarray | pd.DataFrame,
    subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm],
):
    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(DaskArray)
def _subset_dask(
    a: DaskArray, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    if len(subset_idx) > 1 and all(isinstance(x, Iterable) for x in subset_idx):
        if issparse(a._meta) and a._meta.format == "csc":
            return a[:, subset_idx[1]][subset_idx[0], :]
        return a[subset_idx[0], :][:, subset_idx[1]]
    return a[subset_idx]


@_subset.register(CSMatrix)
@_subset.register(CSArray)
def _subset_sparse(
    a: CSMatrix | CSArray,
    subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm],
):
    # Correcting for indexing behaviour of sparse.spmatrix
    if len(subset_idx) > 1 and all(isinstance(x, Iterable) for x in subset_idx):
        first_idx = subset_idx[0]
        if issubclass(first_idx.dtype.type, np.bool_):
            first_idx = np.flatnonzero(first_idx)
        subset_idx = (first_idx.reshape(-1, 1), *subset_idx[1:])
    return a[subset_idx]


@_subset.register(pd.DataFrame)
@_subset.register(Dataset2D)
def _subset_df(
    df: pd.DataFrame | Dataset2D,
    subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm],
):
    return df.iloc[subset_idx]


@_subset.register(AwkArray)
def _subset_awkarray(
    a: AwkArray, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


# Registration for SparseDataset occurs in sparse_dataset.py
@_subset.register(h5py.Dataset)
def _subset_dataset(
    d: h5py.Dataset, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    order: tuple[NDArray[np.integer] | slice, ...]
    inv_order: tuple[NDArray[np.integer] | slice, ...]
    order, inv_order = zip(*map(_index_order_and_inverse, subset_idx), strict=True)
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
def _index_order_and_inverse(axis_idx: slice) -> tuple[slice, slice]: ...
def _index_order_and_inverse(
    axis_idx: Index1DNorm,
) -> tuple[Index1DNorm, NDArray[np.integer] | slice]:
    """Order and get inverse index array."""
    if not isinstance(axis_idx, np.ndarray):
        return axis_idx, slice(None)
    if axis_idx.dtype == bool:
        axis_idx = np.flatnonzero(axis_idx)
    order = np.argsort(axis_idx)
    return axis_idx[order], np.argsort(order)


@overload
def _process_index_for_h5py(
    idx: NDArray[np.integer] | NDArray[np.bool_],
) -> tuple[NDArray[np.integer], NDArray[np.integer]]: ...
@overload
def _process_index_for_h5py(idx: slice) -> tuple[slice, None]: ...
def _process_index_for_h5py(
    idx: Index1DNorm,
) -> tuple[Index1DNorm, NDArray[np.integer] | None]:
    """Process a single index for h5py compatibility, handling sorting and duplicates."""
    if not isinstance(idx, np.ndarray):
        # Not an array (slice, integer, list) - no special processing needed
        return idx, None

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
    dataset: h5py.Dataset,
    subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm],
) -> h5py.Dataset:
    # Handle multi-dimensional indexing of h5py dataset
    # This avoids h5py's limitation with multi-dimensional fancy indexing
    # without loading the entire dataset into memory

    # Convert boolean arrays to integer arrays and handle sorting for h5py
    processed_indices: tuple[NDArray[np.integer] | slice, ...]
    reverse_indices: tuple[NDArray[np.integer] | None, ...]
    processed_indices, reverse_indices = zip(
        *map(_process_index_for_h5py, subset_idx), strict=True
    )

    # First find the index that reduces the size of the dataset the most
    i_min = np.argmin([
        _get_index_size(inds, dataset.shape[i]) / dataset.shape[i]
        for i, inds in enumerate(processed_indices)
    ])

    # Apply the most selective index first to h5py dataset
    first_index = [slice(None)] * len(processed_indices)
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


def _get_index_size(idx: Index1DNorm, dim_size: int) -> int:
    """Get size for any index type."""
    if isinstance(idx, slice):
        return len(range(*idx.indices(dim_size)))
    elif isinstance(idx, int):
        return 1
    else:  # For other types, try to get length
        return len(idx)


def make_slice(idx, dimidx: int, n: int = 2) -> tuple[slice, ...]:
    mut = list(repeat(slice(None), n))
    mut[dimidx] = idx
    return tuple(mut)


def get_vector(adata, k, coldim, idxdim, layer=None):
    # adata could be self if Raw and AnnData shared a parent
    dims = ("obs", "var")
    col = getattr(adata, coldim).columns
    idx = getattr(adata, f"{idxdim}_names")

    in_col = k in col
    in_idx = k in idx

    if (in_col + in_idx) == 2:
        msg = f"Key {k} could be found in both .{idxdim}_names and .{coldim}.columns"
        raise ValueError(msg)
    elif (in_col + in_idx) == 0:
        msg = f"Could not find key {k} in .{idxdim}_names or .{coldim}.columns."
        raise KeyError(msg)
    elif in_col:
        return getattr(adata, coldim)[k].values
    elif in_idx:
        selected_dim = dims.index(idxdim)
        idx = adata._normalize_indices(make_slice(k, selected_dim))
        a = adata._get_X(layer=layer)[idx]
    if issparse(a):
        a = a.toarray()
    return np.ravel(a)
