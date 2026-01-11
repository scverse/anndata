from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import singledispatch
from itertools import repeat
from typing import TYPE_CHECKING, cast, overload

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from ..compat import AwkArray, CSArray, CSMatrix, DaskArray, XDataArray, has_xp
from .xarray import Dataset2D

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from anndata.types import SupportsArrayApi

    from ..compat import Index, Index1D, Index1DNorm


def _normalize_indices(
    index: Index | None | tuple[IndexManager, IndexManager],
    names0: pd.Index,
    names1: pd.Index,
) -> tuple[Index1DNorm | int | np.integer, Index1DNorm | int | np.integer]:
    # deal with tuples of length 1
    if isinstance(index, tuple) and len(index) == 1:
        index = index[0]
    ax0, ax1 = unpack_index(index)
    ax0 = _normalize_index(ax0, names0)
    ax1 = _normalize_index(ax1, names1)
    return ax0, ax1


def _normalize_index(  # noqa: PLR0911, PLR0912
    indexer: Index1D | IndexManager, index: pd.Index
) -> Index1DNorm | int | np.integer:
    # TODO: why is this here? All tests pass without it and it seems at the minimum not strict enough.
    if not isinstance(index, pd.RangeIndex) and index.dtype in (np.float64, np.int64):
        msg = f"Don’t call _normalize_index with non-categorical/string names and non-range index {index}"
        raise TypeError(msg)

    if isinstance(indexer, pd.Index | pd.Series):
        indexer = indexer.array
    if isinstance(indexer, IndexManager):
        indexer = indexer.get_default()

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
    elif (use_xp := has_xp(indexer)) or isinstance(
        indexer,
        Sequence | pd.api.extensions.ExtensionArray | CSMatrix | np.matrix | CSArray,
    ):
        # convert to the 1D if it's accidentally 2D column/row vector
        # convert sparse into dense arrays if needed
        xp = indexer.__array_namespace__() if use_xp else np
        if hasattr(indexer, "shape") and (
            (indexer.shape == (index.shape[0], 1))
            or (indexer.shape == (1, index.shape[0]))
        ):
            if isinstance(indexer, CSMatrix | CSArray):
                indexer = indexer.toarray()
            indexer = xp.ravel(indexer)
        # if it is something else, convert it to numpy
        if (
            not (is_pandas := isinstance(indexer, pd.api.extensions.ExtensionArray))
            and not use_xp
        ):
            indexer = np.array(indexer)
            use_xp = True
            if len(indexer) == 0:
                indexer = indexer.astype(int)
        # https://github.com/numpy/numpy/issues/27545
        is_numpy_string = indexer.dtype == np.dtypes.StringDType()
        # if it is a float array or something along those lines, convert it to integers
        if (
            use_xp
            and not is_numpy_string
            and xp.isdtype(indexer.dtype, "real floating")
        ):
            indexer_int = xp.astype(indexer, xp.int64)
            if xp.all((indexer - indexer_int) != 0):
                msg = f"Indexer {indexer!r} has floating point values."
                raise IndexError(msg)
        if (
            is_pandas
            and (
                issubclass(indexer.dtype.type, np.integer | np.floating)
                or indexer.dtype.kind in "iuf"
            )
        ) or (
            not is_pandas
            and not is_numpy_string
            and xp.isdtype(
                indexer.dtype, ("signed integer", "unsigned integer", "real floating")
            )
        ):
            return (
                np.asarray(indexer) if is_pandas else indexer
            )  # Might not work for range indexes
        elif (
            is_pandas
            and (issubclass(indexer.dtype.type, np.bool_) or indexer.dtype.kind == "b")
        ) or (
            not is_pandas and not is_numpy_string and xp.isdtype(indexer.dtype, "bool")
        ):
            if indexer.shape != index.shape:
                msg = (
                    f"Boolean index does not match AnnData’s shape along this "
                    f"dimension. Boolean index has shape {indexer.shape} while "
                    f"AnnData index has shape {index.shape}."
                )
                raise IndexError(msg)
            return np.asarray(indexer) if is_pandas else indexer
        else:
            positions = index.get_indexer(indexer)
            if xp.any(positions < 0):
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


class IndexManager:
    """Manages index arrays across multiple devices for efficient subsetting.

    This class stores index arrays (for subsetting AnnData views) on different
    devices and can produce the appropriate index array when subsetting arrays
    that live on those devices. For numpy/pandas/sparse operations, use
    ``np.asarray(index_manager)`` or let the ``__array__`` method be invoked.
    """

    _manager: dict[str, SupportsArrayApi]

    def __init__(self, *, device: str, arr: SupportsArrayApi):
        self._manager = {device: arr}

    @classmethod
    def from_array(cls, arr: SupportsArrayApi):
        """Create an IndexManager from an array-api compatible array."""
        return cls(device=arr.device, arr=arr)

    def __array__(
        self,
        dtype: np.dtype | None = None,
        copy: bool | None = None,  # noqa: FBT001
    ) -> np.ndarray:
        """Return numpy array for compatibility with numpy/pandas/sparse operations."""
        if "cpu" not in self:
            arr = self._manager[next(iter(self.keys()))]
            self._manager["cpu"] = np.from_dlpack(arr.to_device("cpu"))
        res = np.from_dlpack(self._manager["cpu"])
        return res.copy() if copy else res

    def __contains__(self, device: str) -> bool:
        """Check if an index array exists for the given device."""
        return device in self._manager

    def __len__(self) -> int:
        """Return the length of the index array."""
        return self.get_default().shape[0]

    def __iter__(self):
        """Iterate over the index values (as numpy for compatibility)."""
        return iter(np.asarray(self))

    def keys(self):
        """Return the devices for which index arrays are available."""
        return self._manager.keys()

    def get_default(self):
        """Returns the first key added i.e., a no-copy index, useful for getting an array-api compatible index on some device."""
        return self._manager[next(iter(self._manager))]

    @property
    def dtype(self):
        """Return the dtype of the index array."""
        arr = self._manager[next(iter(self._manager))]
        return arr.dtype

    def add_array(self, arr: SupportsArrayApi):
        """Add an index array for a specific device."""
        self._manager[arr.device] = arr

    def get_array_with_device(self, device: str) -> SupportsArrayApi:
        """Get the index array for a specific device."""
        return self._manager[device]

    def get_for_array(self, arr: SupportsArrayApi) -> SupportsArrayApi:
        """Get an index array on the same device as the input array.

        If an index doesn't exist for the array's device, it will be
        created by transferring from an existing device. If an index
        exists but has a different array-api implementation, it will
        be converted to match the input array's array-api.
        """
        device = arr.device
        xp = arr.__array_namespace__()

        if device in self._manager:
            existing = self._manager[device]
            existing_xp = existing.__array_namespace__()
            # Check if array-api implementations match
            if existing_xp is xp:
                return existing
            # Convert to matching array-api implementation
            return xp.from_dlpack(existing, device=device)

        # Transfer from an existing device
        src_arr = self._manager[next(iter(self._manager))]
        self._manager[device] = xp.from_dlpack(src_arr.to_device(device))
        return self._manager[device]


def _ensure_numpy_idx(
    subset_idx: tuple,
) -> tuple:
    """Convert IndexManager instances to numpy arrays in a tuple of indices."""
    return tuple(
        np.asarray(idx) if isinstance(idx, IndexManager) else idx for idx in subset_idx
    )


def _prepare_array_api_idx(
    a: SupportsArrayApi,
    subset_idx: tuple,
) -> tuple:
    """Prepare indices for array-api subsetting, implementing np.ix_-like behavior.

    For array-api arrays, this returns indices on the same device as `a` and
    implements coordinate-based indexing similar to np.ix_.
    """
    xp = a.__array_namespace__()

    def get_idx(idx):
        if isinstance(idx, IndexManager):
            return idx.get_for_array(a)
        elif isinstance(idx, slice) or has_xp(idx):
            return idx
        else:
            # Convert numpy/list to array-api array on the target device
            return xp.asarray(idx, device=a.device)

    processed = tuple(get_idx(idx) for idx in subset_idx)

    # Implement np.ix_-like behavior for 2D indexing
    if len(processed) == 2 and all(not isinstance(x, slice) for x in processed):
        # For 2D fancy indexing, we need to reshape indices
        # np.ix_ produces (n,1) and (1,m) shaped arrays for coordinate indexing
        row_idx, col_idx = processed
        # Use outer indexing pattern: first index rows, then index columns
        # This is equivalent to a[row_idx][:, col_idx] but in one operation
        row_idx = xp.reshape(row_idx, (-1, 1))
        col_idx = xp.reshape(col_idx, (1, -1))
        return (row_idx, col_idx)

    return processed


@singledispatch
def _subset(
    a: np.ndarray | pd.DataFrame,
    subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm],
):
    """Select a subset of array `a` using the given indices.

    For numpy arrays with array indices (not slices), this uses np.ix_ for
    coordinate-based indexing. For array-api arrays, it uses device-aware
    indexing with IndexManager support.
    """
    # Check if this is an array-api array (not numpy)
    if has_xp(a) and not isinstance(a, np.ndarray):
        # Use array-api aware indexing
        subset_idx = _prepare_array_api_idx(a, subset_idx)
        return a[subset_idx]

    # For numpy arrays and other types, ensure we have numpy indices
    subset_idx = _ensure_numpy_idx(subset_idx)

    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(DaskArray)
def _subset_dask(
    a: DaskArray, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    # Dask uses numpy-style indexing, convert IndexManager to numpy
    subset_idx = _ensure_numpy_idx(subset_idx)
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
    # Sparse matrices use numpy indexing, convert IndexManager to numpy
    subset_idx = _ensure_numpy_idx(subset_idx)
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
    # DataFrames use numpy indexing, convert IndexManager to numpy
    subset_idx = _ensure_numpy_idx(subset_idx)
    return df.iloc[subset_idx]


@_subset.register(AwkArray)
def _subset_awkarray(
    a: AwkArray, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    # Awkward arrays use numpy indexing, convert IndexManager to numpy
    subset_idx = _ensure_numpy_idx(subset_idx)
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
