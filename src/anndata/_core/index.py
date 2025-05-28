from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import singledispatch
from itertools import repeat
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from ..compat import AwkArray, CSArray, CSMatrix, DaskArray, XDataArray
from .xarray import Dataset2D

if TYPE_CHECKING:
    from ..compat import Index, Index1D


def _normalize_indices(
    index: Index | None, names0: pd.Index, names1: pd.Index
) -> tuple[slice, slice]:
    # deal with tuples of length 1
    if isinstance(index, tuple) and len(index) == 1:
        index = index[0]
    # deal with pd.Series
    if isinstance(index, pd.Series):
        index: Index = index.values
    if isinstance(index, tuple):
        # TODO: The series should probably be aligned first
        index = tuple(i.values if isinstance(i, pd.Series) else i for i in index)
    ax0, ax1 = unpack_index(index)
    ax0 = _normalize_index(ax0, names0)
    ax1 = _normalize_index(ax1, names1)
    return ax0, ax1


def _normalize_index(  # noqa: PLR0911, PLR0912
    indexer: slice
    | np.integer
    | int
    | str
    | Sequence[bool | int | np.integer]
    | np.ndarray
    | pd.Index,
    index: pd.Index,
) -> slice | int | np.ndarray:  # ndarray of int or bool
    # TODO: why is this here? All tests pass without it and it seems at the minimum not strict enough.
    if not isinstance(index, pd.RangeIndex) and index.dtype in (np.float64, np.int64):
        msg = f"Don’t call _normalize_index with non-categorical/string names and non-range index {index}"
        raise TypeError(msg)

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
        indexer, Sequence | np.ndarray | pd.Index | CSMatrix | np.matrix | CSArray
    ):
        if hasattr(indexer, "shape") and (
            (indexer.shape == (index.shape[0], 1))
            or (indexer.shape == (1, index.shape[0]))
        ):
            if isinstance(indexer, CSMatrix | CSArray):
                indexer = indexer.toarray()
            indexer = np.ravel(indexer)
        if not isinstance(indexer, np.ndarray | pd.Index):
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
    raise IndexError()


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
def _subset(a: np.ndarray | pd.DataFrame, subset_idx: Index):
    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(DaskArray)
def _subset_dask(a: DaskArray, subset_idx: Index):
    if len(subset_idx) > 1 and all(isinstance(x, Iterable) for x in subset_idx):
        if issparse(a._meta) and a._meta.format == "csc":
            return a[:, subset_idx[1]][subset_idx[0], :]
        return a[subset_idx[0], :][:, subset_idx[1]]
    return a[subset_idx]


@_subset.register(CSMatrix)
@_subset.register(CSArray)
def _subset_sparse(a: CSMatrix | CSArray, subset_idx: Index):
    # Correcting for indexing behaviour of sparse.spmatrix
    if len(subset_idx) > 1 and all(isinstance(x, Iterable) for x in subset_idx):
        first_idx = subset_idx[0]
        if issubclass(first_idx.dtype.type, np.bool_):
            first_idx = np.where(first_idx)[0]
        subset_idx = (first_idx.reshape(-1, 1), *subset_idx[1:])
    return a[subset_idx]


@_subset.register(pd.DataFrame)
def _subset_df(df: pd.DataFrame, subset_idx: Index):
    return df.iloc[subset_idx]


@_subset.register(AwkArray)
def _subset_awkarray(a: AwkArray, subset_idx: Index):
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(Dataset2D)
def _(a: Dataset2D, subset_idx: Index):
    key = a.index_dim
    # xarray seems to have some code looking for a second entry in tuples
    if isinstance(subset_idx, tuple) and len(subset_idx) == 1:
        subset_idx = subset_idx[0]
    return a.isel(**{key: subset_idx})


# Registration for SparseDataset occurs in sparse_dataset.py
@_subset.register(h5py.Dataset)
def _subset_dataset(d, subset_idx):
    if not isinstance(subset_idx, tuple):
        subset_idx = (subset_idx,)
    ordered = list(subset_idx)
    rev_order = [slice(None) for _ in range(len(subset_idx))]
    for axis, axis_idx in enumerate(ordered.copy()):
        if isinstance(axis_idx, np.ndarray):
            if axis_idx.dtype == bool:
                axis_idx = np.where(axis_idx)[0]
            order = np.argsort(axis_idx)
            ordered[axis] = axis_idx[order]
            rev_order[axis] = np.argsort(order)
    # from hdf5, then to real order
    return d[tuple(ordered)][tuple(rev_order)]


def make_slice(idx, dimidx, n=2):
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
