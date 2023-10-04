from __future__ import annotations

import collections.abc as cabc
from collections.abc import Sequence
from functools import singledispatch
from itertools import repeat

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, issparse, spmatrix

from ..compat import AwkArray, DaskArray, Index, Index1D


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
        if len(index) > 2:
            raise ValueError("AnnData can only be sliced in rows and columns.")
        # deal with pd.Series
        # TODO: The series should probably be aligned first
        if isinstance(index[1], pd.Series):
            index = index[0], index[1].values
        if isinstance(index[0], pd.Series):
            index = index[0].values, index[1]
    ax0, ax1 = unpack_index(index)
    ax0 = _normalize_index(ax0, names0)
    ax1 = _normalize_index(ax1, names1)
    return ax0, ax1


def _normalize_index(
    indexer: slice
    | np.integer
    | int
    | str
    | Sequence[int | np.integer]
    | np.ndarray
    | pd.Index,
    index: pd.Index,
) -> slice | int | np.ndarray:  # ndarray of int
    if not isinstance(index, pd.RangeIndex):
        assert (
            index.dtype != float and index.dtype != int
        ), "Don’t call _normalize_index with non-categorical/string names"

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
    elif isinstance(indexer, (np.integer, int)):
        return indexer
    elif isinstance(indexer, str):
        return index.get_loc(indexer)  # int
    elif isinstance(indexer, (Sequence, np.ndarray, pd.Index, spmatrix, np.matrix)):
        if hasattr(indexer, "shape") and (
            (indexer.shape == (index.shape[0], 1))
            or (indexer.shape == (1, index.shape[0]))
        ):
            if isinstance(indexer, spmatrix):
                indexer = indexer.toarray()
            indexer = np.ravel(indexer)
        if not isinstance(indexer, (np.ndarray, pd.Index)):
            indexer = np.array(indexer)
        if issubclass(indexer.dtype.type, (np.integer, np.floating)):
            return indexer  # Might not work for range indexes
        elif issubclass(indexer.dtype.type, np.bool_):
            if indexer.shape != index.shape:
                raise IndexError(
                    f"Boolean index does not match AnnData’s shape along this "
                    f"dimension. Boolean index has shape {indexer.shape} while "
                    f"AnnData index has shape {index.shape}."
                )
            positions = np.where(indexer)[0]
            return positions  # np.ndarray[int]
        else:  # indexer should be string array
            positions = index.get_indexer(indexer)
            if np.any(positions < 0):
                not_found = indexer[positions < 0]
                raise KeyError(
                    f"Values {list(not_found)}, from {list(indexer)}, "
                    "are not valid obs/ var names or indices."
                )
            return positions  # np.ndarray[int]
    else:
        raise IndexError(f"Unknown indexer {indexer!r} of type {type(indexer)}")


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
        return index, slice(None)
    elif len(index) == 2:
        return index
    elif len(index) == 1:
        return index[0], slice(None)
    else:
        raise IndexError("invalid number of indices")


@singledispatch
def _subset(a: np.ndarray | pd.DataFrame, subset_idx: Index):
    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, cabc.Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(DaskArray)
def _subset_dask(a: DaskArray, subset_idx: Index):
    if all(isinstance(x, cabc.Iterable) for x in subset_idx):
        if isinstance(a._meta, csc_matrix):
            return a[:, subset_idx[1]][subset_idx[0], :]
        elif isinstance(a._meta, spmatrix):
            return a[subset_idx[0], :][:, subset_idx[1]]
        else:
            # TODO: this may have been working for some cases?
            subset_idx = np.ix_(*subset_idx)
            return a.vindex[subset_idx]
    return a[subset_idx]


@_subset.register(spmatrix)
def _subset_spmatrix(a: spmatrix, subset_idx: Index):
    # Correcting for indexing behaviour of sparse.spmatrix
    if len(subset_idx) > 1 and all(isinstance(x, cabc.Iterable) for x in subset_idx):
        subset_idx = (subset_idx[0].reshape(-1, 1), *subset_idx[1:])
    return a[subset_idx]


@_subset.register(pd.DataFrame)
def _subset_df(df: pd.DataFrame, subset_idx: Index):
    return df.iloc[subset_idx]


@_subset.register(AwkArray)
def _subset_awkarray(a: AwkArray, subset_idx: Index):
    if all(isinstance(x, cabc.Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


# Registration for SparseDataset occurs in sparse_dataset.py
@_subset.register(h5py.Dataset)
def _subset_dataset(d, subset_idx):
    if not isinstance(subset_idx, tuple):
        subset_idx = (subset_idx,)
    ordered = list(subset_idx)
    rev_order = [slice(None) for _ in range(len(subset_idx))]
    for axis, axis_idx in enumerate(ordered.copy()):
        if isinstance(axis_idx, np.ndarray) and axis_idx.dtype.type != bool:
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
        raise ValueError(
            f"Key {k} could be found in both .{idxdim}_names and .{coldim}.columns"
        )
    elif (in_col + in_idx) == 0:
        raise KeyError(
            f"Could not find key {k} in .{idxdim}_names or .{coldim}.columns."
        )
    elif in_col:
        return getattr(adata, coldim)[k].values
    elif in_idx:
        selected_dim = dims.index(idxdim)
        idx = adata._normalize_indices(make_slice(k, selected_dim))
        a = adata._get_X(layer=layer)[idx]
    if issparse(a):
        a = a.toarray()
    return np.ravel(a)
