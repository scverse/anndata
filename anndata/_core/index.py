import collections.abc as cabc
from functools import singledispatch
from typing import Union, Sequence, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


Index1D = Union[slice, int, str, np.int64, np.ndarray]
Index = Union[Index1D, Tuple[Index1D, Index1D], spmatrix]


def _normalize_indices(
    index: Optional[Index], names0: pd.Index, names1: pd.Index
) -> Tuple[slice, slice]:
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
    indexer: Union[
        slice,
        np.integer,
        int,
        str,
        Sequence[Union[int, np.integer]],
        np.ndarray,
        pd.Index,
    ],
    index: pd.Index,
) -> Union[slice, int, np.ndarray]:  # ndarray of int
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
    elif isinstance(indexer, (Sequence, np.ndarray, pd.Index)):
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


def unpack_index(index: Index) -> Tuple[Index1D, Index1D]:
    # handle indexing with boolean matrices
    if (
        isinstance(index, (spmatrix, np.ndarray))
        and index.ndim == 2
        and index.dtype.kind == "b"
    ):
        return index.nonzero()

    if not isinstance(index, tuple):
        return index, slice(None)
    elif len(index) == 2:
        return index
    elif len(index) == 1:
        return index[0], slice(None)
    else:
        raise IndexError("invalid number of indices")


@singledispatch
def _subset(a: Union[np.ndarray, spmatrix, pd.DataFrame], subset_idx: Index):
    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, cabc.Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(pd.DataFrame)
def _subset_df(df: pd.DataFrame, subset_idx: Index):
    return df.iloc[subset_idx]
