"""\
This module implements on disk sparse datasets.

This code is based on and uses the conventions of h5sparse_ by `Appier Inc.`_.
See the copyright and license note in this directory source code.

.. _h5sparse: https://github.com/appier/h5sparse
.. _Appier Inc.: https://www.appier.com/
"""

# TODO:
# - think about supporting the COO format
from __future__ import annotations

import collections.abc as cabc
import warnings
from abc import ABC
from itertools import accumulate, chain
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Literal,
    NamedTuple,
)

import h5py
import numpy as np
import scipy.sparse as ss
import zarr
from scipy.sparse import _sparsetools

from anndata._core.index import _fix_slice_bounds
from anndata._core.views import _resolve_idx, as_view
from anndata.compat import H5Group, ZarrGroup

from ..compat import ZarrArray, _read_attr

try:
    # Not really important, just for IDEs to be more helpful
    from scipy.sparse.compressed import _cs_matrix
except ImportError:
    _cs_matrix = ss.spmatrix

from .index import Index, _subset, unpack_index

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence


class BackedFormat(NamedTuple):
    format: str
    backed_type: type[BackedSparseMatrix]
    memory_type: type[ss.spmatrix]


class BackedSparseMatrix(_cs_matrix):
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`,
    since that calls copy on `.data`, `.indices`, and `.indptr`.
    """

    _cached_indptr = None

    def copy(self) -> ss.spmatrix:
        if isinstance(self.data, h5py.Dataset):
            return sparse_dataset(self.data.parent).to_memory()
        if isinstance(self.data, ZarrArray):
            return sparse_dataset(
                zarr.open(
                    store=self.data.store, path=Path(self.data.path).parent, mode="r"
                )
            ).to_memory()
        else:
            return super().copy()

    def _set_many(self, i: Iterable[int], j: Iterable[int], x):
        """\
        Sets value at each (i, j) to x

        Here (i,j) index major and minor respectively,
        and must not contain duplicate entries.
        """
        # Scipy 1.3+ compat
        n_samples = 1 if np.isscalar(x) else len(x)
        offsets = self._offsets(i, j, n_samples)

        if -1 not in offsets:
            # make a list for interaction with h5py
            offsets = list(offsets)
            # only affects existing non-zero cells
            self.data[offsets] = x
            return

        else:
            raise ValueError(
                "You cannot change the sparsity structure of a SparseDataset."
            )
            # replace where possible
            # mask = offsets > -1
            # # offsets[mask]
            # bool_data_mask = np.zeros(len(self.data), dtype=bool)
            # bool_data_mask[offsets[mask]] = True
            # self.data[bool_data_mask] = x[mask]
            # # self.data[offsets[mask]] = x[mask]
            # # only insertions remain
            # mask = ~mask
            # i = i[mask]
            # i[i < 0] += M
            # j = j[mask]
            # j[j < 0] += N
            # self._insert_many(i, j, x[mask])

    def _zero_many(self, i: Sequence[int], j: Sequence[int]):
        """\
        Sets value at each (i, j) to zero, preserving sparsity structure.

        Here (i,j) index major and minor respectively.
        """
        offsets = self._offsets(i, j, len(i))

        # only assign zeros to the existing sparsity structure
        self.data[list(offsets[offsets > -1])] = 0

    def _offsets(
        self, i: Iterable[int], j: Iterable[int], n_samples: int
    ) -> np.ndarray:
        i, j, M, N = self._prepare_indices(i, j)
        offsets = np.empty(n_samples, dtype=self.indices.dtype)
        ret = _sparsetools.csr_sample_offsets(
            M, N, self.indptr, self.indices, n_samples, i, j, offsets
        )
        if ret == 1:
            # rinse and repeat
            self.sum_duplicates()
            _sparsetools.csr_sample_offsets(
                M, N, self.indptr, self.indices, n_samples, i, j, offsets
            )
        return offsets

    @property
    def indptr(self):
        if self._cached_indptr is None:
            self._cached_indptr = self._indptr[:]
        return self._cached_indptr

    @indptr.setter
    def indptr(self, indptr):
        self._indptr = indptr
        self._cached_indptr = None

    def _get_contiguous_compressed_slice(
        self, s: slice
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        new_indptr = self.indptr[s.start : s.stop + 1]

        start = new_indptr[0]
        stop = new_indptr[-1]

        new_indptr -= start

        new_data = self.data[start:stop]
        new_indices = self.indices[start:stop]

        return new_data, new_indices, new_indptr


class backed_csr_matrix(BackedSparseMatrix, ss.csr_matrix):
    def _get_intXslice(self, row: int, col: slice) -> ss.csr_matrix:
        return ss.csr_matrix(
            get_compressed_vector(self, row), shape=(1, self.shape[1])
        )[:, col]

    def _get_sliceXslice(self, row: slice, col: slice) -> ss.csr_matrix:
        row = _fix_slice_bounds(row, self.shape[0])
        col = _fix_slice_bounds(col, self.shape[1])

        out_shape = (
            slice_len(row, self.shape[0]),
            slice_len(col, self.shape[1]),
        )

        if out_shape[0] == 1:
            return self._get_intXslice(slice_as_int(row, self.shape[0]), col)
        elif out_shape[1] == self.shape[1] and out_shape[0] < self.shape[0]:
            if row.step == 1:
                return ss.csr_matrix(
                    self._get_contiguous_compressed_slice(row), shape=out_shape
                )
            return self._get_arrayXslice(np.arange(*row.indices(self.shape[0])), col)
        return super()._get_sliceXslice(row, col)

    def _get_arrayXslice(self, row: Sequence[int], col: slice) -> ss.csr_matrix:
        idxs = np.asarray(row)
        if idxs.dtype == bool:
            idxs = np.where(idxs)
        return ss.csr_matrix(
            get_compressed_vectors(self, idxs), shape=(len(idxs), self.shape[1])
        )[:, col]


class backed_csc_matrix(BackedSparseMatrix, ss.csc_matrix):
    def _get_sliceXint(self, row: slice, col: int) -> ss.csc_matrix:
        return ss.csc_matrix(
            get_compressed_vector(self, col), shape=(self.shape[0], 1)
        )[row, :]

    def _get_sliceXslice(self, row: slice, col: slice) -> ss.csc_matrix:
        row = _fix_slice_bounds(row, self.shape[0])
        col = _fix_slice_bounds(col, self.shape[1])

        out_shape = (
            slice_len(row, self.shape[0]),
            slice_len(col, self.shape[1]),
        )

        if out_shape[1] == 1:
            return self._get_sliceXint(row, slice_as_int(col, self.shape[1]))
        elif out_shape[0] == self.shape[0] and out_shape[1] < self.shape[1]:
            if col.step == 1:
                return ss.csc_matrix(
                    self._get_contiguous_compressed_slice(col), shape=out_shape
                )
            return self._get_sliceXarray(row, np.arange(*col.indices(self.shape[1])))
        return super()._get_sliceXslice(row, col)

    def _get_sliceXarray(self, row: slice, col: Sequence[int]) -> ss.csc_matrix:
        idxs = np.asarray(col)
        if idxs.dtype == bool:
            idxs = np.where(idxs)
        return ss.csc_matrix(
            get_compressed_vectors(self, idxs), shape=(self.shape[0], len(idxs))
        )[row, :]


FORMATS = [
    BackedFormat("csr", backed_csr_matrix, ss.csr_matrix),
    BackedFormat("csc", backed_csc_matrix, ss.csc_matrix),
]


def slice_len(s: slice, l: int) -> int:
    """Returns length of `a[s]` where `len(a) == l`."""
    return len(range(*s.indices(l)))


def slice_as_int(s: slice, l: int) -> int:
    """Converts slices of length 1 to the integer index they’ll access."""
    out = list(range(*s.indices(l)))
    assert len(out) == 1
    return out[0]


def get_compressed_vectors(
    x: BackedSparseMatrix, row_idxs: Iterable[int]
) -> tuple[Sequence, Sequence, Sequence]:
    slices = [slice(*(x.indptr[i : i + 2])) for i in row_idxs]
    data = np.concatenate([x.data[s] for s in slices])
    indices = np.concatenate([x.indices[s] for s in slices])
    indptr = list(accumulate(chain((0,), (s.stop - s.start for s in slices))))
    return data, indices, indptr


def get_compressed_vector(
    x: BackedSparseMatrix, idx: int
) -> tuple[Sequence, Sequence, Sequence]:
    s = slice(*(x.indptr[idx : idx + 2]))
    data = x.data[s]
    indices = x.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


def get_format(data: ss.spmatrix) -> str:
    for fmt, _, memory_class in FORMATS:
        if isinstance(data, memory_class):
            return fmt
    raise ValueError(f"Data type {type(data)} is not supported.")


def get_memory_class(format: str) -> type[ss.spmatrix]:
    for fmt, _, memory_class in FORMATS:
        if format == fmt:
            return memory_class
    raise ValueError(f"Format string {format} is not supported.")


def get_backed_class(format: str) -> type[BackedSparseMatrix]:
    for fmt, backed_class, _ in FORMATS:
        if format == fmt:
            return backed_class
    raise ValueError(f"Format string {format} is not supported.")


def _get_group_format(group) -> str:
    if "h5sparse_format" in group.attrs:
        # TODO: Warn about an old format
        # If this is only just going to be public, I could insist it's not like this
        return _read_attr(group.attrs, "h5sparse_format")
    else:
        # Should this be an extra field?
        return _read_attr(group.attrs, "encoding-type").replace("_matrix", "")


class BaseCompressedSparseDataset(ABC):
    """Analogous to :class:`h5py.Dataset <h5py:Dataset>` or `zarr.Array`, but for sparse matrices."""

    def __init__(self, group: h5py.Group | ZarrGroup):
        type(self)._check_group_format(group)
        self.group = group
        self._row_subset_idx = slice(None, None, None)
        self._col_subset_idx = slice(None, None, None)

    shape: tuple[int, int]
    """Shape of the matrix."""

    @property
    def backend(self) -> Literal["zarr", "hdf5"]:
        if isinstance(self.group, ZarrGroup):
            return "zarr"
        elif isinstance(self.group, H5Group):
            return "hdf5"
        else:
            raise ValueError(f"Unknown group type {type(self.group)}")

    @property
    def row_subset_idx(self):
        if isinstance(self._row_subset_idx, np.ndarray):
            return self._row_subset_idx.flatten()  # why????
        return self._row_subset_idx

    @property
    def has_no_subset_idx(self):
        return self.has_no_col_subset_idx and self.has_no_row_subset_idx

    @property
    def has_no_col_subset_idx(self):
        if isinstance(self.col_subset_idx, slice):
            if self.col_subset_idx == slice(
                None, None, None
            ) or self.col_subset_idx == slice(0, self.get_backing_shape()[1], 1):
                return True
        return False

    @property
    def has_no_row_subset_idx(self):
        if isinstance(self.row_subset_idx, slice):
            if self.row_subset_idx == slice(
                None, None, None
            ) or self.row_subset_idx == slice(0, self.get_backing_shape()[0], 1):
                return True
        return False

    @row_subset_idx.setter
    def row_subset_idx(self, new_idx):
        self._row_subset_idx = (
            new_idx
            if self.row_subset_idx is None
            else _resolve_idx(self.row_subset_idx, new_idx, self.shape[0])
        )

    @property
    def col_subset_idx(self):
        if isinstance(self._col_subset_idx, np.ndarray):
            return self._col_subset_idx.flatten()
        return self._col_subset_idx

    @col_subset_idx.setter
    def col_subset_idx(self, new_idx):
        self._col_subset_idx = (
            new_idx
            if self.col_subset_idx is None
            else _resolve_idx(self.col_subset_idx, new_idx, self.shape[1])
        )

    @property
    def dtype(self) -> np.dtype:
        return self.group["data"].dtype

    @classmethod
    def _check_group_format(cls, group):
        group_format = _get_group_format(group)
        assert group_format == cls.format

    @property
    def format_str(self) -> Literal["csc", "csr"]:
        """DEPRECATED Use .format instead."""
        warnings.warn(
            "The attribute .format_str is deprecated and will be removed in the anndata 0.11.0. "
            "Please use .format instead.",
            FutureWarning,
        )
        return self.format

    @property
    def name(self) -> str:
        return self.group.name

    def get_backing_shape(self) -> tuple[int, int]:
        shape = _read_attr(self.group.attrs, "shape", None)
        if shape is None:
            # TODO warn
            shape = self.group.attrs.get("h5sparse_shape")
        return tuple(shape)

    @property
    def shape(self) -> tuple[int, int]:
        shape = self.get_backing_shape()
        if self.has_no_subset_idx:
            return tuple(shape)
        row_length = 0
        col_length = 0
        if isinstance(self.row_subset_idx, slice):
            if self.row_subset_idx == slice(None, None, None):
                row_length = shape[0]
            else:
                row_length = self.row_subset_idx.stop - self.row_subset_idx.start
        elif np.isscalar(self.row_subset_idx):
            row_length = 1
        else:
            row_length = len(self.row_subset_idx)  # can we assume a flatten method?
        if isinstance(self.col_subset_idx, slice):
            if self.col_subset_idx == slice(None, None, None):
                col_length = shape[1]
            else:
                col_length = self.col_subset_idx.stop - self.col_subset_idx.start
        elif np.isscalar(self.col_subset_idx):
            row_length = 1
        else:
            col_length = len(self.col_subset_idx)  # can we assume a flatten method?
        return (row_length, col_length)

    @property
    def value(self) -> ss.spmatrix:
        """DEPRECATED Use .to_memory() instead."""
        warnings.warn(
            "The .value attribute is deprecated and will be removed in the anndata 0.11.0. "
            "Please use .to_memory() instead.",
            FutureWarning,
        )
        return self.to_memory()

    def __repr__(self) -> str:
        return f"{type(self).__name__}: backend {self.backend}, shape {self.shape}, data_dtype {self.dtype}"

    def __getitem__(self, index: Index | tuple[()]) -> float | ss.spmatrix:
        row, col = self._normalize_index(index)
        new_mtx = sparse_dataset(self.group)
        new_mtx.row_subset_idx = (
            self.row_subset_idx
        )  # subset idx of subset idx requires original row_subset_idx and new row
        new_mtx.row_subset_idx = row
        new_mtx.col_subset_idx = self.col_subset_idx
        new_mtx.col_subset_idx = col
        return new_mtx

    def _normalize_index(
        self, index: Index | tuple[()]
    ) -> tuple[np.ndarray, np.ndarray]:
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, cabc.Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        return row, col

    def __setitem__(self, index: Index | tuple[()], value):
        warnings.warn(
            "__setitem__ will likely be removed in the near future. We do not recommend relying on its stability.",
            PendingDeprecationWarning,
        )
        row, col = self._normalize_index(index)
        mock_matrix = self.to_backed()
        mock_matrix[row, col] = value

    # TODO: split to other classes?
    def append(self, sparse_matrix: ss.spmatrix):
        # Prep variables
        shape = self.shape
        if isinstance(sparse_matrix, BaseCompressedSparseDataset):
            sparse_matrix = sparse_matrix.to_backed()

        # Check input
        if not ss.isspmatrix(sparse_matrix):
            raise NotImplementedError(
                "Currently, only sparse matrices of equivalent format can be "
                "appended to a SparseDataset."
            )
        if self.format not in {"csr", "csc"}:
            raise NotImplementedError(
                f"The append method for format {self.format} " f"is not implemented."
            )
        if self.format != get_format(sparse_matrix):
            raise ValueError(
                f"Matrices must have same format. Currently are "
                f"{self.format!r} and {get_format(sparse_matrix)!r}"
            )

        # shape
        if self.format == "csr":
            assert (
                shape[1] == sparse_matrix.shape[1]
            ), "CSR matrices must have same size of dimension 1 to be appended."
            new_shape = (shape[0] + sparse_matrix.shape[0], shape[1])
        elif self.format == "csc":
            assert (
                shape[0] == sparse_matrix.shape[0]
            ), "CSC matrices must have same size of dimension 0 to be appended."
            new_shape = (shape[0], shape[1] + sparse_matrix.shape[1])
        else:
            assert False, "We forgot to update this branching to a new format"
        if "h5sparse_shape" in self.group.attrs:
            del self.group.attrs["h5sparse_shape"]
        self.group.attrs["shape"] = new_shape

        # data
        data = self.group["data"]
        orig_data_size = data.shape[0]
        data.resize((orig_data_size + sparse_matrix.data.shape[0],))
        data[orig_data_size:] = sparse_matrix.data

        # indptr
        indptr = self.group["indptr"]
        orig_data_size = indptr.shape[0]
        append_offset = indptr[-1]
        indptr.resize((orig_data_size + sparse_matrix.indptr.shape[0] - 1,))
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + append_offset
        )

        # indices
        indices = self.group["indices"]
        orig_data_size = indices.shape[0]
        indices.resize((orig_data_size + sparse_matrix.indices.shape[0],))
        indices[orig_data_size:] = sparse_matrix.indices

    def to_backed(self) -> BackedSparseMatrix:
        format_class = get_backed_class(self.format)
        mtx = format_class(self.get_backing_shape(), dtype=self.dtype)
        mtx.data = self.group["data"]
        mtx.indices = self.group["indices"]
        mtx.indptr = self.group["indptr"]
        return mtx

    def to_memory(self) -> ss.spmatrix:
        # Could not get row idx with csc and vice versa working without reading into memory but shouldn't matter
        if (self.format == "csr" and self.has_no_row_subset_idx) or (
            self.format == "csc" and self.has_no_col_subset_idx
        ):
            format_class = get_memory_class(self.format)
            mtx = format_class(self.get_backing_shape(), dtype=self.dtype)
            mtx.data = self.group["data"][...]
            mtx.indices = self.group["indices"][...]
            mtx.indptr = self.group["indptr"][...]
            if self.has_no_subset_idx:
                return mtx
        else:
            mtx = self.to_backed()
        if self.format == "csr":
            return mtx[self.row_subset_idx, :][:, self.col_subset_idx]
        return mtx[:, self.col_subset_idx][self.row_subset_idx, :]

    def toarray(self) -> np.ndarray:
        return self.to_memory().toarray()


_sparse_dataset_doc = """\
    On disk {format} sparse matrix.

    Parameters
    ----------
    group
        The backing group store.
"""


class CSRDataset(BaseCompressedSparseDataset):
    __doc__ = _sparse_dataset_doc.format(
        format="CSR",
    )
    format = "csr"


class CSCDataset(BaseCompressedSparseDataset):
    __doc__ = _sparse_dataset_doc.format(
        format="CSC",
    )
    format = "csc"


def sparse_dataset(group: ZarrGroup | H5Group) -> CSRDataset | CSCDataset:
    """Generates a backed mode-compatible sparse dataset class.

    Parameters
    ----------
    group
        The backing group store.

    Returns
    -------
        Sparse dataset class.

    Example
    -------

    >>> import zarr
    >>> from anndata.experimental import sparse_dataset
    >>> group = zarr.open_group('./my_test_store.zarr')
    >>> group['data'] = [10, 20, 30, 40, 50, 60, 70, 80]
    >>> group['indices'] = [0, 1, 1, 3, 2, 3, 4, 5]
    >>> group['indptr'] = [0, 2, 4, 7, 8]
    >>> group.attrs['shape'] = (4, 6)
    >>> group.attrs['encoding-type'] = 'csr_matrix'
    >>> sparse_dataset(group)
    CSRDataset: backend zarr, shape (4, 6), data_dtype int64
    """
    encoding_type = _get_group_format(group)
    if encoding_type == "csr":
        return CSRDataset(group)
    elif encoding_type == "csc":
        return CSCDataset(group)


@_subset.register(BaseCompressedSparseDataset)
def subset_sparsedataset(d, subset_idx):
    return d[subset_idx]


@as_view.register(BaseCompressedSparseDataset)
def _view_masked(a: BaseCompressedSparseDataset, view_args):
    return a


## Backwards compat

_sparsedataset_depr_msg = """\
SparseDataset is deprecated and will be removed in late 2024. It has been replaced by the public classes CSRDataset and CSCDataset.

For instance checks, use `isinstance(X, (anndata.experimental.CSRDataset, anndata.experimental.CSCDataset))` instead.

For creation, use `anndata.experimental.sparse_dataset(X)` instead.
"""


class SparseDataset(ABC):
    """DEPRECATED.

    Use CSRDataset, CSCDataset, and sparse_dataset from anndata.experimental instead.
    """

    def __new__(cls, group):
        warnings.warn(FutureWarning(_sparsedataset_depr_msg), stacklevel=2)
        return sparse_dataset(group)

    @classmethod
    def __subclasshook__(cls, C):
        warnings.warn(FutureWarning(_sparsedataset_depr_msg), stacklevel=3)
        return issubclass(C, (CSRDataset, CSCDataset))
