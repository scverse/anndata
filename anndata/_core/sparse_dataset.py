"""\
This module implements on disk sparse datasets.

This code is based on and uses the conventions of h5sparse_ by `Appier Inc.`_.
See the copyright and license note in this directory source code.

.. _h5sparse: https://github.com/appier/h5sparse
.. _Appier Inc.: https://www.appier.com/
"""

# TODO:
# - think about supporting the COO format
import collections.abc as cabc
from itertools import accumulate, chain
from typing import Union, NamedTuple, Tuple, Sequence, Iterable, Type
from warnings import warn

import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

try:
    # Not really important, just for IDEs to be more helpful
    from scipy.sparse.compressed import _cs_matrix
except ImportError:
    _cs_matrix = ss.spmatrix

from .index import unpack_index, Index


class BackedFormat(NamedTuple):
    format_str: str
    backed_type: Type["BackedSparseMatrix"]
    memory_type: Type[ss.spmatrix]


class BackedSparseMatrix(_cs_matrix):
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`,
    since that calls copy on `.data`, `.indices`, and `.indptr`.
    """

    def copy(self) -> ss.spmatrix:
        if isinstance(self.data, h5py.Dataset):
            return SparseDataset(self.data.parent).to_memory()
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


class backed_csr_matrix(BackedSparseMatrix, ss.csr_matrix):
    def _get_intXslice(self, row: int, col: slice) -> ss.csr_matrix:
        return ss.csr_matrix(
            get_compressed_vector(self, row), shape=(1, self.shape[1])
        )[:, col]

    def _get_sliceXslice(self, row: slice, col: slice) -> ss.csr_matrix:
        out_shape = (
            slice_len(row, self.shape[0]),
            slice_len(col, self.shape[1]),
        )
        if out_shape[0] == 1:
            return self._get_intXslice(slice_as_int(row, self.shape[0]), col)
        elif out_shape[1] == self.shape[1] and out_shape[0] < self.shape[0]:
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
        out_shape = (
            slice_len(row, self.shape[0]),
            slice_len(col, self.shape[1]),
        )
        if out_shape[1] == 1:
            return self._get_sliceXint(row, slice_as_int(col, self.shape[1]))
        elif out_shape[0] == self.shape[0] and out_shape[1] < self.shape[1]:
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
) -> Tuple[Sequence, Sequence, Sequence]:
    slices = [slice(*(x.indptr[i : i + 2])) for i in row_idxs]
    data = np.concatenate([x.data[s] for s in slices])
    indices = np.concatenate([x.indices[s] for s in slices])
    indptr = list(accumulate(chain((0,), (s.stop - s.start for s in slices))))
    return data, indices, indptr


def get_compressed_vector(
    x: BackedSparseMatrix, idx: int
) -> Tuple[Sequence, Sequence, Sequence]:
    s = slice(*(x.indptr[idx : idx + 2]))
    data = x.data[s]
    indices = x.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


def get_format_str(data: ss.spmatrix) -> str:
    for fmt, _, memory_class in FORMATS:
        if isinstance(data, memory_class):
            return fmt
    raise ValueError(f"Data type {type(data)} is not supported.")


def get_memory_class(format_str: str) -> Type[ss.spmatrix]:
    for fmt, _, memory_class in FORMATS:
        if format_str == fmt:
            return memory_class
    raise ValueError(f"Format string {format_str} is not supported.")


def get_backed_class(format_str: str) -> Type[BackedSparseMatrix]:
    for fmt, backed_class, _ in FORMATS:
        if format_str == fmt:
            return backed_class
    raise ValueError(f"Format string {format_str} is not supported.")


class SparseDataset:
    """Analogous to :class:`h5py.Dataset <h5py:Dataset>`, but for sparse matrices."""

    def __init__(self, group: h5py.Group):
        self.group = group

    @property
    def dtype(self) -> np.dtype:
        return self.group["data"].dtype

    @property
    def format_str(self) -> str:
        if "h5sparse_format" in self.group.attrs:
            return self.group.attrs["h5sparse_format"]
        else:
            # Should this be an extra field?
            return self.group.attrs["encoding-type"].replace("_matrix", "")

    @property
    def h5py_group(self) -> h5py.Group:
        warn(
            "Attribute `h5py_group` of SparseDatasets is deprecated. "
            "Use `group` instead.",
            DeprecationWarning,
        )
        return self.group

    @property
    def name(self) -> str:
        return self.group.name

    @property
    def file(self) -> h5py.File:
        return self.group.file

    @property
    def shape(self) -> Tuple[int, int]:
        shape = self.group.attrs.get("h5sparse_shape")
        return tuple(self.group.attrs["shape"] if shape is None else shape)

    @property
    def value(self) -> ss.spmatrix:
        return self.to_memory()

    def __repr__(self) -> str:
        return (
            f"<HDF5 sparse dataset: format {self.format_str!r}, "
            f"shape {self.shape}, "
            f'type {self.group["data"].dtype.str!r}>'
        )

    def __getitem__(self, index: Union[Index, Tuple[()]]) -> Union[float, ss.spmatrix]:
        row, col = self._normalize_index(index)
        mtx = self.to_backed()
        return mtx[row, col]

    def __setitem__(self, index: Union[Index, Tuple[()]], value):
        row, col = self._normalize_index(index)
        mock_matrix = self.to_backed()
        mock_matrix[row, col] = value

    def _normalize_index(
        self, index: Union[Index, Tuple[()]]
    ) -> Tuple[np.ndarray, np.ndarray]:
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, cabc.Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        return row, col

    def append(self, sparse_matrix: ss.spmatrix):
        # Prep variables
        shape = self.shape
        if isinstance(sparse_matrix, SparseDataset):
            sparse_matrix = sparse_matrix.to_backed()

        # Check input
        if not ss.isspmatrix(sparse_matrix):
            raise NotImplementedError(
                "Currently, only sparse matrices of equivalent format can be "
                "appended to a SparseDataset."
            )
        if self.format_str not in {"csr", "csc"}:
            raise NotImplementedError(
                f"The append method for format {self.format_str} "
                f"is not implemented."
            )
        if self.format_str != get_format_str(sparse_matrix):
            raise ValueError(
                f"Matrices must have same format. Currently are "
                f"{self.format_str!r} and {get_format_str(sparse_matrix)!r}"
            )

        # shape
        if self.format_str == "csr":
            assert (
                shape[1] == sparse_matrix.shape[1]
            ), "CSR matrices must have same size of dimension 1 to be appended."
            new_shape = (shape[0] + sparse_matrix.shape[0], shape[1])
        elif self.format_str == "csc":
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
        format_class = get_backed_class(self.format_str)
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group["data"]
        mtx.indices = self.group["indices"]
        mtx.indptr = self.group["indptr"]
        return mtx

    def to_memory(self) -> ss.spmatrix:
        format_class = get_memory_class(self.format_str)
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group["data"][...]
        mtx.indices = self.group["indices"][...]
        mtx.indptr = self.group["indptr"][...]
        return mtx
