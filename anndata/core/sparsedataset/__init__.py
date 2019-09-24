"""This module implements on disk sparse datasets.

This code was originally based on and uses the conventions of
`h5sparse <https://github.com/appier/h5sparse>`_ by
`Appier Inc. <https://www.appier.com/>`_.
See the copyright and license note in this directory source code.
"""

# TODO:
# - think about supporting the COO format
from collections.abc import Iterable
from itertools import accumulate, chain
from typing import Union, NamedTuple, Tuple, Sequence
from warnings import warn

import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

from ...utils import unpack_index


class BackedFormat(NamedTuple):
    format_str: str
    backed_type: type
    memory_type: type


class BackedSparseMatrixMixin:
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`, since that calls
    copy on `.data`, `.indices`, and `.indptr`.
    """

    def copy(self):
        if isinstance(self.data, h5py.Dataset):
            return SparseDataset(self.data.parent).tomemory()
        else:
            return super().copy()


class backed_csr_matrix(BackedSparseMatrixMixin, ss.csr_matrix):
    pass


class backed_csc_matrix(BackedSparseMatrixMixin, ss.csc_matrix):
    pass


backed_sparse_matrix = Union[backed_csc_matrix, backed_csr_matrix]


FORMATS = [
    BackedFormat("csr", backed_csr_matrix, ss.csr_matrix),
    BackedFormat("csc", backed_csc_matrix, ss.csc_matrix),
]


def get_format_str(data):
    for fmt, backed_class, memory_class in FORMATS:
        if isinstance(data, memory_class):
            return fmt
    raise ValueError(f"Data type {type(data)} is not supported.")


def get_memory_class(format_str):
    for fmt, backed_class, memory_class in FORMATS:
        if format_str == fmt:
            return memory_class
    raise ValueError(f"Format string {format_str} is not supported.")


def get_backed_class(format_str):
    for fmt, backed_class, memory_class in FORMATS:
        if format_str == fmt:
            return backed_class
    raise ValueError(f"Format string {format_str} is not supported.")


def _set_many(self, i, j, x):
    """\
    Sets value at each (i, j) to x

    Here (i,j) index major and minor respectively, and must not contain
    duplicate entries.
    """
    i, j, M, N = self._prepare_indices(i, j)

    if np.isscalar(x):  # Scipy 1.3+ compat
        n_samples = 1
    else:
        n_samples = len(x)
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

    if -1 not in offsets:
        # make a list for interaction with h5py
        offsets = list(offsets)
        # only affects existing non-zero cells
        self.data[offsets] = x
        return

    else:
        raise ValueError(
            'Currently, you cannot change the sparsity structure of a SparseDataset.'
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


backed_csr_matrix._set_many = _set_many
backed_csc_matrix._set_many = _set_many


def _zero_many(self, i, j):
    """\
    Sets value at each (i, j) to zero, preserving sparsity structure.

    Here (i,j) index major and minor respectively.
    """
    i, j, M, N = self._prepare_indices(i, j)

    n_samples = len(i)
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

    # only assign zeros to the existing sparsity structure
    self.data[list(offsets[offsets > -1])] = 0


backed_csr_matrix._zero_many = _zero_many
backed_csc_matrix._zero_many = _zero_many


def get_compressed_vectors(X: backed_sparse_matrix, row_idxs: "np.ndarray[int]") -> Tuple[Sequence, Sequence, Sequence]:
    slices = [slice(*(X.indptr[i:i+2])) for i in row_idxs]
    data = np.concatenate([X.data[s] for s in slices])
    indices = np.concatenate([X.indices[s] for s in slices])
    indptr = list(accumulate(chain((0,), (s.stop - s.start for s in slices))))
    return data, indices, indptr


def get_csc_cols(X: backed_csc_matrix, idxs: np.ndarray) -> ss.csc_matrix:
    idxs = np.asarray(idxs)
    if idxs.dtype == bool:
        idxs = np.where(idxs)
    return ss.csc_matrix(
        get_compressed_vectors(X, idxs),
        shape=(X.shape[0], len(idxs))
    )


def get_csr_rows(X: backed_csr_matrix, idxs: np.ndarray) -> ss.csr_matrix:
    idxs = np.asarray(idxs)
    if idxs.dtype == bool:
        idxs = np.where(idxs)
    return ss.csr_matrix(
        get_compressed_vectors(X, idxs),
        shape=(len(idxs), X.shape[1])
    )


def get_compressed_vector(X, idx: int) -> Tuple[Sequence, Sequence, Sequence]:
    s = slice(*(X.indptr[idx:idx+2]))
    data = X.data[s]
    indices = X.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


def get_csc_col(X: backed_sparse_matrix, idx: int) -> ss.csc_matrix:
    return ss.csc_matrix(
        get_compressed_vector(X, idx),
        shape=(X.shape[0], 1)
    )


def get_csr_row(X, idx: int) -> ss.csr_matrix:
    return ss.csr_matrix(
        get_compressed_vector(X, idx),
        shape=(1, X.shape[1])
    )


def csc_get_sliceXint(X, row: slice, col: int):
    return get_csc_col(X, col)[row, :]


def csr_get_intXslice(X, row: int, col: slice):
    return get_csr_row(X, row)[:, col]


backed_csc_matrix._get_sliceXint = csc_get_sliceXint
backed_csr_matrix._get_intXslice = csr_get_intXslice

backed_csc_matrix._get_sliceXarray = lambda x, row, col: get_csc_cols(x, col)[row, :]
backed_csr_matrix._get_arrayXslice = lambda x, row, col: get_csr_rows(x, row)[:, col]


class SparseDataset:
    """\
    Analogous to :class:`h5py.Dataset <h5py:Dataset>`, but for sparse matrices.
    """

    def __init__(self, group):
        self.group = group

    @property
    def dtype(self):
        return self.group['data'].dtype

    @property
    def format_str(self):
        if "h5sparse_format" in self.group.attrs:
            return self.group.attrs['h5sparse_format']
        else:
            return self.group.attrs["encoding-type"].replace("_matrix", "")  # Should this be an extra field?

    @property
    def h5py_group(self):
        warn(DeprecationWarning("Attribute `h5py_group` of SparseDatasets is deprecated. Use `group` instead."))
        return self.group

    @property
    def name(self):
        return self.group.name

    @property
    def shape(self):
        if "h5sparse_shape" in self.group.attrs:
            return tuple(self.group.attrs['h5sparse_shape'])
        else:
            return tuple(self.group.attrs["shape"])

    @property
    def value(self):
        return self.tomemory()

    def __repr__(self):
        return (
            f'<HDF5 sparse dataset: format {self.format_str!r}, '
            f'shape {self.shape}, '
            f'type {self.group["data"].dtype.str!r}>'
        )

    def __getitem__(self, index):
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        mtx = self.tobacked()
        return mtx[row, col]

    def __setitem__(self, index, value):
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        mock_matrix = self.tobacked()
        mock_matrix[row, col] = value

    def append(self, sparse_matrix):
        # Prep variables
        shape = self.shape
        if isinstance(sparse_matrix, SparseDataset):
            sparse_matrix = sparse_matrix.tobacked()

        # Check input
        if self.format_str not in {'csr', 'csc'}:
            raise NotImplementedError(
                f"The append method for format {self.format_str} "
                f"is not implemented."
            )
        if not isinstance(sparse_matrix, ss.spmatrix):
            raise NotImplementedError(
                "Currently, only sparse matrices of equivalent format can be "
                "appended to a SparseDataset."
            )
        if self.format_str != get_format_str(sparse_matrix):
            raise ValueError(
                f"Matrices must have same format. Currently are "
                f"'{self.format_str}' and '{get_format_str(sparse_matrix)}'"
            )

        # data
        data = self.group['data']
        orig_data_size = data.shape[0]
        new_shape = (orig_data_size + sparse_matrix.data.shape[0],)
        data.resize(new_shape)
        data[orig_data_size:] = sparse_matrix.data

        # indptr
        indptr = self.group['indptr']
        orig_data_size = indptr.shape[0]
        append_offset = indptr[-1]
        new_shape = (orig_data_size + sparse_matrix.indptr.shape[0] - 1,)
        indptr.resize(new_shape)
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + append_offset
        )

        # indices
        indices = self.group['indices']
        orig_data_size = indices.shape[0]
        new_shape = (orig_data_size + sparse_matrix.indices.shape[0],)
        indices.resize(new_shape)
        indices[orig_data_size:] = sparse_matrix.indices

        # shape
        if "h5sparse_shape" in self.group.attrs:
            del self.group.attrs["h5sparse_shape"]

        # TODO: Do we want to allow different sizes on the unaligned axis?
        if self.format_str == "csr":
            new_shape = (
                shape[0] + sparse_matrix.shape[0],
                max(shape[1], sparse_matrix.shape[1]),
            )
        elif self.format_str == "csc":
            new_shape = (
                max(shape[0], sparse_matrix.shape[0]),
                shape[1] + sparse_matrix.shape[1],
            )
        self.group.attrs['shape'] = new_shape

    def tobacked(self):
        format_class = get_backed_class(self.format_str)
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group['data']
        mtx.indices = self.group['indices']
        mtx.indptr = self.group['indptr']
        return mtx

    def tomemory(self):
        format_class = get_memory_class(self.format_str)
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group['data'][...]
        mtx.indices = self.group['indices'][...]
        mtx.indptr = self.group['indptr'][...]
        return mtx
