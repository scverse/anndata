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

import warnings
from abc import ABC
from collections.abc import Iterable
from functools import cached_property
from itertools import accumulate, chain
from math import floor
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

from .. import abc
from .._settings import settings
from ..compat import H5Group, SpArray, ZarrArray, ZarrGroup, _read_attr
from .index import _fix_slice_bounds, _subset, unpack_index

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from scipy.sparse._compressed import _cs_matrix

    from .._types import GroupStorageType
    from .index import Index
else:
    from scipy.sparse import spmatrix as _cs_matrix


class BackedFormat(NamedTuple):
    format: Literal["csr", "csc"]
    backed_type: type[BackedSparseMatrix]
    memory_type: type[_cs_matrix]


class BackedSparseMatrix(_cs_matrix):
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`,
    since that calls copy on `.data`, `.indices`, and `.indptr`.
    """

    data: GroupStorageType
    indices: GroupStorageType
    indptr: np.ndarray

    def copy(self) -> ss.csr_matrix | ss.csc_matrix:
        if isinstance(self.data, h5py.Dataset):
            return sparse_dataset(self.data.parent).to_memory()
        if isinstance(self.data, ZarrArray):
            import zarr

            return sparse_dataset(
                zarr.open(
                    store=self.data.store,
                    mode="r",
                    chunk_store=self.data.chunk_store,  # chunk_store is needed, not clear why
                )[Path(self.data.path).parent]
            ).to_memory()
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

    def _get_contiguous_compressed_slice(
        self, s: slice
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        new_indptr = self.indptr[s.start : s.stop + 1].copy()

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
        if row.step != 1:
            return self._get_arrayXslice(np.arange(*row.indices(self.shape[0])), col)
        res = ss.csr_matrix(
            self._get_contiguous_compressed_slice(row),
            shape=(out_shape[0], self.shape[1]),
        )
        return res if out_shape[1] == self.shape[1] else res[:, col]

    def _get_arrayXslice(self, row: Sequence[int], col: slice) -> ss.csr_matrix:
        idxs = np.asarray(row)
        if len(idxs) == 0:
            return ss.csr_matrix((0, self.shape[1]))
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
        if col.step != 1:
            return self._get_sliceXarray(row, np.arange(*col.indices(self.shape[1])))
        res = ss.csc_matrix(
            self._get_contiguous_compressed_slice(col),
            shape=(self.shape[0], out_shape[1]),
        )
        return res if out_shape[0] == self.shape[0] else res[row, :]

    def _get_sliceXarray(self, row: slice, col: Sequence[int]) -> ss.csc_matrix:
        idxs = np.asarray(col)
        if len(idxs) == 0:
            return ss.csc_matrix((self.shape[0], 0))
        if idxs.dtype == bool:
            idxs = np.where(idxs)
        return ss.csc_matrix(
            get_compressed_vectors(self, idxs), shape=(self.shape[0], len(idxs))
        )[row, :]


FORMATS = [
    BackedFormat("csr", backed_csr_matrix, ss.csr_matrix),
    BackedFormat("csc", backed_csc_matrix, ss.csc_matrix),
    BackedFormat("csr", backed_csr_matrix, ss.csr_array),
    BackedFormat("csc", backed_csc_matrix, ss.csc_array),
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


def get_compressed_vectors_for_slices(
    x: BackedSparseMatrix, slices: Iterable[slice]
) -> tuple[Sequence, Sequence, Sequence]:
    indptr_sels = [x.indptr[slice(s.start, s.stop + 1)] for s in slices]
    data = np.concatenate([x.data[s[0] : s[-1]] for s in indptr_sels])
    indices = np.concatenate([x.indices[s[0] : s[-1]] for s in indptr_sels])
    # Need to track the size of the gaps in the slices to each indptr subselection
    total = indptr_sels[0][0]
    offsets = [total]
    for i, sel in enumerate(indptr_sels[1:]):
        total = (sel[0] - indptr_sels[i][-1]) + total
        offsets.append(total)
    start_indptr = indptr_sels[0] - offsets[0]
    if len(slices) < 2:  # there is only one slice so no need to concatenate
        return data, indices, start_indptr
    end_indptr = np.concatenate(
        [s[1:] - offsets[i + 1] for i, s in enumerate(indptr_sels[1:])]
    )
    indptr = np.concatenate([start_indptr, end_indptr])
    return data, indices, indptr


def get_compressed_vector(
    x: BackedSparseMatrix, idx: int
) -> tuple[Sequence, Sequence, Sequence]:
    s = slice(*(x.indptr[idx : idx + 2]))
    data = x.data[s]
    indices = x.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


def subset_by_major_axis_mask(
    mtx: _cs_matrix, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    slices = np.ma.extras._ezclump(mask)

    def mean_slice_length(slices):
        return floor(sum(s.stop - s.start for s in slices) / len(slices))

    # heuristic for whether slicing should be optimized
    if len(slices) > 0:
        if mean_slice_length(slices) <= 7:
            return get_compressed_vectors(mtx, np.where(mask)[0])
        else:
            return get_compressed_vectors_for_slices(mtx, slices)
    return [], [], [0]


def get_memory_class(
    format: Literal["csr", "csc"], *, use_sparray_in_io: bool = False
) -> type[_cs_matrix]:
    for fmt, _, memory_class in FORMATS:
        if format == fmt:
            if use_sparray_in_io and issubclass(memory_class, SpArray):
                return memory_class
            elif not use_sparray_in_io and issubclass(memory_class, ss.spmatrix):
                return memory_class
    raise ValueError(f"Format string {format} is not supported.")


def get_backed_class(
    format: Literal["csr", "csc"], *, use_sparray_in_io: bool = False
) -> type[BackedSparseMatrix]:
    for fmt, backed_class, _ in FORMATS:
        if format == fmt:
            if use_sparray_in_io and issubclass(backed_class, SpArray):
                return backed_class
            elif not use_sparray_in_io and issubclass(backed_class, ss.spmatrix):
                return backed_class
    raise ValueError(f"Format string {format} is not supported.")


def _get_group_format(group: GroupStorageType) -> str:
    if "h5sparse_format" in group.attrs:
        # TODO: Warn about an old format
        # If this is only just going to be public, I could insist it's not like this
        return _read_attr(group.attrs, "h5sparse_format")
    else:
        # Should this be an extra field?
        return _read_attr(group.attrs, "encoding-type").replace("_matrix", "")


# Check for the overridden few methods above in our BackedSparseMatrix subclasses
def is_sparse_indexing_overridden(format: Literal["csr", "csc"], row, col):
    major_indexer, minor_indexer = (row, col) if format == "csr" else (col, row)
    return isinstance(minor_indexer, slice) and (
        (isinstance(major_indexer, int | np.integer))
        or (isinstance(major_indexer, slice))
        or (isinstance(major_indexer, np.ndarray) and major_indexer.ndim == 1)
    )


class BaseCompressedSparseDataset(abc._AbstractCSDataset, ABC):
    _group: GroupStorageType

    def __init__(self, group: GroupStorageType):
        type(self)._check_group_format(group)
        self._group = group

    @property
    def group(self) -> GroupStorageType:
        """The group underlying the backed matrix."""
        return self._group

    @group.setter
    def group(self, val):
        raise AttributeError(
            f"Do not reset group on a {type(self)} with {val}.  Instead use `sparse_dataset` to make a new class."
        )

    @property
    def backend(self) -> Literal["zarr", "hdf5"]:
        """Which file type is used on-disk."""
        if isinstance(self.group, ZarrGroup):
            return "zarr"
        elif isinstance(self.group, H5Group):
            return "hdf5"
        else:
            raise ValueError(f"Unknown group type {type(self.group)}")

    @property
    def dtype(self) -> np.dtype:
        """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""
        return self.group["data"].dtype

    @classmethod
    def _check_group_format(cls, group):
        group_format = _get_group_format(group)
        assert group_format == cls.format

    @property
    def _name(self) -> str:
        """Name of the group."""
        return self.group.name

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the matrix read off disk."""
        shape = _read_attr(self.group.attrs, "shape", None)
        if shape is None:
            # TODO warn
            shape = self.group.attrs.get("h5sparse_shape")
        return tuple(map(int, shape))

    def __repr__(self) -> str:
        name = type(self).__name__.removeprefix("_")
        return f"{name}: backend {self.backend}, shape {self.shape}, data_dtype {self.dtype}"

    def __getitem__(
        self, index: Index | tuple[()]
    ) -> float | ss.csr_matrix | ss.csc_matrix | SpArray:
        indices = self._normalize_index(index)
        row, col = indices
        mtx = self._to_backed()
        row_sp_matrix_validated, col_sp_matrix_validated = mtx._validate_indices(
            (row, col)
        )

        # Handle masked indexing along major axis
        if self.format == "csr" and np.array(row).dtype == bool:
            sub = ss.csr_matrix(
                subset_by_major_axis_mask(mtx, row), shape=(row.sum(), mtx.shape[1])
            )[:, col]
        elif self.format == "csc" and np.array(col).dtype == bool:
            sub = ss.csc_matrix(
                subset_by_major_axis_mask(mtx, col), shape=(mtx.shape[0], col.sum())
            )[row, :]
        # read into memory data if we do not override access methods
        elif not is_sparse_indexing_overridden(
            self.format, row_sp_matrix_validated, col_sp_matrix_validated
        ):
            sub = self.to_memory()[row_sp_matrix_validated, col_sp_matrix_validated]
        else:
            sub = mtx[row, col]

        # If indexing is array x array it returns a backed_sparse_matrix
        # Not sure what the performance is on that operation
        # Also need to check if memory format is not matrix
        mtx_fmt = get_memory_class(
            self.format, use_sparray_in_io=settings.use_sparse_array_on_read
        )
        must_convert_to_array = issubclass(mtx_fmt, SpArray) and not isinstance(
            sub, SpArray
        )
        if isinstance(sub, BackedSparseMatrix) or must_convert_to_array:
            return mtx_fmt(sub)
        else:
            return sub

    def _normalize_index(
        self, index: Index | tuple[()]
    ) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(index, tuple) and not len(index):
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        return row, col

    def __setitem__(self, index: Index | tuple[()], value) -> None:
        warnings.warn(
            "__setitem__ will likely be removed in the near future. We do not recommend relying on its stability.",
            PendingDeprecationWarning,
        )
        row, col = self._normalize_index(index)
        mock_matrix = self._to_backed()
        mock_matrix[row, col] = value

    # TODO: split to other classes?
    def append(self, sparse_matrix: ss.csr_matrix | ss.csc_matrix | SpArray) -> None:
        """Append an in-memory or on-disk sparse matrix to the current object's store.

        Parameters
        ----------
        sparse_matrix
            The matrix to append.

        Raises
        ------
        NotImplementedError
            If the matrix to append is not one of :class:`~scipy.sparse.csr_array`, :class:`~scipy.sparse.csc_array`, :class:`~scipy.sparse.csr_matrix`, or :class:`~scipy.sparse.csc_matrix`.
        ValueError
            If both the on-disk and to-append matrices are not of the same format i.e., `csr` or `csc`.
        OverflowError
            If the underlying data store has a 32 bit indptr, and the new matrix is too large to fit in it i.e., would cause a 64 bit `indptr` to be written.
        AssertionError
            If the on-disk data does not have `csc` or `csr` format.
        """
        # Prep variables
        shape = self.shape
        if isinstance(sparse_matrix, BaseCompressedSparseDataset):
            sparse_matrix = sparse_matrix._to_backed()

        # Check input
        if not ss.issparse(sparse_matrix):
            raise NotImplementedError(
                "Currently, only sparse matrices of equivalent format can be "
                "appended to a SparseDataset."
            )
        if self.format not in {"csr", "csc"}:
            raise NotImplementedError(
                f"The append method for format {self.format} " f"is not implemented."
            )
        if self.format != sparse_matrix.format:
            raise ValueError(
                f"Matrices must have same format. Currently are "
                f"{self.format!r} and {sparse_matrix.format!r}"
            )
        indptr_offset = len(self.group["indices"])
        if self.group["indptr"].dtype == np.int32:
            new_nnz = indptr_offset + len(sparse_matrix.indices)
            if new_nnz >= np.iinfo(np.int32).max:
                raise OverflowError(
                    "This array was written with a 32 bit intptr, but is now large "
                    "enough to require 64 bit values. Please recreate the array with "
                    "a 64 bit indptr."
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
            raise AssertionError("We forgot to update this branching to a new format")
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
        indptr.resize((orig_data_size + sparse_matrix.indptr.shape[0] - 1,))
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + indptr_offset
        )
        # Clear cached property
        if hasattr(self, "indptr"):
            del self._indptr

        # indices
        indices = self.group["indices"]
        orig_data_size = indices.shape[0]
        indices.resize((orig_data_size + sparse_matrix.indices.shape[0],))
        indices[orig_data_size:] = sparse_matrix.indices

    @cached_property
    def _indptr(self) -> np.ndarray:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        arr = self.group["indptr"][...]
        return arr

    def _to_backed(self) -> BackedSparseMatrix:
        format_class = get_backed_class(self.format)
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group["data"]
        mtx.indices = self.group["indices"]
        mtx.indptr = self._indptr
        return mtx

    def to_memory(self) -> ss.csr_matrix | ss.csc_matrix | SpArray:
        format_class = get_memory_class(
            self.format, use_sparray_in_io=settings.use_sparse_array_on_read
        )
        mtx = format_class(self.shape, dtype=self.dtype)
        mtx.data = self.group["data"][...]
        mtx.indices = self.group["indices"][...]
        mtx.indptr = self._indptr
        return mtx


class _CSRDataset(BaseCompressedSparseDataset, abc.CSRDataset):
    """Internal concrete version of :class:`anndata.abc.CSRDataset`."""


class _CSCDataset(BaseCompressedSparseDataset, abc.CSCDataset):
    """Internal concrete version of :class:`anndata.abc.CSRDataset`."""


def sparse_dataset(group: GroupStorageType) -> abc.CSRDataset | abc.CSCDataset:
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

    First we'll need a stored dataset:

    >>> import scanpy as sc
    >>> import h5py
    >>> from anndata.io import sparse_dataset
    >>> from anndata.io import read_elem
    >>> sc.datasets.pbmc68k_reduced().raw.to_adata().write_h5ad("pbmc.h5ad")

    Initialize a sparse dataset from storage

    >>> f = h5py.File("pbmc.h5ad")
    >>> X = sparse_dataset(f["X"])
    >>> X
    CSRDataset: backend hdf5, shape (700, 765), data_dtype float32

    Indexing returns sparse matrices

    >>> X[100:200]  # doctest: +ELLIPSIS
    <...sparse matrix of...float32...with 25003 stored elements...>

    These can also be used inside of an AnnData object, no need for backed mode

    >>> from anndata import AnnData
    >>> adata = AnnData(
    ...     layers={"backed": X}, obs=read_elem(f["obs"]), var=read_elem(f["var"])
    ... )
    >>> adata.layers["backed"]
    CSRDataset: backend hdf5, shape (700, 765), data_dtype float32

    Indexing access (i.e., from views) brings selection into memory

    >>> adata[adata.obs["bulk_labels"] == "CD56+ NK"].layers[
    ...     "backed"
    ... ]  # doctest: +ELLIPSIS
    <...sparse matrix of...float32...with 7340 stored elements...>
    """
    encoding_type = _get_group_format(group)
    if encoding_type == "csr":
        return _CSRDataset(group)
    elif encoding_type == "csc":
        return _CSCDataset(group)
    raise ValueError(f"Unknown encoding type {encoding_type}")


@_subset.register(BaseCompressedSparseDataset)
def subset_sparsedataset(d, subset_idx):
    return d[subset_idx]
