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

from abc import ABC
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cached_property, singledispatchmethod
from itertools import accumulate, chain, pairwise
from math import floor
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, cast

import h5py
import numpy as np
import scipy.sparse as ss
import zarr

from testing.anndata._doctest import doctest_filterwarnings

from .. import abc
from .._settings import settings
from ..compat import (
    CSArray,
    CSMatrix,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    _read_attr,
)
from .index import (
    _ensure_numpy_idx,
    _fix_slice_bounds,
    _subset_dispatch,
    unpack_index,
)

if TYPE_CHECKING:
    from collections.abc import Mapping
    from types import EllipsisType, ModuleType
    from typing import Any, Literal

    from .._types import _ArrayStorageType, _GroupStorageType
    from ..typing import Index, Index1D
    from .index import NumpySubsetIdx


type DenseType = np.ndarray | CupyArray
type SparseMatrixType = CupyCSCMatrix | CupyCSRMatrix | CSMatrix | CSArray


def is_gpu() -> bool:
    config: Mapping[str, str] = getattr(zarr, "config", {"ndbuffer": "cpu"})
    return "gpu" in config.get("ndbuffer", "cpu")


def get_np_module() -> ModuleType:
    if is_gpu():
        import cupy as cp

        return cp
    return np


class CompressedVectors[DenseT: DenseType](NamedTuple):
    data: DenseT
    indices: DenseT
    indptr: DenseT

    @classmethod
    def from_buffers(
        cls,
        data: DenseT,
        indices: DenseT,
        indptr: DenseType,
    ) -> CompressedVectors[DenseT]:
        return cls(data, indices, get_np_module().array(indptr))


def _read_dense(
    arr: DenseType | _ArrayStorageType, idx: slice | np.ndarray | EllipsisType
) -> DenseType:
    return arr[idx]


def _index_in_memory(
    mtx: SparseMatrixType, row: Index1D, col: Index1D
) -> SparseMatrixType | float:
    # scipy’s stubs describe fewer index types than its sparse classes accept
    return cast("Any", mtx)[row, col]


def _group_array(group: _GroupStorageType, key: str) -> _ArrayStorageType:
    if not isinstance(arr := group[key], zarr.Array | h5py.Dataset):
        msg = f"Expected {key!r} of {group} to be an array, got {type(arr)}"
        raise ValueError(msg)
    return arr


def slice_len(s: slice, l: int) -> int:
    """Returns length of `a[s]` where `len(a) == l`."""
    return len(range(*s.indices(l)))


def slice_as_int(s: slice, l: int) -> int:
    """Converts slices of length 1 to the integer index they’ll access."""
    out = list(range(*s.indices(l)))
    assert len(out) == 1
    return out[0]


@dataclass
class BackedSparseMatrix[ArrayT: _ArrayStorageType]:
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`,
    since that calls copy on `.data`, `.indices`, and `.indptr`.
    """

    data: ArrayT
    indices: ArrayT
    indptr: DenseType | ArrayT
    format: Literal["csr", "csc"]
    shape: tuple[int, int]

    @property
    def memory_format(self) -> type[SparseMatrixType]:
        if is_gpu():
            return CupyCSRMatrix if self.format == "csr" else CupyCSCMatrix
        if settings.use_sparse_array_on_read:
            return ss.csr_array if self.format == "csr" else ss.csc_array
        return ss.csr_matrix if self.format == "csr" else ss.csc_matrix

    @property
    def np_module(self) -> ModuleType:
        return get_np_module()

    @property
    def major_axis(self) -> Literal[0, 1]:
        return 1 if self.format == "csc" else 0

    @property
    def minor_axis(self) -> Literal[0, 1]:
        return 1 if self.format == "csr" else 0

    @property
    def major_axis_size(self) -> int:
        return self.shape[self.major_axis]

    @property
    def minor_axis_size(self) -> int:
        return self.shape[self.minor_axis]

    def copy(self) -> SparseMatrixType:
        if isinstance(self.data, h5py.Dataset):
            return sparse_dataset(self.data.parent).to_memory()
        if isinstance(self.data, zarr.Array):
            anndata_group = zarr.open_group(store=self.data.store, mode="r")
            sparse_group = anndata_group[
                str(
                    Path(str(self.data.store_path))
                    .relative_to(str(anndata_group.store_path))
                    .parent
                )
            ]
            if not isinstance(sparse_group, zarr.Group):
                msg = f"Expected a group backing {self.data.store_path}, got {type(sparse_group)}"
                raise ValueError(msg)
            return sparse_dataset(sparse_group).to_memory()
        msg = f"Unsupported array types {type(self.data)}"
        raise ValueError(msg)

    def _get_contiguous_compressed_slice(
        self, s: slice[int, int, Literal[1] | None]
    ) -> CompressedVectors:
        new_indptr: DenseType = _read_dense(
            self.indptr, slice(s.start, s.stop + 1)
        ).copy()

        start = new_indptr[0].item()
        stop = new_indptr[-1].item()

        new_indptr -= start

        new_data = _read_dense(self.data, slice(start, stop))
        new_indices = _read_dense(self.indices, slice(start, stop))

        return CompressedVectors.from_buffers(new_data, new_indices, new_indptr)

    def get_compressed_vectors(self, row_idxs: Iterable[int]) -> CompressedVectors:
        indptr_slices = [
            slice(*_read_dense(self.indptr, slice(i, i + 2))) for i in row_idxs
        ]
        data: DenseType
        indices: DenseType
        # HDF5 cannot handle out-of-order integer indexing
        if isinstance(self.data, zarr.Array):
            as_np_indptr = self.np_module.concatenate([
                self.np_module.arange(s.start, s.stop) for s in indptr_slices
            ])
            data = _read_dense(self.data, as_np_indptr)
            indices = _read_dense(self.indices, as_np_indptr)
        else:
            data = np.concatenate([_read_dense(self.data, s) for s in indptr_slices])
            indices = np.concatenate([
                _read_dense(self.indices, s) for s in indptr_slices
            ])
        indptr = self.np_module.array(
            list(accumulate(chain((0,), (s.stop - s.start for s in indptr_slices))))
        )
        return CompressedVectors.from_buffers(data, indices, indptr)

    def get_compressed_vectors_for_slices(
        self, slices: Sequence[slice]
    ) -> CompressedVectors:
        indptr_indices = [
            _read_dense(self.indptr, slice(s.start, s.stop + 1)) for s in slices
        ]
        indptr_limits = [slice(i[0].item(), i[-1].item()) for i in indptr_indices]
        data: DenseType
        indices: DenseType
        # HDF5 cannot handle out-of-order integer indexing
        if isinstance(self.data, zarr.Array):
            indptr_int = self.np_module.concatenate([
                self.np_module.arange(s.start, s.stop) for s in indptr_limits
            ])
            data = _read_dense(self.data, indptr_int)
            indices = _read_dense(self.indices, indptr_int)
        else:
            data = np.concatenate([_read_dense(self.data, s) for s in indptr_limits])
            indices = np.concatenate([
                _read_dense(self.indices, s) for s in indptr_limits
            ])
        # Need to track the size of the gaps in the slices to each indptr subselection
        gaps = (s1.start - s0.stop for s0, s1 in pairwise(indptr_limits))
        offsets = accumulate(chain([indptr_limits[0].start], gaps))
        start_indptr = indptr_indices[0] - next(offsets)
        if len(slices) < 2:  # there is only one slice so no need to concatenate
            return CompressedVectors.from_buffers(data, indices, start_indptr)
        end_indptr = self.np_module.concatenate([
            s[1:] - o for s, o in zip(indptr_indices[1:], offsets, strict=False)
        ])
        indptr = self.np_module.concatenate([start_indptr, end_indptr])
        return CompressedVectors.from_buffers(data, indices, indptr)

    def get_compressed_vector(self, idx: int) -> CompressedVectors:
        s = slice(*_read_dense(self.indptr, slice(idx, idx + 2)))
        data = _read_dense(self.data, s)
        indices = _read_dense(self.indices, s)
        indptr: DenseType = [0, len(data)]
        return CompressedVectors.from_buffers(data, indices, indptr)

    def __getitem__(self, key) -> SparseMatrixType:
        row, col = key if isinstance(key, tuple) else (key, slice(None))
        major_index, minor_index = (row, col) if self.format == "csr" else (col, row)
        return self._get(major_index, minor_index)

    def _gen_maj_min_tuple[Ma, Mi](
        self, major: Ma, minor: Mi
    ) -> tuple[Ma, Mi] | tuple[Mi, Ma]:
        return (major, minor) if self.format == "csr" else (minor, major)

    @singledispatchmethod
    def _get(
        self, major_index: int | slice | Sequence | np.ndarray, minor_index: slice
    ) -> SparseMatrixType:
        vectors: CompressedVectors = CompressedVectors(
            _read_dense(self.data, ...),
            _read_dense(self.indices, ...),
            _read_dense(self.indptr, ...),
        )
        key = self._gen_maj_min_tuple(major_index, minor_index)
        # scipy-stubs’ __getitem__ overloads only cover fixed index shapes,
        # not this dispatch fallback’s generic major_index/minor_index union.
        return self.memory_format(vectors, shape=self.shape)[key]  # type: ignore[index]

    @_get.register
    def _get_intXslice(self, major_index: int, minor_index: slice) -> SparseMatrixType:
        return self.memory_format(
            self.get_compressed_vector(major_index),
            shape=self._gen_maj_min_tuple(1, self.minor_axis_size),
        )[self._gen_maj_min_tuple(slice(None), minor_index)]

    @_get.register
    def _get_sliceXslice(
        self, major_index: slice, minor_index: slice
    ) -> SparseMatrixType:
        major_index = _fix_slice_bounds(major_index, self.shape[self.major_axis])
        minor_index = _fix_slice_bounds(minor_index, self.shape[self.minor_axis])

        major_index_size = slice_len(major_index, self.shape[self.major_axis])
        if major_index_size == 1:
            return self._get(
                slice_as_int(major_index, self.shape[self.major_axis]), minor_index
            )
        if major_index.step != 1:
            return self._get(
                np.arange(*major_index.indices(self.shape[self.major_axis])),
                minor_index,
            )
        compressed_vectors = self._get_contiguous_compressed_slice(
            slice(major_index.start, major_index.stop)
        )
        return self.memory_format(
            compressed_vectors,
            shape=self._gen_maj_min_tuple(major_index_size, self.minor_axis_size),
        )[self._gen_maj_min_tuple(slice(None), minor_index)]

    @_get.register
    def _get_arrayXslice(
        self, major_index: Sequence | np.ndarray, minor_index: slice
    ) -> SparseMatrixType:
        if len(major_index) == 0:
            return self.memory_format(
                (0, self.minor_axis_size)
                if self.format == "csr"
                else (self.minor_axis_size, 0)
            )
        if isinstance(major_index, np.ndarray) and major_index.dtype == bool:
            major_index = np.flatnonzero(major_index)
        out_shape = self._gen_maj_min_tuple(len(major_index), self.minor_axis_size)
        return self.memory_format(
            self.get_compressed_vectors(major_index), shape=out_shape
        )[self._gen_maj_min_tuple(slice(None), minor_index)]

    def subset_by_major_axis_mask(
        self: BackedSparseMatrix, mask: np.ndarray
    ) -> CompressedVectors:
        slices: list[slice] = np.ma.clump_unmasked(np.ma.masked_array(mask, mask=~mask))

        def mean_slice_length(slices) -> int:
            return floor(sum(s.stop - s.start for s in slices) / len(slices))

        # heuristic for whether slicing should be optimized
        if len(slices) > 0:
            if mean_slice_length(slices) <= 7:
                return self.get_compressed_vectors(np.where(mask)[0])
            else:
                return self.get_compressed_vectors_for_slices(slices)
        return CompressedVectors.from_buffers(
            self.np_module.array([]),
            self.np_module.array([]),
            self.np_module.array([0]),
        )


def _get_group_format(group: _GroupStorageType) -> str:
    if "h5sparse_format" in group.attrs:
        # TODO: Warn about an old format
        # If this is only just going to be public, I could insist it's not like this
        fmt = _read_attr(group.attrs, "h5sparse_format")
    else:
        # Should this be an extra field?
        fmt = _read_attr(group.attrs, "encoding-type")
    if not isinstance(fmt, str):
        msg = f"Expected a string sparse format for {group}, got {fmt!r}"
        raise ValueError(msg)
    return fmt.removesuffix("_matrix")


# Check for the overridden few methods above in our BackedSparseMatrix subclasses
def is_sparse_indexing_overridden(
    format: Literal["csr", "csc"], row: Index1D, col: Index1D
) -> bool:
    major_indexer, minor_indexer = (row, col) if format == "csr" else (col, row)
    return isinstance(minor_indexer, slice) and (
        isinstance(major_indexer, int | np.integer | slice)
        or (isinstance(major_indexer, np.ndarray) and major_indexer.ndim == 1)
    )


class BaseCompressedSparseDataset[GroupT: _GroupStorageType](
    abc._AbstractCSDataset, ABC
):
    _group: GroupT
    _should_cache_indptr: bool

    def __init__(self, group: GroupT, *, should_cache_indptr: bool = True):
        type(self)._check_group_format(group)
        self._group = group
        self._should_cache_indptr = should_cache_indptr

    @property
    def group(self) -> GroupT:
        """The group underlying the backed matrix."""
        return self._group

    @group.setter
    def group(self, val):
        msg = f"Do not reset group on a {type(self)} with {val}.  Instead use `sparse_dataset` to make a new class."
        raise AttributeError(msg)

    @property
    def backend(self) -> Literal["zarr", "hdf5"]:
        """Which file type is used on-disk."""
        if isinstance(self.group, zarr.Group):
            return "zarr"
        elif isinstance(self.group, h5py.Group):
            return "hdf5"
        else:
            msg = f"Unknown group type {type(self.group)}"
            raise ValueError(msg)

    @property
    def dtype(self) -> np.dtype:
        """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""
        return self._data.dtype

    @classmethod
    def _check_group_format(cls, group: GroupT) -> None:
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
            shape = _read_attr(self.group.attrs, "h5sparse_shape")
        n_obs, n_vars = cast("Iterable[int]", shape)
        return int(n_obs), int(n_vars)

    def __repr__(self) -> str:
        name = type(self).__name__.removeprefix("_")
        return f"{name}: backend {self.backend}, shape {self.shape}, data_dtype {self.dtype}"

    def __getitem__(self, index: Index | tuple[()]) -> float | SparseMatrixType:
        row, col = self._normalize_index(index)
        mtx = self._to_backed()

        sub: SparseMatrixType | float
        # Handle masked indexing along major axis
        major_index = np.asarray(row if self.format == "csr" else col)
        if major_index.dtype == bool:
            vectors = mtx.subset_by_major_axis_mask(major_index)
            major_index_size = int(major_index.sum())
            if self.format == "csr":
                shape = (major_index_size, mtx.shape[1])
                sub = _index_in_memory(
                    ss.csr_matrix(vectors, shape=shape), slice(None), col
                )
            else:
                shape = (mtx.shape[0], major_index_size)
                sub = _index_in_memory(
                    ss.csc_matrix(vectors, shape=shape), row, slice(None)
                )
        # read into memory data if we do not override access methods
        elif not is_sparse_indexing_overridden(self.format, row, col):
            sub = _index_in_memory(self.to_memory(), row, col)
        else:
            sub = mtx[row, col]

        # Scalars are returned as-is, everything else in the configured memory format
        mtx_fmt = mtx.memory_format
        if issubclass(mtx_fmt, CSArray) and isinstance(sub, CSMatrix):
            return mtx_fmt(sub)
        return sub

    def _normalize_index(self, index: Index | tuple[()]) -> tuple[Index1D, Index1D]:
        if isinstance(index, tuple) and not len(index):
            index = slice(None)
        row, col = unpack_index(index)
        if isinstance(row, Iterable) and isinstance(col, Iterable):
            row, col = np.ix_(np.asarray(row), np.asarray(col))
        return row, col

    # TODO: split to other classes?
    def append(  # noqa: PLR0912, PLR0915
        self,
        sparse_matrix: CSMatrix | CSArray | abc._AbstractCSDataset | BackedSparseMatrix,
    ) -> None:
        # Prep variables
        shape = self.shape
        if isinstance(sparse_matrix, BaseCompressedSparseDataset):
            sparse_matrix = sparse_matrix._to_backed()

        # Check input
        if not (
            ss.issparse(sparse_matrix) or isinstance(sparse_matrix, BackedSparseMatrix)
        ):
            msg = (
                "Currently, only sparse matrices of equivalent format can be "
                f"appended to a SparseDataset, not {sparse_matrix}."
            )
            raise NotImplementedError(msg)
        if self.format not in {"csr", "csc"}:
            msg = f"The append method for format {self.format} is not implemented."
            raise NotImplementedError(msg)
        if self.format != sparse_matrix.format:
            msg = (
                f"Matrices must have same format. Currently are "
                f"{self.format!r} and {sparse_matrix.format!r}"
            )
            raise ValueError(msg)
        [indptr_offset] = _group_array(self.group, "indices").shape
        if _group_array(self.group, "indptr").dtype == np.int32:
            new_nnz = indptr_offset + sparse_matrix.indices.shape[0]
            if new_nnz >= np.iinfo(np.int32).max:
                msg = (
                    "This array was written with a 32 bit intptr, but is now large "
                    "enough to require 64 bit values. Please recreate the array with "
                    "a 64 bit indptr."
                )
                raise OverflowError(msg)

        # shape
        if self.format == "csr":
            assert shape[1] == sparse_matrix.shape[1], (
                "CSR matrices must have same size of dimension 1 to be appended."
            )
            new_shape = (shape[0] + sparse_matrix.shape[0], shape[1])
        elif self.format == "csc":
            assert shape[0] == sparse_matrix.shape[0], (
                "CSC matrices must have same size of dimension 0 to be appended."
            )
            new_shape = (shape[0], shape[1] + sparse_matrix.shape[1])
        else:
            msg = "We forgot to update this branching to a new format"
            raise AssertionError(msg)
        if "h5sparse_shape" in self.group.attrs:
            del self.group.attrs["h5sparse_shape"]
        self.group.attrs["shape"] = new_shape

        # data
        data = _group_array(self.group, "data")
        orig_data_size = data.shape[0]
        data.resize((orig_data_size + sparse_matrix.data.shape[0],))
        # see https://github.com/zarr-developers/zarr-python/discussions/2712 for why we need to read first
        append_data = sparse_matrix.data
        append_indices = sparse_matrix.indices
        if isinstance(sparse_matrix.data, zarr.Array):
            data[orig_data_size:] = append_data[...]
        else:
            data[orig_data_size:] = append_data
        # indptr
        indptr = _group_array(self.group, "indptr")
        orig_data_size = indptr.shape[0]
        indptr.resize((orig_data_size + sparse_matrix.indptr.shape[0] - 1,))
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + indptr_offset
        )

        # indices
        if isinstance(sparse_matrix.data, zarr.Array):
            append_indices = append_indices[...]
        indices = _group_array(self.group, "indices")
        orig_data_size = indices.shape[0]
        indices.resize((orig_data_size + sparse_matrix.indices.shape[0],))
        indices[orig_data_size:] = append_indices

        # Clear cached property
        for attr in ["_indptr", "_indices", "_data"]:
            if hasattr(self, attr):
                delattr(self, attr)

    @cached_property
    def _indptr(self) -> DenseType | _ArrayStorageType:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        indptr = _group_array(self.group, "indptr")
        if self._should_cache_indptr:
            return _read_dense(indptr, ...)
        return indptr

    @cached_property
    def _indices(self) -> _ArrayStorageType:
        """\
        Cache access to the indices to prevent unnecessary reads of the zarray
        """
        return _group_array(self.group, "indices")

    @cached_property
    def _data(self) -> _ArrayStorageType:
        """\
        Cache access to the data to prevent unnecessary reads of the zarray
        """
        return _group_array(self.group, "data")

    def _to_backed(self) -> BackedSparseMatrix:
        mtx = BackedSparseMatrix(
            format=self.format,
            data=self._data,
            indices=self._indices,
            indptr=self._indptr,
            shape=self.shape,
        )
        return mtx

    def to_memory(self) -> SparseMatrixType:
        backed_class = BackedSparseMatrix(
            format=self.format,
            data=self._data,
            indices=self._indices,
            indptr=self._indptr,
            shape=self.shape,
        )
        mtx = backed_class.memory_format(self.shape, dtype=self.dtype)
        mtx.data = _read_dense(self._data, ...)
        mtx.indices = _read_dense(self._indices, ...)
        mtx.indptr = _read_dense(self._indptr, ...)
        return mtx


class _CSRDataset(BaseCompressedSparseDataset, abc.CSRDataset):
    """Internal concrete version of :class:`anndata.abc.CSRDataset`."""


class _CSCDataset(BaseCompressedSparseDataset, abc.CSCDataset):
    """Internal concrete version of :class:`anndata.abc.CSRDataset`."""


@doctest_filterwarnings("ignore", r"Moving element.*uns.*to.*obsp", FutureWarning)
def sparse_dataset(
    group: _GroupStorageType,
    *,
    should_cache_indptr: bool = True,
) -> abc.CSRDataset | abc.CSCDataset:
    """Generates a backed mode-compatible sparse dataset class.

    Parameters
    ----------
    group
        The backing group store.
    should_cache_indptr
        Whether or not to cache the indptr for repeated reuse as a :class:`numpy.ndarray`.
        The default is `True` but one might set it to false if the dataset is repeatedly reopened
        using this command, and then only a subset is read in before closing again.
        See https://github.com/scverse/anndata/blob/3c489b979086c39c59d3eb5dad90ebacce3b9a80/src/anndata/_io/specs/lazy_methods.py#L85-L95
        for the target use-case.

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
        return _CSRDataset(group, should_cache_indptr=should_cache_indptr)
    elif encoding_type == "csc":
        return _CSCDataset(group, should_cache_indptr=should_cache_indptr)
    msg = f"Unknown encoding type {encoding_type}"
    raise ValueError(msg)


# registered on the abstract base so the public `abc.CS[RC]Dataset` types dispatch too
@_subset_dispatch.register(abc._AbstractCSDataset)
@_ensure_numpy_idx
def subset_sparsedataset(
    d: abc._AbstractCSDataset, subset_idx: NumpySubsetIdx
) -> float | CSMatrix | CSArray:
    return d[subset_idx]
