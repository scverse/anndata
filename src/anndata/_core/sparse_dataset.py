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

import asyncio
from abc import ABC
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cached_property, singledispatchmethod
from itertools import accumulate, chain, pairwise
from math import floor
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import h5py
import numpy as np
import scipy
import scipy.sparse as ss
from packaging.version import Version

from .. import abc
from .._settings import settings
from ..compat import CSArray, H5Group, ZarrArray, ZarrGroup, _read_attr, is_zarr_v2
from .index import _fix_slice_bounds, _subset, unpack_index

if TYPE_CHECKING:
    from typing import Any, Literal

    from anndata._types import ArrayStorageType

    from .._types import GroupStorageType
    from ..compat import CSMatrix, H5Array
    from .index import Index, Index1D

SCIPY_1_15 = Version(scipy.__version__) >= Version("1.15rc0")


async def index_array(elem: ArrayStorageType | np.ndarray, selection) -> np.ndarray:
    from ..compat import ZarrAsyncArray

    if isinstance(elem, ZarrAsyncArray):
        return await elem.getitem(selection)
    return elem[selection]


async def get(elem: GroupStorageType, key) -> ArrayStorageType | GroupStorageType:
    from ..compat import ZarrAsyncGroup

    if isinstance(elem, ZarrAsyncGroup):
        return await elem.get(key)
    return elem[key]


class CompressedVectors(NamedTuple):
    data: np.ndarray
    indices: np.ndarray
    indptr: np.ndarray


def slice_len(s: slice, l: int) -> int:
    """Returns length of `a[s]` where `len(a) == l`."""
    return len(range(*s.indices(l)))


def slice_as_int(s: slice, l: int) -> int:
    """Converts slices of length 1 to the integer index they’ll access."""
    out = list(range(*s.indices(l)))
    assert len(out) == 1
    return out[0]


async def select_many(arr: ArrayStorageType, slices) -> np.ndarray:
    return np.concatenate(await asyncio.gather(*[index_array(arr, s) for s in slices]))


@dataclass
class BackedSparseMatrix:
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`,
    since that calls copy on `.data`, `.indices`, and `.indptr`.
    """

    data: ArrayStorageType
    indices: ArrayStorageType
    indptr: np.ndarray
    format: Literal["csr", "csc"]
    shape: tuple[int, int]

    @cached_property
    def memory_format(self):
        if self.format == "csr":
            return ss.csr_array if settings.use_sparse_array_on_read else ss.csr_matrix
        return ss.csc_array if settings.use_sparse_array_on_read else ss.csc_matrix

    @property
    def major_axis(self):
        return self.format == "csc"

    @property
    def minor_axis(self):
        return self.format == "csr"

    @property
    def major_axis_size(self):
        return self.shape[self.major_axis]

    @property
    def minor_axis_size(self):
        return self.shape[self.minor_axis]

    def copy(self) -> CSMatrix:
        if isinstance(self.data, h5py.Dataset):
            return sparse_dataset(self.data.parent).to_memory()
        if isinstance(self.data, ZarrArray):
            import zarr

            if is_zarr_v2():
                sparse_group = zarr.open(
                    store=self.data.store,
                    mode="r",
                    chunk_store=self.data.chunk_store,  # chunk_store is needed, not clear why
                )[Path(self.data.path).parent]
            else:
                anndata_group = zarr.open_group(store=self.data.store, mode="r")
                sparse_group = anndata_group[
                    str(
                        Path(str(self.data.store_path))
                        .relative_to(str(anndata_group.store_path))
                        .parent
                    )
                ]
            return sparse_dataset(sparse_group).to_memory()
        msg = f"Unsupported array types {type(self.data)}"
        raise ValueError(msg)

    async def _get_contiguous_compressed_slice(self, s: slice) -> CompressedVectors:
        new_indptr: np.ndarray = self.indptr[s.start : s.stop + 1].copy()

        start = new_indptr[0]
        stop = new_indptr[-1]

        new_data, new_indices = await asyncio.gather(
            index_array(self.data, slice(start, stop)),
            index_array(self.indices, slice(start, stop)),
        )

        new_indptr -= start

        return CompressedVectors(new_data, new_indices, new_indptr)

    async def get_compressed_vectors(
        self, row_idxs: Iterable[int]
    ) -> CompressedVectors:
        indptr_slices = [slice(*(self.indptr[i : i + 2])) for i in row_idxs]
        # HDF5 cannot handle out-of-order integer indexing
        if isinstance(self.data, ZarrArray):
            as_np_indptr = np.concatenate(
                [np.arange(s.start, s.stop) for s in indptr_slices]
            )
            data, indices = await asyncio.gather(
                index_array(self.data, as_np_indptr),
                index_array(self.indices, as_np_indptr),
            )
        else:
            data, indices = await asyncio.gather(
                select_many(self.data, indptr_slices),
                select_many(self.indices, indptr_slices),
            )
        indptr = np.array(
            list(accumulate(chain((0,), (s.stop - s.start for s in indptr_slices))))
        )
        return CompressedVectors(data, indices, indptr)

    async def get_compressed_vectors_for_slices(
        self, slices: Iterable[slice]
    ) -> CompressedVectors:
        indptr_indices = [self.indptr[slice(s.start, s.stop + 1)] for s in slices]
        indptr_limits = [slice(i[0], i[-1]) for i in indptr_indices]
        # HDF5 cannot handle out-of-order integer indexing
        if isinstance(self.data, ZarrArray):
            indptr_int = np.concatenate(
                [np.arange(s.start, s.stop) for s in indptr_limits]
            )
            data, indices = await asyncio.gather(
                index_array(self.data, indptr_int),
                index_array(self.indices, indptr_int),
            )
        else:
            data, indices = await asyncio.gather(
                select_many(self.data, indptr_limits),
                select_many(self.indices, indptr_limits),
            )
        # Need to track the size of the gaps in the slices to each indptr subselection
        gaps = (s1.start - s0.stop for s0, s1 in pairwise(indptr_limits))
        offsets = accumulate(chain([indptr_limits[0].start], gaps))
        start_indptr = indptr_indices[0] - next(offsets)
        if len(slices) < 2:  # there is only one slice so no need to concatenate
            return CompressedVectors(data, indices, start_indptr)
        end_indptr = np.concatenate(
            [s[1:] - o for s, o in zip(indptr_indices[1:], offsets)]
        )
        indptr = np.concatenate([start_indptr, end_indptr])
        return CompressedVectors(data, indices, indptr)

    async def get_compressed_vector(self, idx: int) -> CompressedVectors:
        s = slice(*(self.indptr[idx : idx + 2]))
        data, indices = await asyncio.gather(
            index_array(self.data, s), index_array(self.indices, s)
        )
        indptr: np.ndarray = [0, len(data)]
        return CompressedVectors(data, indices, indptr)

    async def getitem(self, key):
        if isinstance(key, tuple):
            row, col = key
        else:
            row = key
            col = slice(None)
        major_index, minor_index = (row, col) if self.format == "csr" else (col, row)
        return await self._get(major_index, minor_index)

    def _gen_index(self, major_index: Any, minor_index: Any):
        return (
            (major_index, minor_index)
            if self.format == "csr"
            else (minor_index, major_index)
        )

    @singledispatchmethod
    async def _get(self, major_index: Any, minor_index: slice) -> CSMatrix | CSArray:
        indices, data = await asyncio.gather(
            index_array(self.indices, Ellipsis), index_array(self.data, Ellipsis)
        )
        return self.memory_format((data, indices, self.indptr))[
            self._gen_index(major_index, minor_index)
        ]

    @_get.register
    async def _get_intXslice(
        self, major_index: int, minor_index: slice
    ) -> ss.csr_matrix:
        return self.memory_format(
            await self.get_compressed_vector(major_index),
            shape=(1, self.minor_axis_size)
            if self.format == "csr"
            else (self.minor_axis_size, 1),
        )[self._gen_index(slice(None), minor_index)]

    @_get.register
    async def _get_sliceXslice(
        self, major_index: slice, minor_index: slice
    ) -> ss.csr_matrix:
        major_index = _fix_slice_bounds(major_index, self.shape[self.major_axis])
        minor_index = _fix_slice_bounds(minor_index, self.shape[self.minor_axis])

        major_index_size = slice_len(major_index, self.shape[self.major_axis])
        minor_index_size = slice_len(minor_index, self.shape[self.minor_axis])

        out_shape = (
            (major_index_size, minor_index_size)
            if self.format == "csr"
            else (minor_index_size, major_index_size)
        )
        if out_shape[self.major_axis] == 1:
            return await self._get_intXslice(
                slice_as_int(major_index, self.shape[self.major_axis]), minor_index
            )
        if major_index.step != 1:
            return await self._get_arrayXslice(
                np.arange(*major_index.indices(self.shape[self.major_axis])),
                minor_index,
            )
        compressed_vectors = await self._get_contiguous_compressed_slice(major_index)
        return self.memory_format(
            compressed_vectors,
            shape=(out_shape[self.major_axis], self.minor_axis_size)
            if self.format == "csr"
            else (self.minor_axis_size, out_shape[self.major_axis]),
        )[self._gen_index(slice(None), minor_index)]

    @_get.register
    async def _get_arrayXslice(
        self, major_index: Sequence | np.ndarray, minor_index: slice
    ) -> ss.csr_matrix:
        idxs = np.asarray(major_index)
        if len(idxs) == 0:
            return self.memory_format(
                (0, self.minor_axis_size)
                if self.format == "csr"
                else (self.minor_axis_size, 0)
            )
        if idxs.dtype == bool:
            idxs = np.where(idxs)
        out_shape = (
            (len(idxs), self.minor_axis_size)
            if self.format == "csr"
            else (self.minor_axis_size, len(idxs))
        )
        return self.memory_format(
            await self.get_compressed_vectors(idxs), shape=out_shape
        )[self._gen_index(slice(None), minor_index)]

    async def subset_by_major_axis_mask(
        self: BackedSparseMatrix, mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        slices = np.ma.extras._ezclump(mask)

        def mean_slice_length(slices):
            return floor(sum(s.stop - s.start for s in slices) / len(slices))

        # heuristic for whether slicing should be optimized
        if len(slices) > 0:
            if mean_slice_length(slices) <= 7:
                return await self.get_compressed_vectors(np.where(mask)[0])
            else:
                return await self.get_compressed_vectors_for_slices(slices)
        return [], [], [0]


def _get_group_format(group: GroupStorageType) -> str:
    if "h5sparse_format" in group.attrs:
        # TODO: Warn about an old format
        # If this is only just going to be public, I could insist it's not like this
        return _read_attr(group.attrs, "h5sparse_format")
    else:
        # Should this be an extra field?
        return _read_attr(group.attrs, "encoding-type").replace("_matrix", "")


# Check for the overridden few methods above in our BackedSparseMatrix subclasses
def is_sparse_indexing_overridden(
    format: Literal["csr", "csc"], row: Index1D, col: Index1D
):
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
        msg = f"Do not reset group on a {type(self)} with {val}.  Instead use `sparse_dataset` to make a new class."
        raise AttributeError(msg)

    @property
    def backend(self) -> Literal["zarr", "hdf5"]:
        """Which file type is used on-disk."""
        if isinstance(self.group, ZarrGroup):
            return "zarr"
        elif isinstance(self.group, H5Group):
            return "hdf5"
        else:
            msg = f"Unknown group type {type(self.group)}"
            raise ValueError(msg)

    @property
    def dtype(self) -> np.dtype:
        """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""
        return self._data.dtype

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

    async def getitem(self, index: Index) -> float | CSMatrix | CSArray:
        indices = self._normalize_index(index)
        row, col = indices
        mtx = await self._to_backed_async()

        # Handle masked indexing along major axis
        if self.format == "csr" and np.array(row).dtype == bool:
            sub = ss.csr_matrix(
                await mtx.subset_by_major_axis_mask(row),
                shape=(row.sum(), mtx.shape[1]),
            )[:, col]
        elif self.format == "csc" and np.array(col).dtype == bool:
            sub = ss.csc_matrix(
                await mtx.subset_by_major_axis_mask(col),
                shape=(mtx.shape[0], col.sum()),
            )[row, :]
        # read into memory data if we do not override access methods
        elif not is_sparse_indexing_overridden(self.format, row, col):
            sub = (await mtx.getitem((slice(None), slice(None))))[row, col]
        else:
            sub = await mtx.getitem((row, col))

        # If indexing is array x array it returns a backed_sparse_matrix
        # Not sure what the performance is on that operation
        # Also need to check if memory format is not matrix
        mtx_fmt = mtx.memory_format
        must_convert_to_array = issubclass(mtx_fmt, CSArray) and not isinstance(
            sub, CSArray
        )
        if isinstance(sub, BackedSparseMatrix) or must_convert_to_array:
            return mtx_fmt(sub)
        else:
            return sub

    def __getitem__(self, index: Index) -> float | CSArray:
        return asyncio.run(self.getitem(index))

    def _normalize_index(
        self, index: Index | tuple[()]
    ) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(index, tuple) and not len(index):
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        return row, col

    # TODO: split to other classes?
    def append(self, sparse_matrix: CSMatrix | CSArray) -> None:
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
        [indptr_offset] = self.group["indices"].shape
        if self.group["indptr"].dtype == np.int32:
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
        data = self.group["data"]
        orig_data_size = data.shape[0]
        data.resize((orig_data_size + sparse_matrix.data.shape[0],))
        # see https://github.com/zarr-developers/zarr-python/discussions/2712 for why we need to read first
        append_data = sparse_matrix.data
        append_indices = sparse_matrix.indices
        if isinstance(sparse_matrix.data, ZarrArray) and not is_zarr_v2():
            data[orig_data_size:] = append_data[...]
        else:
            data[orig_data_size:] = append_data
        # indptr
        indptr = self.group["indptr"]
        orig_data_size = indptr.shape[0]
        indptr.resize((orig_data_size + sparse_matrix.indptr.shape[0] - 1,))
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + indptr_offset
        )

        # indices
        if isinstance(sparse_matrix.data, ZarrArray) and not is_zarr_v2():
            append_indices = append_indices[...]
        indices = self.group["indices"]
        orig_data_size = indices.shape[0]
        indices.resize((orig_data_size + sparse_matrix.indices.shape[0],))
        indices[orig_data_size:] = append_indices

        # Clear cached property
        for attr in ["_indptr", "_indices", "_data"]:
            if hasattr(self, attr):
                delattr(self, attr)

    @cached_property
    def _indptr(self) -> np.ndarray:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        arr = self.group["indptr"][...]
        return arr

    async def _indptr_async(self) -> np.ndarray:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        if not hasattr(self, "_indptr_cache_async") or self._indptr_cache_async is None:
            self._indptr_cache_async = await index_array(
                await get(self.group, "indptr"), Ellipsis
            )
        return self._indptr_cache_async

    async def _indices_async(self) -> ArrayStorageType:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        if (
            not hasattr(self, "_indices_cache_async")
            or self._indices_cache_async is None
        ):
            self._indices_cache_async = await get(self.group, "indices")
        return self._indices_cache_async

    async def _data_async(self) -> ArrayStorageType:
        """\
        Other than `data` and `indices`, this is only as long as the major axis

        It should therefore fit into memory, so we cache it for faster access.
        """
        if not hasattr(self, "_data_cache_async") or self._data_cache_async is None:
            self._data_cache_async = await get(self.group, "data")
        return self._data_cache_async

    @cached_property
    def _indices(self) -> H5Array | ZarrArray:
        """\
        Cache access to the indices to prevent unnecessary reads of the zarray
        """
        return self.group["indices"]

    @cached_property
    def _data(self) -> H5Array | ZarrArray:
        """\
        Cache access to the data to prevent unnecessary reads of the zarray
        """
        return self.group["data"]

    def _to_backed(self) -> BackedSparseMatrix:
        mtx = BackedSparseMatrix(
            format=self.format,
            data=self._data,
            indices=self._indices,
            indptr=self._indptr,
            shape=self.shape,
        )
        return mtx

    async def _to_backed_async(self) -> BackedSparseMatrix:
        indptr, indices, data = await asyncio.gather(
            self._indptr_async(), self._indices_async(), self._data_async()
        )
        mtx = BackedSparseMatrix(
            format=self.format,
            data=data,
            indices=indices,
            indptr=indptr,
            shape=self.shape,
        )
        return mtx

    def to_memory(self) -> CSMatrix | CSArray:
        backed_class = BackedSparseMatrix(
            format=self.format,
            data=self._data,
            indices=self._indices,
            indptr=self._indptr,
            shape=self.shape,
        )
        mtx = backed_class.memory_format(self.shape, dtype=self.dtype)
        mtx.data = self._data[...]
        mtx.indices = self._indices[...]
        mtx.indptr = self._indptr
        return mtx

    async def to_memory_async(self) -> CSMatrix | CSArray:
        return await self.getitem(())


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
    msg = f"Unknown encoding type {encoding_type}"
    raise ValueError(msg)


@_subset.register(BaseCompressedSparseDataset)
def subset_sparsedataset(d, subset_idx):
    return d[subset_idx]
