# TODO:
# - think about making all of the below subclasses
# - think about supporting the COO format
from collections.abc import Iterable, Mapping
from itertools import accumulate, chain
from os import PathLike
from typing import Optional, Union, KeysView, NamedTuple, Tuple, Sequence

import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

from ..utils import unpack_index
from ..compat import Literal
from .utils import _chunked_rows


class BackedFormat(NamedTuple):
    format_str: str
    backed_type: Type['BackedSparseMatrixMixin']
    memory_type: Type[ss.spmatrix]


class BackedSparseMatrixMixin:
    """\
    Mixin class for backed sparse matrices.

    Largely needed for the case `backed_sparse_csr(...)[:]`, since that calls
    copy on `.data`, `.indices`, and `.indptr`.
    """

    def copy(self):
        if isinstance(self.data, h5py.Dataset):
            return SparseDataset(self.data.parent).value
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


def _load_h5_dataset_as_sparse(sds, chunk_size=6000):
    # efficient for csr, not so for csc (but still better than loompy it seems)
    if not isinstance(sds, h5py.Dataset):
        raise ValueError('sds should be a h5py Dataset')

    if 'sparse_format' in sds.attrs:
        sparse_class = get_memory_class(sds.attrs['sparse_format'])
    else:
        sparse_class = ss.csr_matrix

    data = None

    for chunk, _, _ in _chunked_rows(sds, chunk_size, True):
        data = (
            sparse_class(chunk)
            if data is None
            else ss.vstack([data, sparse_class(chunk)])
        )

    return data


class Group(Mapping):
    """\
    Like :class:`h5py.Group <h5py:Group>`, but able to handle sparse matrices.
    """

    def __init__(self, h5py_group, force_dense=False):
        self.h5py_group = h5py_group
        self.force_dense = force_dense

    def __iter__(self):
        for k in self.keys():
            yield k

    def __len__(self):
        return len(self.h5py_group)

    def __getitem__(
        self, key: str
    ) -> Union[h5py.Group, h5py.Dataset, 'SparseDataset']:
        h5py_item = self.h5py_group[key]
        if isinstance(h5py_item, h5py.Group):
            if 'h5sparse_format' in h5py_item.attrs:
                # detect the sparse matrix
                return SparseDataset(h5py_item)
            else:
                return Group(h5py_item)
        elif isinstance(h5py_item, h5py.Dataset):
            return h5py_item
        else:
            raise ValueError("Unexpected item type.")

    @property
    def attrs(self):
        return self.h5py_group.attrs

    @property
    def name(self):
        return self.h5py_group.name

    def __delitem__(self, name):
        self.h5py_group.__delitem__(name)

    def __setitem__(
        self, key: str, value: Union[h5py.Group, h5py.Dataset, 'SparseDataset']
    ):
        self.h5py_group.__setitem__(key, value)

    def keys(self) -> KeysView[str]:
        return self.h5py_group.keys()

    def create_dataset(self, name, data=None, chunk_size=6000, **kwargs):
        if data is None:
            raise NotImplementedError(
                "`create_dataset` is only supported if `data` is passed."
            )
        if not isinstance(data, SparseDataset) and not ss.issparse(data):
            return self.h5py_group.create_dataset(
                name=name, data=data, **kwargs
            )
        if self.force_dense:
            sds = self.h5py_group.create_dataset(
                name=name, shape=data.shape, dtype=data.dtype, **kwargs
            )
            for chunk, start, end in _chunked_rows(data, chunk_size):
                sds[start:end] = chunk.toarray()
            sds.attrs['sparse_format'] = (
                data.format_str
                if isinstance(data, SparseDataset)
                else get_format_str(data)
            )
            return sds
        # SparseDataset or spmatrix
        group = self.h5py_group.create_group(name)
        if isinstance(data, SparseDataset):
            for attr in ['h5sparse_format', 'h5sparse_shape']:
                group.attrs[attr] = data.h5py_group.attrs[attr]
            get_dataset = data.h5py_group.__getitem__
        else:  # ss.issparse(data):
            group.attrs['h5sparse_format'] = get_format_str(data)
            group.attrs['h5sparse_shape'] = data.shape
            get_dataset = lambda d: getattr(data, d)
        for dataset in ['data', 'indices', 'indptr']:
            group.create_dataset(
                dataset, data=get_dataset(dataset), maxshape=(None,), **kwargs
            )
        return SparseDataset(group)

    def create_group(self, name, track_order=None):
        return Group(
            self.h5py_group.create_group(name, track_order=track_order)
        )


Group.create_dataset.__doc__ = h5py.Group.create_dataset.__doc__
Group.create_group.__doc__ = h5py.Group.create_group.__doc__


class File(Group):
    """\
    Like :class:`h5py.File <h5py:File>`, but able to handle sparse matrices.
    """

    def __init__(
        self,
        name: PathLike,
        mode: Optional[Literal['r', 'r+']] = None,
        driver: Optional[str] = None,
        libver: Optional[str] = None,
        userblock_size: Optional[int] = None,
        swmr: bool = False,
        force_dense: bool = False,
        **kwds,
    ):
        self.h5f = h5py.File(
            name,
            mode=mode,
            driver=driver,
            libver=libver,
            userblock_size=userblock_size,
            swmr=swmr,
            **kwds,
        )
        super().__init__(self.h5f, force_dense)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5f.__exit__(exc_type, exc_value, traceback)

    def close(self):
        self.h5f.close()

    @property
    def id(self):
        return self.h5f.id

    @property
    def filename(self):
        return self.h5f.filename


File.__init__.__doc__ = h5py.File.__init__.__doc__


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


def get_compressed_vectors(
    X: backed_sparse_matrix, row_idxs: "np.ndarray[int]"
) -> Tuple[Sequence, Sequence, Sequence]:
    slices = [slice(*(X.indptr[i : i + 2])) for i in row_idxs]
    data = np.concatenate([X.data[s] for s in slices])
    indices = np.concatenate([X.indices[s] for s in slices])
    indptr = list(accumulate(chain((0,), (s.stop - s.start for s in slices))))
    return data, indices, indptr


def get_csc_cols(X: backed_csc_matrix, idxs: np.ndarray) -> ss.csc_matrix:
    idxs = np.asarray(idxs)
    if idxs.dtype == bool:
        idxs = np.where(idxs)
    return ss.csc_matrix(
        get_compressed_vectors(X, idxs), shape=(X.shape[0], len(idxs))
    )


def get_csr_rows(X: backed_csr_matrix, idxs: np.ndarray) -> ss.csr_matrix:
    idxs = np.asarray(idxs)
    if idxs.dtype == bool:
        idxs = np.where(idxs)
    return ss.csr_matrix(
        get_compressed_vectors(X, idxs), shape=(len(idxs), X.shape[1])
    )


def get_compressed_vector(X, idx: int) -> Tuple[Sequence, Sequence, Sequence]:
    s = slice(*(X.indptr[idx : idx + 2]))
    data = X.data[s]
    indices = X.indices[s]
    indptr = [0, len(data)]
    return data, indices, indptr


def get_csc_col(X: backed_sparse_matrix, idx: int) -> ss.csc_matrix:
    return ss.csc_matrix(get_compressed_vector(X, idx), shape=(X.shape[0], 1))


def get_csr_row(X, idx: int) -> ss.csr_matrix:
    return ss.csr_matrix(get_compressed_vector(X, idx), shape=(1, X.shape[1]))


def csc_get_sliceXint(X, row: slice, col: int):
    return get_csc_col(X, col)[row, :]


def csr_get_intXslice(X, row: int, col: slice):
    return get_csr_row(X, row)[:, col]


backed_csc_matrix._get_sliceXint = csc_get_sliceXint
backed_csr_matrix._get_intXslice = csr_get_intXslice

backed_csc_matrix._get_sliceXarray = lambda x, row, col: get_csc_cols(x, col)[
    row, :
]
backed_csr_matrix._get_arrayXslice = lambda x, row, col: get_csr_rows(x, row)[
    :, col
]


class SparseDataset:
    """\
    Analogous to :class:`h5py.Dataset <h5py:Dataset>`, but for sparse matrices.
    """

    def __init__(self, h5py_group):
        self.h5py_group = h5py_group
        self.file = h5py_group.file

    @property
    def name(self):
        return self.h5py_group.name

    def __repr__(self):
        return (
            f'<HDF5 sparse dataset: format {self.format_str!r}, '
            f'shape {self.shape}, '
            f'type {self.h5py_group["data"].dtype.str!r}>'
        )

    @property
    def format_str(self):
        return self.h5py_group.attrs['h5sparse_format']

    def __getitem__(self, index):
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        format_class = get_backed_class(self.format_str)
        mock_matrix = format_class(self.shape, dtype=self.dtype)
        mock_matrix.data = self.h5py_group['data']
        mock_matrix.indices = self.h5py_group['indices']
        mock_matrix.indptr = self.h5py_group['indptr']
        return mock_matrix[row, col]

    def __setitem__(self, index, value):
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        format_class = get_backed_class(self.format_str)
        mock_matrix = format_class(self.shape, dtype=self.dtype)
        mock_matrix.data = self.h5py_group['data']
        mock_matrix.indices = self.h5py_group['indices']
        mock_matrix.indptr = self.h5py_group['indptr']
        mock_matrix[row, col] = value

    @property
    def shape(self):
        return tuple(self.h5py_group.attrs['h5sparse_shape'])

    @property
    def dtype(self):
        return self.h5py_group['data'].dtype

    @property
    def value(self):
        format_class = get_memory_class(self.format_str)
        mtx = format_class(self.shape, dtype=self.dtype)
        group = self.h5py_group
        # read_direct had trouble with empty matrices, seems to take same amount of time
        mtx.data = group["data"][...]
        mtx.indices = group["indices"][...]
        mtx.indptr = group["indptr"][...]
        return mtx

    def append(self, sparse_matrix):
        shape = self.h5py_group.attrs['h5sparse_shape']

        if self.format_str != get_format_str(sparse_matrix):
            raise ValueError("Format not the same.")

        if self.format_str != 'csr':
            raise NotImplementedError(
                f"The append method for format {self.format_str} "
                f"is not implemented."
            )

        # data
        data = self.h5py_group['data']
        orig_data_size = data.shape[0]
        new_shape = (orig_data_size + sparse_matrix.data.shape[0],)
        data.resize(new_shape)
        data[orig_data_size:] = sparse_matrix.data

        # indptr
        indptr = self.h5py_group['indptr']
        orig_data_size = indptr.shape[0]
        append_offset = indptr[-1]
        new_shape = (orig_data_size + sparse_matrix.indptr.shape[0] - 1,)
        indptr.resize(new_shape)
        indptr[orig_data_size:] = (
            sparse_matrix.indptr[1:].astype(np.int64) + append_offset
        )

        # indices
        indices = self.h5py_group['indices']
        orig_data_size = indices.shape[0]
        new_shape = (orig_data_size + sparse_matrix.indices.shape[0],)
        indices.resize(new_shape)
        indices[orig_data_size:] = sparse_matrix.indices

        # shape
        self.h5py_group.attrs['h5sparse_shape'] = (
            shape[0] + sparse_matrix.shape[0],
            max(shape[1], sparse_matrix.shape[1]),
        )
