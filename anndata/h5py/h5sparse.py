# TODO:
# - think about making all of the below subclasses
# - think about supporting the COO format
from typing import Optional, Union, KeysView, NamedTuple

import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

from ..utils import unpack_index
from ..compat import PathLike

from .utils import _chunked_rows


class BackedFormat(NamedTuple):
    format_str: str
    backed_type: type
    memory_type: type


class backed_csr_matrix(ss.csr_matrix):
    pass


class backed_csc_matrix(ss.csc_matrix):
    pass


FORMATS = [
    BackedFormat("csr", backed_csr_matrix, ss.csr_matrix),
    BackedFormat("csc", backed_csc_matrix, ss.csc_matrix),
]


def get_format_str(data):
    for fmt, backed_class, memory_class in FORMATS:
        if isinstance(data, memory_class):
            return fmt
    raise ValueError("Data type {} is not supported.".format(type(data)))


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
        data = sparse_class(chunk) if data is None else ss.vstack([data, sparse_class(chunk)])

    return data


class Group:
    """Like :class:`h5py.Group <h5py:Group>`, but able to handle sparse matrices.
    """

    def __init__(self, h5py_group, force_dense=False):
        self.h5py_group = h5py_group
        self.force_dense = force_dense

    def __getitem__(self, key: str) -> Union[h5py.Group, h5py.Dataset, 'SparseDataset']:
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

    def __delitem__(self, name):
        self.h5py_group.__delitem__(name)

    def __setitem__(self, key: str, value: Union[h5py.Group, h5py.Dataset, 'SparseDataset']):
        self.h5py_group.__setitem__(key, value)

    def keys(self) -> KeysView[str]:
        return self.h5py_group.keys()

    def create_dataset(self, name, data=None, chunk_size=6000, **kwargs):
        if data is None:
            raise NotImplementedError("Only support create_dataset with "
                                      "if `data` is passed.")

        elif (isinstance(data, SparseDataset) or ss.issparse(data)) and self.force_dense:
            sds = self.h5py_group.create_dataset(name=name, shape=data.shape, dtype=data.dtype, **kwargs)
            for chunk, start, end in _chunked_rows(data, chunk_size):
                sds[start:end] = chunk.toarray()
            sds.attrs['sparse_format'] = data.format_str if isinstance(data, SparseDataset) else get_format_str(data)
            return sds
        elif isinstance(data, SparseDataset) and not self.force_dense:
            group = self.h5py_group.create_group(name)
            group.attrs['h5sparse_format'] = data.h5py_group.attrs['h5sparse_format']
            group.attrs['h5sparse_shape'] = data.h5py_group.attrs['h5sparse_shape']
            group.create_dataset('data', data=data.h5py_group['data'], maxshape=(None,), **kwargs)
            group.create_dataset('indices', data=data.h5py_group['indices'], maxshape=(None,), **kwargs)
            group.create_dataset('indptr', data=data.h5py_group['indptr'], maxshape=(None,), **kwargs)
            return SparseDataset(group)
        elif ss.issparse(data) and not self.force_dense:
            group = self.h5py_group.create_group(name)
            group.attrs['h5sparse_format'] = get_format_str(data)
            group.attrs['h5sparse_shape'] = data.shape
            group.create_dataset('data', data=data.data, maxshape=(None,), **kwargs)
            group.create_dataset('indices', data=data.indices, maxshape=(None,), **kwargs)
            group.create_dataset('indptr', data=data.indptr, maxshape=(None,), **kwargs)
            return SparseDataset(group)
        else:
            return self.h5py_group.create_dataset(
                name=name, data=data, **kwargs)


Group.create_dataset.__doc__ = h5py.Group.create_dataset.__doc__


class File(Group):
    """Like :class:`h5py.File <h5py:File>`, but able to handle sparse matrices.
    """

    def __init__(
        self, name: PathLike,
        mode: Optional[str] = None,
        driver: Optional[str] = None,
        libver: Optional[str] = None,
        userblock_size: Optional[int] = None,
        swmr: bool = False,
        force_dense: bool = False,
        **kwds  # Python 3.5 can’t handle trailing commas here
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
    """Sets value at each (i, j) to x

    Here (i,j) index major and minor respectively, and must not contain
    duplicate entries.
    """
    i, j, M, N = self._prepare_indices(i, j)

    if np.isscalar(x):  # Scipy 1.3+ compat
        n_samples = 1
    else:
        n_samples = len(x)
    offsets = np.empty(n_samples, dtype=self.indices.dtype)
    ret = _sparsetools.csr_sample_offsets(M, N, self.indptr, self.indices,
                                          n_samples, i, j, offsets)
    if ret == 1:
        # rinse and repeat
        self.sum_duplicates()
        _sparsetools.csr_sample_offsets(M, N, self.indptr,
                                        self.indices, n_samples, i, j,
                                        offsets)

    if -1 not in offsets:
        # make a list for interaction with h5py
        offsets = list(offsets)
        # only affects existing non-zero cells
        self.data[offsets] = x
        return

    else:
        raise ValueError(
            'Currently, you cannot change the sparsity structure of a SparseDataset.')
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
    """Sets value at each (i, j) to zero, preserving sparsity structure.

    Here (i,j) index major and minor respectively.
    """
    i, j, M, N = self._prepare_indices(i, j)

    n_samples = len(i)
    offsets = np.empty(n_samples, dtype=self.indices.dtype)
    ret = _sparsetools.csr_sample_offsets(M, N, self.indptr, self.indices,
                                          n_samples, i, j, offsets)
    if ret == 1:
        # rinse and repeat
        self.sum_duplicates()
        _sparsetools.csr_sample_offsets(M, N, self.indptr,
                                        self.indices, n_samples, i, j,
                                        offsets)

    # only assign zeros to the existing sparsity structure
    self.data[list(offsets[offsets > -1])] = 0


backed_csr_matrix._zero_many = _zero_many
backed_csc_matrix._zero_many = _zero_many


class SparseDataset:
    """Analogous to :class:`h5py.Dataset <h5py:Dataset>`, but for sparse matrices.
    """

    def __init__(self, h5py_group):
        self.h5py_group = h5py_group

    def __repr__(self):
        return ('<HDF5 sparse dataset: format \'{}\', shape {}, type \'{}\'>'
                .format(
                    self.format_str,
                    self.shape, self.h5py_group['data'].dtype.str))

    @property
    def format_str(self):
        return self.h5py_group.attrs['h5sparse_format']

    def __getitem__(self, index):
        if index == (): index = slice(None)
        row, col = unpack_index(index)
        format_class = get_backed_class(self.format_str)
        mock_matrix = format_class(self.shape, dtype=self.dtype)
        mock_matrix.data = self.h5py_group['data']
        mock_matrix.indices = self.h5py_group['indices']
        mock_matrix.indptr = self.h5py_group['indptr']
        return mock_matrix[row, col]

    def __setitem__(self, index, value):
        if index == (): index = slice(None)
        row, col = unpack_index(index)
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
        object = self.h5py_group
        data_array = format_class(self.shape, dtype=self.dtype)
        data_array.data = np.empty(object['data'].shape, object['data'].dtype)
        data_array.indices = np.empty(object['indices'].shape, object['indices'].dtype)
        data_array.indptr = np.empty(object['indptr'].shape, object['indptr'].dtype)
        object['data'].read_direct(data_array.data)
        object['indices'].read_direct(data_array.indices)
        object['indptr'].read_direct(data_array.indptr)
        return data_array

    def append(self, sparse_matrix):
        shape = self.h5py_group.attrs['h5sparse_shape']

        if self.format_str != get_format_str(sparse_matrix):
            raise ValueError("Format not the same.")

        if self.format_str == 'csr':
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
            indptr[orig_data_size:] = (sparse_matrix.indptr[1:].astype(np.int64)
                                       + append_offset)

            # indices
            indices = self.h5py_group['indices']
            orig_data_size = indices.shape[0]
            new_shape = (orig_data_size + sparse_matrix.indices.shape[0],)
            indices.resize(new_shape)
            indices[orig_data_size:] = sparse_matrix.indices

            # shape
            self.h5py_group.attrs['h5sparse_shape'] = (
                shape[0] + sparse_matrix.shape[0],
                max(shape[1], sparse_matrix.shape[1]))
        else:
            raise NotImplementedError("The append method for format {} is not "
                                      "implemented.".format(self.format_str))
