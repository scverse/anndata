# TODO:
# - think about making all of the below subclasses
# - think about supporting the COO format

import six
import h5py
import numpy as np
import scipy.sparse as ss
from scipy.sparse.sputils import IndexMixin

FORMAT_DICT = {
    'csr': ss.csr_matrix,
    'csc': ss.csc_matrix,
}


def get_format_str(data):
    for format_str, format_class in six.viewitems(FORMAT_DICT):
        if isinstance(data, format_class):
            return format_str
    raise ValueError("Data type {} is not supported.".format(type(data)))


def get_format_class(format_str):
    format_class = FORMAT_DICT.get(format_str, None)
    if format_class is None:
        raise ValueError("Format string {} is not supported."
                         .format(format_str))
    return format_class


class Group(object):
    """Like `h5py.Group`, but able to handle sparse matrices.
    """

    def __init__(self, h5py_group):
        self.h5py_group = h5py_group

    def __getitem__(self, key):
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

    def __setitem__(self, key, value):
        self.h5py_group.__setitem__(key, value)

    def keys(self):
        return self.h5py_group.keys()

    def create_dataset(self, name, data=None, **kwargs):
        if data is None:
            raise NotImplementedError("Only support create_dataset with "
                                      "if `data` is passed.")
        elif isinstance(data, SparseDataset):
            group = self.h5py_group.create_group(name)
            group.attrs['h5sparse_format'] = data.h5py_group.attrs['h5sparse_format']
            group.attrs['h5sparse_shape'] = data.h5py_group.attrs['h5sparse_shape']
            group.create_dataset('data', data=data.h5py_group['data'], maxshape=(None,), **kwargs)
            group.create_dataset('indices', data=data.h5py_group['indices'], maxshape=(None,), **kwargs)
            group.create_dataset('indptr', data=data.h5py_group['indptr'], maxshape=(None,), **kwargs)
            return SparseDataset(group)
        elif ss.issparse(data):
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


class File(Group):
    """Like `h5py.File`, but able to handle sparse matrices.
    """

    def __init__(self, *args, **kwargs):
        self.h5f = h5py.File(*args, **kwargs)
        self.h5py_group = self.h5f

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


from scipy.sparse.compressed import _cs_matrix
from scipy.sparse import _sparsetools


def _set_many(self, i, j, x):
    """Sets value at each (i, j) to x

    Here (i,j) index major and minor respectively, and must not contain
    duplicate entries.
    """
    i, j, M, N = self._prepare_indices(i, j)

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
        # raise ValueError(
        #     'Currently, you cannot change the sparsity structure of a SparseDataset.')
        # replace where possible
        mask = offsets > -1
        self.data[offsets[mask]] = x[mask]
        # only insertions remain
        mask = ~mask
        i = i[mask]
        i[i < 0] += M
        j = j[mask]
        j[j < 0] += N
        self._insert_many(i, j, x[mask])

_cs_matrix._set_many = _set_many


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
        
_cs_matrix._zero_many = _zero_many

    
class SparseDataset(IndexMixin):
    """Analogous to `h5py.Dataset`, but for sparse matrices.
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
        row, col = self._unpack_index(index)
        format_class = get_format_class(self.format_str)
        mock_matrix = format_class(self.shape, dtype=self.dtype)
        mock_matrix.data = self.h5py_group['data']
        mock_matrix.indices = self.h5py_group['indices']
        mock_matrix.indptr = self.h5py_group['indptr']
        return mock_matrix[row, col]

    def __setitem__(self, index, value):
        if index == (): index = slice(None)
        row, col = self._unpack_index(index)
        format_class = get_format_class(self.format_str)
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
        data = self.h5py_group['data'].value
        indices = self.h5py_group['indices'].value
        indptr = self.h5py_group['indptr'].value
        shape = self.h5py_group.attrs['h5sparse_shape']
        format_class = get_format_class(self.format_str)
        return format_class((data, indices, indptr), shape=shape)

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
