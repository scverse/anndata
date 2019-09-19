# TODO:
# - think about supporting the COO format
from collections.abc import Iterable, Mapping
from typing import Union, KeysView, NamedTuple

import numpy as np
import scipy.sparse as ss
from scipy.sparse import _sparsetools

from ..h5py.utils import _chunked_rows
from ..utils import unpack_index


def create_formats(dataset_class, sparse_dataset_wrapper_class):
    class BackedSparseMatrixMixin:
        """\
        Mixin class for backed sparse matrices.

        Largely needed for the case `backed_sparse_csr(...)[:]`, since that calls
        copy on `.data`, `.indices`, and `.indptr`.
        """

        def copy(self):
            if isinstance(self.data, dataset_class):
                return sparse_dataset_wrapper_class(self.data.parent).value
            else:
                return super().copy()

    class backed_csr_matrix(BackedSparseMatrixMixin, ss.csr_matrix):
        pass

    class backed_csc_matrix(BackedSparseMatrixMixin, ss.csc_matrix):
        pass

    backed_csr_matrix._set_many = _set_many
    backed_csc_matrix._set_many = _set_many
    backed_csr_matrix._zero_many = _zero_many
    backed_csc_matrix._zero_many = _zero_many

    return [
        BackedFormat("csr_matrix", backed_csr_matrix, ss.csr_matrix),
        BackedFormat("csr", backed_csr_matrix, ss.csr_matrix),
        BackedFormat("csc", backed_csc_matrix, ss.csc_matrix),
        BackedFormat("csc_matrix", backed_csc_matrix, ss.csc_matrix),
    ]


class BackedFormat(NamedTuple):
    format_str: str
    backed_type: type
    memory_type: type


def add_mapping_impl(
    dataset_impl_class,
    group_impl_class,
    X_encoding_attr_key,
    X_shape_attr_key,
    group_wrapper_class,
    sparse_dataset_wrapper_class,
    additional_wrapper_classes,
):
    formats = create_formats(dataset_impl_class, sparse_dataset_wrapper_class)
    wrapper_classes = [
        sparse_dataset_wrapper_class,
        group_wrapper_class,
    ] + additional_wrapper_classes
    for c in wrapper_classes:
        c.FORMATS = formats
        c.X_ENCODING_FORMAT_ATTR = X_encoding_attr_key
        c.X_SHAPE_ATTR = X_shape_attr_key
        c.GROUP_WRAPPER_CLASS = group_wrapper_class
        c.DATASET_CLASS = dataset_impl_class
        c.GROUP_CLASS = group_impl_class
        c.SPARSE_DATASET_WRAPPER_CLASS = sparse_dataset_wrapper_class


def get_format_str(data, formats):
    for fmt, backed_class, memory_class in formats:
        if isinstance(data, memory_class):
            return fmt
    raise ValueError(f"Data type {type(data)} is not supported.")


def get_memory_class(format_str, formats):
    for fmt, backed_class, memory_class in formats:
        if format_str == fmt:
            return memory_class
    raise ValueError(f"Format string {format_str} is not supported.")


def get_backed_class(format_str, formats):
    for fmt, backed_class, memory_class in formats:
        if format_str == fmt:
            return backed_class
    raise ValueError(f"Format string {format_str} is not supported.")


# def _load_dataset_as_sparse(sds, chunk_size=6000):
#     # efficient for csr, not so for csc (but still better than loompy it seems)
#     if 'sparse_format' in sds.attrs:
#         sparse_class = get_memory_class(sds.attrs['sparse_format'])
#     else:
#         sparse_class = ss.csr_matrix
#
#     data = None
#
#     for chunk, _, _ in _chunked_rows(sds, chunk_size, True):
#         data = (
#                 sparse_class(chunk)
#                 if data is None
#                 else ss.vstack([data, sparse_class(chunk)])
#         )
#
#     return data


class Group(Mapping):
    """\
    Like :class:`h5py.Group <h5py:Group>`, but able to handle sparse matrices.
    """

    def __init__(self, group, force_dense=False):
        self.group = group
        self.force_dense = force_dense

    def __iter__(self):
        for k in self.keys():
            yield k

    def __len__(self):
        return len(self.group)

    def __getitem__(
        self, key: str
    ) -> Union[
        'h5py.Group',
        'h5py.Dataset',
        'zarr.core.Array',
        'zarr.hierarchy.Group',
        'SparseDataset',
    ]:
        item = self.group[key]
        if isinstance(item, self.__class__.GROUP_CLASS):
            encoding = item.attrs.get(self.__class__.X_ENCODING_FORMAT_ATTR, '')
            if (
                encoding == 'csr_matrix'
                or encoding == 'csc_matrix'
                or encoding == 'csr'
                or encoding == 'csc'
            ):
                return self.__class__.SPARSE_DATASET_WRAPPER_CLASS(item)
            else:
                return self.__class__.GROUP_WRAPPER_CLASS(item)
        elif isinstance(item, self.__class__.DATASET_CLASS):
            return item
        else:
            raise ValueError("Unexpected item type.")

    def create_group(self, name):
        return self.__class__.GROUP_WRAPPER_CLASS(self.group.create_group(name))

    @property
    def attrs(self):
        return self.group.attrs

    @property
    def name(self):
        return self.group.name

    def __delitem__(self, name):
        self.group.__delitem__(name)

    def __setitem__(
        self,
        key: str,
        value: Union[
            'h5py.Group',
            'h5py.Dataset',
            'zarr.hierarchy.Group',
            'zarr.core.Array',
            'SparseDataset',
        ],
    ):
        self.group.__setitem__(key, value)

    def keys(self) -> KeysView[str]:
        return self.group.keys()

    def create_dataset(self, name, data=None, chunk_size=6000, **kwargs):
        if data is None:
            raise NotImplementedError(
                "`create_dataset` is only supported if `data` is passed."
            )
        if not isinstance(data, SparseDataset) and not ss.issparse(data):
            return self.group.create_dataset(name=name, data=data, **kwargs)
        if self.force_dense:
            sds = self.group.create_dataset(
                name=name, shape=data.shape, dtype=data.dtype, **kwargs
            )
            for chunk, start, end in _chunked_rows(data, chunk_size):
                sds[start:end] = chunk.toarray()
            sds.attrs['sparse_format'] = (
                data.format_str
                if isinstance(data, SparseDataset)
                else get_format_str(data, self.__class__.FORMATS)
            )
            return sds
        # SparseDataset or spmatrix
        group = self.group.create_group(name)
        if isinstance(data, SparseDataset):
            for attr in [
                self.__class__.X_ENCODING_FORMAT_ATTR,
                self.__class__.X_SHAPE_ATTR,
            ]:
                group.attrs[attr] = data.group.attrs[attr]
            get_dataset = data.group.__getitem__
        else:  # ss.issparse(data):
            group.attrs[self.__class__.X_ENCODING_FORMAT_ATTR] = get_format_str(
                data, self.__class__.FORMATS
            )
            group.attrs[self.__class__.X_SHAPE_ATTR] = data.shape
            get_dataset = lambda d: getattr(data, d)
        for dataset in ['data', 'indices', 'indptr']:
            group.create_dataset(
                dataset, data=get_dataset(dataset), maxshape=(None,), **kwargs
            )
        return self.__class__.SPARSE_DATASET_WRAPPER_CLASS(group)


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
            'Currently, you cannot change the sparsity structure of a AbstractSparseDataset.'
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


class SparseDataset:
    """\
    Analogous to :class:`h5py.Dataset <h5py:Dataset>`, but for sparse matrices.
    """

    def __init__(self, group):
        self.group = group

    @property
    def value(self):
        pass

    @property
    def format_str(self):
        return self.group.attrs[self.__class__.X_ENCODING_FORMAT_ATTR]

    def _get_mock_matrix(self, index):
        if index == ():
            index = slice(None)
        row, col = unpack_index(index)
        if all(isinstance(x, Iterable) for x in (row, col)):
            row, col = np.ix_(row, col)
        format_class = get_backed_class(self.format_str, self.__class__.FORMATS)
        mock_matrix = format_class(self.shape, dtype=self.dtype)
        mock_matrix.data = self.group['data']
        mock_matrix.indices = self.group['indices']
        mock_matrix.indptr = self.group['indptr']
        return mock_matrix, row, col

    def __getitem__(self, index):
        mock_matrix, row, col = self._get_mock_matrix(index)
        return mock_matrix[row, col]

    def __setitem__(self, index, value):
        mock_matrix, row, col = self._get_mock_matrix(index)
        mock_matrix[row, col] = value

    @property
    def shape(self):
        return tuple(self.group.attrs[self.__class__.X_SHAPE_ATTR])

    @property
    def dtype(self):
        return self.group['data'].dtype

    def append(self, sparse_matrix):
        shape = self.group.attrs[self.__class__.X_SHAPE_ATTR]

        if self.format_str != get_format_str(
            sparse_matrix, self.__class__.FORMATS
        ):
            raise ValueError("Format not the same.")

        if self.format_str != 'csr' and self.format_str != 'csr_matrix':
            raise NotImplementedError(
                f"The append method for format {self.format_str} "
                f"is not implemented."
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
        self.group.attrs[self.__class__.X_SHAPE_ATTR] = (
            shape[0] + sparse_matrix.shape[0],
            max(shape[1], sparse_matrix.shape[1]),
        )
