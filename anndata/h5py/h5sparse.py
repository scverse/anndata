from typing import Optional

import six
import h5py
import numpy as np
import scipy.sparse as ss

from .utils import _chunked_rows
from ..compat import PathLike

from .. import backwards

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


class Group:
    """Like :ref:`h5py.Group <h5py:Group>`, but able to handle sparse matrices.
    """

    def __init__(self, h5py_group):
        self.h5py_group = h5py_group

    def __getitem__(self, key):
        h5py_item = self.h5py_group[key]
        if isinstance(h5py_item, h5py.Group):
            if 'h5sparse_format' in h5py_item.attrs:
                # detect the sparse matrix
                return backwards.SparseDataset(h5py_item)
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

    def create_dataset(self, name, data=None, chunk_size=6000, **kwargs):
        if data is None:

            raise NotImplementedError("Only support create_dataset with "
                                      "if `data` is passed.")
        elif ss.issparse(data):

            sds = self.h5py_group.create_dataset(name=name, shape=data.shape, dtype=data.dtype, **kwargs)
            for chunk, start, end in _chunked_rows(data, chunk_size):
                sds[start:end] = chunk.toarray()
            sds.attrs['sparse_format'] = get_format_str(data)
            return sds

        else:
            return self.h5py_group.create_dataset(
                name=name, data=data, **kwargs)


Group.create_dataset.__doc__ = h5py.Group.create_dataset.__doc__


class File(Group):
    """Like :ref:`h5py.File <h5py:File>`, but able to handle sparse matrices.
    """

    def __init__(
        self, name: PathLike,
        mode: Optional[str] = None,
        driver: Optional[str] = None,
        libver: Optional[str] = None,
        userblock_size: Optional[int] = None,
        swmr: bool = False,
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
        super().__init__(self.h5f)

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


def _load_h5_dataset_as_sparse(sds, chunk_size=6000):
    #efficient for csr, not so for csc (but still better than loompy it seems)
    if not isinstance(sds, h5py.Dataset):
        raise ValueError('sds should be a h5py Dataset')

    if 'sparse_format' in sds.attrs:
        sparse_class = get_format_class(sds.attrs['sparse_format'])
    else:
        sparse_class = ss.csr_matrix

    data = None

    for chunk, _, _ in _chunked_rows(sds, chunk_size):
        data = sparse_class(chunk) if data is None else ss.vstack([data, sparse_class(chunk)])

    return data

def _load_h5_dataset_as_sparse2(sds, chunk_size=512):
    #Like loompy does, not efficient it seems
    if not isinstance(sds, h5py.Dataset):
        raise ValueError('sds should be a h5py Dataset')

    elif 'sparse_format' not in sds.attrs:
        raise ValueError('sds should have the sparse_format attribute')

    else:
        sparse_class = get_format_class(sds.attrs['sparse_format'])

        data = []
        row = []
        col = []

        for vals, start, _ in _chunked_rows(sds, chunk_size):
            nonzeros = np.where(vals > 0)
            data.append(vals[nonzeros])
            row.append(nonzeros[0] + start)
            col.append(nonzeros[1])

        return sparse_class((np.concatenate(data), (np.concatenate(row), np.concatenate(col))),
                            shape=sds.shape)
