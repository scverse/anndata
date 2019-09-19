from os import PathLike
from typing import Optional

import h5py
import numpy as np

import anndata.mapping as adm

Dataset = h5py.Dataset


class Group(adm.Group):

    def __init__(self, group, force_dense=False):
        super().__init__(group, force_dense)


class SparseDataset(adm.SparseDataset):
    def __init__(self, group):
        super().__init__(group)

    @property
    def value(self):
        format_class = adm.get_memory_class(self.format_str, self.__class__.FORMATS)
        object = self.group
        data_array = format_class(self.shape, dtype=self.dtype)
        data_array.data = np.empty(object['data'].shape, object['data'].dtype)
        data_array.indices = np.empty(
            object['indices'].shape, object['indices'].dtype
        )
        data_array.indptr = np.empty(
            object['indptr'].shape, object['indptr'].dtype
        )
        object['data'].read_direct(data_array.data)
        object['indices'].read_direct(data_array.indices)
        object['indptr'].read_direct(data_array.indptr)
        return data_array


class File(Group):
    """\
    Like :class:`h5py.File <h5py:File>`, but able to handle sparse matrices.
    """

    def __init__(
            self,
            name: PathLike,
            mode: Optional[str] = None,
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
            swmr=swmr
        )
        super().__init__(self.h5f, force_dense)

    def __enter__(self):
        return self

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        if self._file is not None:
            self._file.close()

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

adm.add_mapping_impl(h5py.Dataset, h5py.Group, 'h5sparse_format', 'h5sparse_shape', Group, SparseDataset, [File,
                                                                                                           Dataset])

Group.create_dataset.__doc__ = h5py.Group.create_dataset.__doc__
Group.create_group.__doc__ = h5py.Group.create_group.__doc__
