from os import PathLike
from typing import Optional

import zarr
from .. import mapping as adm


Dataset = zarr.core.Array


class Group(adm.Group):
    def __init__(self, group, force_dense=False):
        super().__init__(group, force_dense)


class SparseDataset(adm.SparseDataset):
    def __init__(self, group):
        super().__init__(group)

    @property
    def value(self):
        format_class = adm.get_memory_class(
            self.format_str, self.__class__.FORMATS
        )
        g = self.group
        data_array = format_class(self.shape, dtype=self.dtype)
        data_array.data = g['data'][()]
        data_array.indices = g['indices'][()]
        data_array.indptr = g['indptr'][()]
        return data_array


class File(Group):
    def __init__(
        self,
        name: PathLike,
        mode: Optional[str] = None,
        force_dense: bool = False,
        **kwds,
    ):
        self.root = zarr.open(str(name), mode=mode)
        self.fake_id = 1
        super().__init__(self.root, force_dense)

    def __enter__(self):
        return self

    def close(self):
        self.fake_id = None

    @property
    def id(self):
        return self.fake_id

    @property
    def filename(self):
        return self.root.store.path


adm.add_mapping_impl(
    zarr.core.Array,
    zarr.hierarchy.Group,
    'encoding-type',
    'shape',
    Group,
    SparseDataset,
    [File, Dataset],
)
