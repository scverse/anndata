import numpy as np
import scipy.sparse as sparse


from .. import mapping as adm
from .zarrsparse import Dataset, File
from ..compat import _from_fixed_length_strings


def read_zarr(filename, backed=None, chunk_size=None):
    return Reader().read(filename, backed=backed)


class Reader(adm.MappingReader):
    def __init__(self):
        super().__init__(Dataset)

    def read_csr(self, group: 'zarr.Group') -> sparse.csr_matrix:
        return sparse.csr_matrix(
            (group["data"], group["indices"], group["indptr"]),
            shape=group.attrs["shape"],
        )

    def read_group(self, group):
        if "encoding-type" in group.attrs:
            enctype = group.attrs["encoding-type"]
            if enctype == "dataframe":
                return self.read_dataframe(group)
            elif enctype == "csr_matrix":
                return self.read_csr(group)
            elif enctype == "csc_matrix":
                return self.read_csc(group)

        return {k: self.read_attribute(group[k]) for k in group.keys()}

    def read_csc(self, group: 'zarr.Group') -> sparse.csc_matrix:
        return sparse.csc_matrix(
            (group["data"], group["indices"], group["indptr"]),
            shape=group.attrs["shape"],
        )

    def read_dataset(self, dataset):
        value = dataset[...]
        if not hasattr(value, "dtype"):
            return value
        elif isinstance(value.dtype, str):
            pass
        elif issubclass(value.dtype.type, np.str_):
            value = value.astype(object)
        elif issubclass(value.dtype.type, np.string_):
            value = value.astype(str).astype(
                object
            )  # bytestring -> unicode -> str
        elif len(value.dtype.descr) > 1:  # Compound dtype
            # For backwards compat, now strings are written as variable length
            value = _from_fixed_length_strings(value)
        if value.shape == ():
            value = value[()]
        return value

    def create_file(self, filename, mode):
        return File(filename, mode=mode)

    def close_file(self, f):
        pass
