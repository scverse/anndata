import anndata.mapping as adm

import anndata.h5py as ah5py
import os


def read_h5ad(filename, backed=None, chunk_size=None):
    return Reader().read(filename, backed=backed)


class Reader(adm.MappingReader):

    def __init__(self):
        super().__init__(ah5py.Dataset)

    def create_file(self, filename, mode):
        if mode is None:
            mode = 'w' if os.path.exists(filename) else 'r'
        return ah5py.File(filename, mode=mode)

    def close_file(self, f):
        f.close()
