#thats just start

from .base import AnnData
from .h5py import _load_h5_dataset_as_sparse
import scipy.sparse as ss

# l = AnnDataLayer(adata, 'obsm', 'principal_components')
class AnnDataLayer():

    def __init__(adata: AnnData, atype: str, name: str):
        self._adata = adata
        self._type = atype
        self._name = name

    def __getitem__(self, slice):
        if self._adata.isbacked:
            return self._adata.file._file['/'+self._type+'/'+self._name][slice]
        else:
            return getattr(self._adata, self._type)[self._name][slice]

    def __setitem__(self, slice, data):
        if self._adata.isbacked:
            self._adata.file._file['/'+self._type+'/'+self._name][slice] = data
        else:
            getattr(self._adata, self._type)[self._name][slice] = data

    def sparse(self):
        if self._adata.isbacked
            return _load_h5_dataset_as_sparse(self._adata.file._file['/'+self._type+'/'+self._name])
        else:
            return ss.csr_matrix(getattr(self._adata, self._type)[self._name])
