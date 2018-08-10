#thats just start
import numpy as np
from collections import OrderedDict
from scipy.sparse import issparse

class AnnDataLayers():

    def __init__(self, adata, layers=None, dtype='float32', adata_ref=None, oidx=None, vidx=None):

        self._adata = adata
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx

        self._layers = OrderedDict()

        if layers is not None:
            for key in layers.keys():
                if layers[key].shape != self._adata.shape:
                    raise ValueError('Shape does not fit.')
                if layers[key].dtype != np.dtype(dtype):
                    self._layers[key] = layers[key].astype(dtype, copy=False)
                else:
                    self._layers[key] = layers[key]

    def __getitem__(self, key):
        if self.isview:
            return self._adata_ref.layers[key][self._oidx, self._vidx]
        else:
            return self._layers[key]

    def __setitem__(self, key, value):
        if not isinstance(value, np.ndarray) and not issparse(value):
            raise ValueError('Value should be numpy array or sparse matrix')
        if self.isview:
            if key not in self._adata_ref.layers.keys():
                # can't load _init_actual_AnnData from .base - circular
                if self._adata.isbacked:
                    raise ValueError(
                        'You cannot modify elements of an AnnData view, '
                        'but need a copy of the subset.\n\n'
                        'Call `adata_subset = adata[index].copy(filename=...)`.')
                self._adata._init_as_actual(self._adata.copy())
                self._adata.layers[key] = value
            else:
                self._adata_ref.layers[key][self._oidx, self._vidx] = value
        else:
            if value.shape != self._adata.shape:
                raise ValueError('Shape does not fit.')
            self._layers[key] = value

    def __delitem__(self, key):
        self.__delattr__(key)

    def __delattr__(self, key):
        if self.isview:
            adata = self._adata.copy()
            adata.layers.__delattr__(key)
            self._adata._init_as_actual(adata)
        else:
            del self._layers[key]

    def keys(self):
        if self.isview:
            return self._adata_ref.layers.keys()
        else:
            return list(self._layers.keys())

    def items(self):
        if self.isview:
            return [(k, v[self._oidx, self._vidx].copy()) for (k, v) in self._adata_ref.layers.items()]
        else:
            return [(k, v.copy()) for (k, v) in self._layers.items()]

    def __len__(self):
        if self.isview:
            return len(self._adata_ref.layers)
        else:
            return len(self._layers)

    def as_dict(self):
        return {k:v for (k, v) in self.items()}

    @property
    def isview(self):
        return self._adata_ref is not None
