#thats just start
import numpy as np
from collections import OrderedDict

from scipy.sparse import issparse

class AnnDataLayers():

    def __init__(self, adata, layers=None, adata_ref=None, oidx=None, vidx=None):

        self._adata = adata
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx

        if layers is not None:
            for key in layers.keys():
                if layers[key].shape != self._adata.shape:
                    raise ValueError('Shape does not fit.')

        self._layers = OrderedDict(layers) if layers is not None else OrderedDict()

    def __getitem__(self, key):
        if self.isview:
            return self._adata_ref.layers_X[key][self._oidx, self._vidx]
        else:
            return self._layers[key]

    def __setitem__(self, key, value):
        if not isinstance(value, np.ndarray) and not issparse(value):
            raise ValueError('Value should be numpy array or sparse matrix')
        if self.isview:
            if key not in self._adata_ref.layers_X.keys():
                # can't load _init_actual_AnnData from .base - circular
                if self._adata.isbacked:
                    raise ValueError(
                        'You cannot modify elements of an AnnData view, '
                        'but need a copy of the subset.\n\n'
                        'Call `adata_subset = adata[index].copy(filename=...)`.')
                self._adata._init_as_actual(self._adata.copy())
                self._adata.layers_X[key] = value
            else:
                self._adata_ref.layers_X[key][self._oidx, self._vidx] = value
        else:
            if value.shape != self._adata.shape:
                raise ValueError('Shape does not fit.')
            self._layers[key] = value

    def keys(self):
        if self.isview:
            return self._adata_ref.layers_X.keys()
        else:
            return self._layers.keys()

    def items(self):
        if self.isview:
            return [(k, v[self._oidx, self._vidx].copy()) for (k, v) in self._adata_ref.layers_X.items()]
        else:
            return [(k, v.copy()) for (k, v) in self._layers.items()]

    def __len__(self):
        if self.isview:
            return len(self._adata_ref.layers_X)
        else:
            return len(self._layers)

    def as_dict(self):
        return {k:v for (k, v) in self.items()}

    @property
    def isview(self):
        return self._adata_ref is not None
