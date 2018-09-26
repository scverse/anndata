from typing import Mapping, Optional, Union

import numpy as np
from collections import OrderedDict
from scipy.sparse import issparse

if False:
    from .base import AnnData, Index  # noqa


class AnnDataLayers:
    def __init__(
        self,
        adata: 'AnnData',
        layers: Optional[Mapping[str, np.ndarray]] = None,
        dtype: Union[str, np.dtype] = 'float32',
        adata_ref: 'AnnData' = None,
        oidx: 'Index' = None,
        vidx: 'Index' = None,
    ):
        self._adata = adata
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx

        self._layers = OrderedDict()

        for key, layer in (layers or {}).items():
            if not isinstance(layer, np.ndarray) and not issparse(layer):
                raise TypeError(
                    'Layer {} is not a numpy.ndarray or scipy.sparse.spmatrix but a {}'
                    .format(key, type(layer).__name__)
                )
            if layer.shape != self._adata.shape:
                raise ValueError('Shape does not fit.')
            if layer.dtype != np.dtype(dtype):
                self._layers[key] = layer.astype(dtype, copy=False)
            else:
                self._layers[key] = layer

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
        return (self._adata_ref.layers if self.isview else self._layers).keys()

    def items(self, copy=True):
        if self.isview:
            pairs = [(k, v[self._oidx, self._vidx]) for (k, v) in self._adata_ref.layers.items()]
        else:
            pairs = self._layers.items()
        return [(k, v.copy()) for (k, v) in pairs] if copy else pairs

    def as_dict(self, copy=True):
        return dict(self.items(copy))

    def __len__(self):
        return len(self._adata_ref.layers if self.isview else self._layers)

    @property
    def isview(self):
        return self._adata_ref is not None
