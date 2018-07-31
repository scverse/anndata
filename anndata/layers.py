#thats just start
import numpy as np
from collections import OrderedDict

class AnnDataLayers():

    def __init__(self, adata, layers = None):

        self._adata = adata

        if layers is not None:
            for key in layers.keys():
                if layers[key].shape != self._adata.shape:
                    raise ValueError('Shape does not fit.')

        self._layers = OrderedDict(layers) if layers is not None else OrderedDict()

    def __getitem__(self, key):
        if key not in self._layers.keys():
            raise KeyError
        return self._layers[key]

    def __setitem__(self, key, value):
        if value.shape != self._adata.shape:
            raise ValueError('Shape does not fit.')
        self._layers[key] = value

    def keys(self):
        return self._layers.keys()
