from abc import ABC, abstractmethod
from collections import OrderedDict
from collections.abc import MutableMapping
from typing import Mapping, Optional, Union, Tuple

import numpy as np
from scipy.sparse import issparse

if False:
    from .base import AnnData, Index  # noqa


class AlignedMapping(MutableMapping, ABC):

    @property
    @abstractmethod
    def isview(self):
        pass

    @abstractmethod # This could probably just get implemented here
    def _validate_value(self, val, key) -> "NoneType":
        """Raises an error if value is invalid"""
        pass

    # TODO
    # I need to figure out what these idx value are
    # def _validate_value(self, val, key):
    #     for curidx, vallen in zip(self.idx, value.shape):
    #         if curlen == slice(None):
    #             continue
    #         curl


    @abstractmethod
    def view(self, parent, idx):
        """Return a view of the object"""
        pass

    @property
    @abstractmethod
    def parent(self) -> "AnnData":
        """The AnnData this is connected to"""
        pass

    @property
    @abstractmethod
    def idx(self) -> Tuple:
        """How to slice it"""
        pass

    @property
    @abstractmethod
    def attrname(self) -> str:
        """What attr for the AnnData is this?"""
        pass

    def __repr__(self):
        return f"{self.__class__.__name__}: with keys: {', '.join(self._data.keys())}"

    def as_dict(self):
        return dict(self.items())

class AlignedViewMixin:
    @property
    def isview(self):
        return True

    @property
    def parent_mapping(self):
        return getattr(self.parent._adata_ref, self.attrname)

    def __getitem__(self, key):
        return self.parent_mapping[key][self.idx]

    def __delitem__(self, key):
        adata = self.parent.copy()
        getattr(adata, self.attrname).__delitem__(key)
        self.parent._init_as_actual(adata)

    def __setitem__(self, key, value):
        adata = self.parent.copy()
        getattr(adata, self.attrname).__setitem__(key, value)
        self.parent._init_as_actual(adata)

    def view(self, parent: "AnnData", idx):
        """Returns a subsetted view of this object"""
        return self.__class__(self, parent, idx)

class AlignedActualMixin:
    @property
    def isview(self):
        return False

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._validate_value(value, key)
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

class LayersBase(AlignedMapping):
    @property
    def attrname(self):
        return "layers"

class Layers(LayersBase, AlignedActualMixin):
    pass

class LayersView(LayersBase, AlignedViewMixin):
    pass

class AxisArraysBase(AlignedMapping):
    _dimnames = ("obs", "var")

    @property
    def dim(self):
        return self._dimnames[self._axis]

    @property
    def attrname(self):
        return "{}m".format(self.dim)

class AxisArrays(AxisArraysBase, AlignedActualMixin):
    pass

class AxisArraysView(AxisArraysBase, AlignedViewMixin):
    pass

# class LayersBase(AlignedMapping):
#     # def __init__(
#         # self,
#         # adata
#     # )
#     # def view(self, subset):


#     def __iter__(self):
#         return self._data.__iter__()

#     def __len__(self):
#         return self._data.__len__()


# class LayersView(LayersBase):
#     def __init__(self, parent, oidx, vidx):
#         self.parent = parent
#         self.parent_layer = parent_layer
#         self._oidx = oidx
#         self._vidx = vidx

#     def __getitem__(self, key):
#         return self._adata_ref.layers[key][self._oidx, self._vidx]

# class Layers(LayersBase):
#     def __init__(
#         self,
#         parent: "AnnData",
#         layers: Mapping = {}
#     ):
#         self.parent = parent
#         self.obs_names = parent.var_names
#         self.var_names = parent.obs_names
#         self._data = {}
#         self.update(layers)

#     def __getitem__(self, key):
#         return self._data[key]

#     def __setitem__(self, key, value):
#         if value.shape != self.parent.shape:
#             raise ValueError()
#         self._data[key] = value

#     def __delitem__(self, key):
#         del self._data[key]


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

    def __contains__(self, key):
        return key in self._layers

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
