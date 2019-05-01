from abc import ABC, abstractmethod
from collections.abc import MutableMapping
from typing import Mapping, Optional, Tuple

import numpy as np
import pandas as pd

if False:
    from .base import AnnData, Index  # noqa


class AlignedMapping(MutableMapping, ABC):
    """
    Attributes
    ----------
    axes: Tuple[int]
        Which axes is this aligned along?
    """

    @property
    @abstractmethod
    def isview(self):
        pass

    def _validate_value(self, val, key) -> None:
        """Raises an error if value is invalid"""
        for i, axis in enumerate(self.axes):
            if self.parent.shape[axis] != val.shape[i]:
                raise ValueError(
                    "Value passed for key '{key}' is of incorrect shape. Values of"
                    " {attrname} must match dimensions {axes} of parent."
                    "value had shape {wrong_shape}.".format(
                        key=key, attrname=self.attrname, axes=self.axes, wrong_shape=val.shape)
                )
        try: # TODO: Handle objects with indices
            # Could probably also re-order index if it's contained
            if not (val.index == self.dim_names).all():
                raise IndexError()  # Maybe not index error
        except AttributeError:
            pass

    @property
    @abstractmethod
    def axes(self) -> Tuple:
        pass

    @abstractmethod
    def view(self, parent, idx):
        """Returns a subetted copy-on-write view of the object"""
        pass

    @property
    def parent(self) -> "AnnData":
        return self._parent

    @property
    @abstractmethod
    def attrname(self) -> str:
        """What attr for the AnnData is this?"""
        pass
    
    @abstractmethod
    def copy(self):
        pass

    def __repr__(self):
        return f"{self.__class__.__name__} with keys: {', '.join(self.keys())}"

    def as_dict(self):
        return dict(self.items())


class AlignedViewMixin:
    """
    Attributes
    ----------
    parent
        Reference to parent AnnData view
    attrname
        What attribute in the parent is this?
    parent_mapping
        The object this is a view of.
    """
    isview = True

    def __getitem__(self, key):
        return self.parent_mapping[key][self.subset_idx]
    
    def __iter__(self):
        return self.parent_mapping.keys().__iter__()
    
    def __len__(self):
        return len(self.parent_mapping)

    def __delitem__(self, key):
        adata = self.parent.copy()
        getattr(adata, self.attrname).__delitem__(key)
        self.parent._init_as_actual(adata)

    def __setitem__(self, key, value):
        self._validate_value(value, key)  # Validate before mutating
        self.parent._init_as_actual(self.parent.copy())
        new_mapping = getattr(self.parent, self.attrname)
        new_mapping[key] = value

    def view(self, parent: "AnnData", subset_idx):
        """Returns a subsetted view of this object"""
        return self.__class__(self, parent, subset_idx)

class AlignedActualMixin:
    """
    Attributes
    _data: dict
        Underlying mapping to the data
    """
    isview = False

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._validate_value(value, key)
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]
    
    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return self._data.__len__()
    
    def __contains__(self, k):
        return self._data.__contains__(k)

class AxisArraysBase(AlignedMapping):
    """
    Mapping of key: array-like, where array-like is aligned to an axis of parent AnnData.
    """
    _dimnames = ("obs", "var")

    @property
    def axes(self):
        """Axes of the parent this is aligned to"""
        return (self._axis,)

    @property
    def dim(self):
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]

    @property
    def attrname(self):
        return "{}m".format(self.dim)

    def flipped(self):
        """transpose"""
        new = self.copy()
        new.dimension = abs(self._axis - 1)
        return new

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=self.dim_names)
        for key in self.keys():
            value = self[key]
            for icolumn, column in enumerate(value.T):
                df['{}{}'.format(key, icolumn + 1)] = column
        return df
    
    def copy(self):
        d = AxisArrays(self.parent, self._axis)
        for k, v in self.items():
            d[k] = v.copy()
        return d


class AxisArrays(AlignedActualMixin, AxisArraysBase):
    def __init__(self, parent, axis: int):
        self._parent = parent
        if axis not in (0, 1):
            raise ValueError()
        self._axis = axis
        self.dim_names = (parent.obs_names, parent.var_names)[self._axis]
        self._data = dict()

    def view(self, parent, subset):
        return AxisArraysView(self, parent, subset)


class AxisArraysView(AlignedViewMixin, AxisArraysBase):
    def __init__(self, parent_mapping, parent_view, subset_idx):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
        self._axis = parent_mapping._axis
        self.dim_names = parent_mapping.dim_names[subset_idx]


class LayersBase(AlignedMapping):
    """
    Mapping of key: array-like, where array-like is aligned to both axes of the parent anndata.
    """
    attrname = "layers"
    axes = (0, 1)

    def copy(self):
        return Layers(
            self.parent,
            vals={k: v.copy() for k, v in self.items()}
        )

class Layers(AlignedActualMixin, LayersBase):
    def __init__(self, parent: 'AnnData', vals: Optional[Mapping] = None):
        self._parent = parent
        self._data = {}
        if vals is not None:
            self._data.update(vals)

    def view(self, parent, subset_idx):
        return LayersView(self, parent, subset_idx)


class LayersView(AlignedViewMixin, LayersBase):
    def __init__(self, parent_mapping, parent_view, subset_idx):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
