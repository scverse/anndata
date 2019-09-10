from abc import ABC, abstractmethod
import collections.abc as cabc
from functools import singledispatch
from typing import (
    Mapping,
    Optional,
    Tuple,
    Iterable,
    Sequence,
    Union,
    List,
    TypeVar,
    Iterator,
    ClassVar,
)

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix

from ..utils import deprecated
from .views import asview, ViewArgs


OneDIdx = Sequence[Union[int, bool]]
TwoDIdx = Tuple[OneDIdx, OneDIdx]

I = TypeVar("I", OneDIdx, TwoDIdx)
V = TypeVar("V", pd.DataFrame, spmatrix, np.ndarray)


@singledispatch
def _subset(a: V, subset_idx: I):
    # Select as combination of indexes, not coordinates
    # Correcting for indexing behaviour of np.ndarray
    if all(isinstance(x, Iterable) for x in subset_idx):
        subset_idx = np.ix_(*subset_idx)
    return a[subset_idx]


@_subset.register(pd.DataFrame)
def _subset_df(df: pd.DataFrame, subset_idx: OneDIdx):
    return df.iloc[subset_idx]


class AlignedMapping(cabc.MutableMapping, ABC):
    """
    Attributes
    ----------
    axes
        Which axes is this aligned along?
    _view_class
        The view class for this aligned mapping.
    _actual_class
        The actual class (which has it's own data) for this aligned mapping.
    """

    _view_class: ClassVar["AlignedViewMixin"]
    _actual_class: ClassVar["AlignedActualMixin"]

    def __repr__(self):
        return f"{self.__class__.__name__,} with keys: {self.keys()}"

    def _ipython_key_completions_(self) -> List[str]:
        return list(self.keys())

    def _validate_value(self, val: V, key: str) -> None:
        """Raises an error if value is invalid"""
        for i, axis in enumerate(self.axes):
            if self.parent.shape[axis] != val.shape[i]:
                right_shape = tuple(self.parent.shape[a] for a in self.axes)
                raise ValueError(
                    f"Value passed for key '{key}' is of incorrect shape. "
                    f"Values of {self.attrname} must match dimensions "
                    f"{self.axes} of parent. Value had shape {val.shape} while "
                    f"it should have had {right_shape}."
                )
        try:  # TODO: Handle objects with indices
            # Could probably also re-order index if it's contained
            if not (val.index == self.dim_names).all():
                raise IndexError()  # Maybe not index error
        except AttributeError:
            pass

    @property
    @abstractmethod
    def attrname(self) -> str:
        """What attr for the AnnData is this?"""
        pass

    @property
    @abstractmethod
    def axes(self) -> Tuple[int, ...]:
        """Which axes of the parent is this aligned to?"""
        pass

    @property
    @abstractmethod
    def isview(self) -> bool:
        pass

    @property
    def parent(self) -> "AnnData":
        return self._parent

    def copy(self):
        d = self._actual_class(self.parent, self._axis)
        for k, v in self.items():
            d[k] = v.copy()
        return d

    def _view(self, parent: "AnnData", subset_idx: I):
        """Returns a subset copy-on-write view of the object."""
        return self._view_class(self, parent, subset_idx)

    @deprecated("dict(obj)")
    def as_dict(self) -> dict:
        return dict(self)


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

    def __getitem__(self, key: str) -> V:
        return asview(
            _subset(self.parent_mapping[key], self.subset_idx),
            ViewArgs(self.parent, self.attrname, (key,)),
        )

    def __setitem__(self, key: str, value: V):
        self._validate_value(value, key)  # Validate before mutating
        adata = self.parent.copy()
        new_mapping = getattr(adata, self.attrname)
        new_mapping[key] = value
        self.parent._init_as_actual(adata)

    def __delitem__(self, key: str):
        self[key]  # Make sure it exists before bothering with a copy
        adata = self.parent.copy()
        new_mapping = getattr(adata, self.attrname)
        del new_mapping[key]
        self.parent._init_as_actual(adata)

    def __contains__(self, key: str) -> bool:
        return key in self.parent_mapping

    def __iter__(self) -> Iterator[str]:
        return iter(self.parent_mapping)

    def __len__(self) -> int:
        return len(self.parent_mapping)


class AlignedActualMixin:
    """
    Attributes
    _data: dict
        Underlying mapping to the data
    """

    isview = False

    def __getitem__(self, key: str) -> V:
        return self._data[key]

    def __setitem__(self, key: str, value: V):
        self._validate_value(value, key)
        self._data[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self._data

    def __delitem__(self, key: str):
        del self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)


class AxisArraysBase(AlignedMapping):
    """
    Mapping of key→array-like,
    where array-like is aligned to an axis of parent AnnData.
    """

    _dimnames = ("obs", "var")

    @property
    def attrname(self):
        return f"{self.dim}m"

    @property
    def axes(self) -> Tuple[int]:
        """Axes of the parent this is aligned to"""
        return (self._axis,)

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]

    def flipped(self) -> "AxisArraysBase":
        """Transpose."""
        new = self.copy()
        new.dimension = abs(self._axis - 1)
        return new

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=self.dim_names)
        for key in self.keys():
            value = self[key]
            for icolumn, column in enumerate(value.T):
                df[f'{key}{icolumn + 1}'] = column
        return df


class AxisArrays(AlignedActualMixin, AxisArraysBase):
    def __init__(
        self, parent: "AnnData", axis: int, vals: Optional[Mapping] = None
    ):
        self._parent = parent
        if axis not in (0, 1):
            raise ValueError()
        self._axis = axis
        self.dim_names = (parent.obs_names, parent.var_names)[self._axis]
        self._data = dict()
        if vals is not None:
            self.update(vals)


class AxisArraysView(AlignedViewMixin, AxisArraysBase):
    def __init__(
        self,
        parent_mapping: AxisArraysBase,
        parent_view: "AnnData",
        subset_idx: OneDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
        self._axis = parent_mapping._axis
        self.dim_names = parent_mapping.dim_names[subset_idx]


AxisArraysBase._view_class = AxisArraysView
AxisArraysBase._actual_class = AxisArrays


class LayersBase(AlignedMapping):
    """
    Mapping of key: array-like, where array-like is aligned to both axes of the
    parent anndata.
    """

    attrname = "layers"
    axes = (0, 1)

    # TODO: I thought I had a more elegant solution to overiding this...
    def copy(self) -> "Layers":
        d = self._actual_class(self.parent)
        for k, v in self.items():
            d[k] = v.copy()
        return d


class Layers(AlignedActualMixin, LayersBase):
    def __init__(self, parent: "AnnData", vals: Optional[Mapping] = None):
        self._parent = parent
        self._data = dict()
        if vals is not None:
            self.update(vals)


class LayersView(AlignedViewMixin, LayersBase):
    def __init__(
        self,
        parent_mapping: LayersBase,
        parent_view: "AnnData",
        subset_idx: TwoDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx


LayersBase._view_class = LayersView
LayersBase._actual_class = Layers


class PairwiseArraysBase(AlignedMapping):
    """
    Mapping of key: array-like, where both axes of array-like are aligned to
    one axis of the parent anndata.
    """

    _dimnames = ("obs", "var")

    @property
    def attrname(self) -> str:
        return f"{self.dim}p"

    @property
    def axes(self) -> Tuple[int, int]:
        """Axes of the parent this is aligned to"""
        return (self._axis, self._axis)

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]


class PairwiseArrays(AlignedActualMixin, PairwiseArraysBase):
    def __init__(
        self, parent: "AnnData", axis: int, vals: Optional[Mapping] = None
    ):
        self._parent = parent
        if axis not in (0, 1):
            raise ValueError()
        self._axis = axis
        self.dim_names = (parent.obs_names, parent.var_names)[self._axis]
        self._data = dict()
        if vals is not None:
            self.update(vals)


class PairwiseArraysView(AlignedViewMixin, PairwiseArraysBase):
    def __init__(
        self,
        parent_mapping: PairwiseArraysBase,
        parent_view: "AnnData",
        subset_idx: OneDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = (subset_idx, subset_idx)
        self._axis = parent_mapping._axis
        self.dim_names = parent_mapping.dim_names[subset_idx]


PairwiseArraysBase._view_class = PairwiseArraysView
PairwiseArraysBase._actual_class = PairwiseArrays
