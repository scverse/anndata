from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from collections.abc import MutableMapping, Sequence
from copy import copy
from dataclasses import dataclass
from typing import TYPE_CHECKING, Generic, TypeVar

import numpy as np
import pandas as pd

from .._warnings import ExperimentalFeatureWarning, ImplicitModificationWarning
from ..compat import AwkArray, CSArray, CSMatrix, CupyArray, XDataset
from ..utils import (
    axis_len,
    convert_to_dict,
    deprecated,
    raise_value_error_if_multiindex_columns,
    warn_once,
)
from .access import ElementRef
from .index import _subset
from .storage import coerce_array
from .views import as_view, view_update
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Mapping
    from typing import ClassVar, Literal, Self

    from .anndata import AnnData
    from .raw import Raw


OneDIdx = Sequence[int] | Sequence[bool] | slice
TwoDIdx = tuple[OneDIdx, OneDIdx]
# TODO: pd.DataFrame only allowed in AxisArrays?
Value = pd.DataFrame | CSMatrix | CSArray | np.ndarray

P = TypeVar("P", bound="AlignedMappingBase")
"""Parent mapping an AlignedView is based on."""
I = TypeVar("I", OneDIdx, TwoDIdx)


class AlignedMappingBase(MutableMapping[str, Value], ABC):
    """\
    An abstract base class for Mappings containing array-like values aligned
    to either one or both AnnData axes.
    """

    _allow_df: ClassVar[bool]
    """If this mapping supports heterogeneous DataFrames"""

    _view_class: ClassVar[type[AlignedView]]
    """The view class for this aligned mapping."""

    _actual_class: ClassVar[type[AlignedActual]]
    """The actual class (which has it’s own data) for this aligned mapping."""

    _parent: AnnData | Raw
    """The parent object that this mapping is aligned to."""

    def __repr__(self):
        return f"{type(self).__name__} with keys: {', '.join(self.keys())}"

    def _ipython_key_completions_(self) -> list[str]:
        return list(self.keys())

    def _validate_value(self, val: Value, key: str) -> Value:
        """Raises an error if value is invalid"""
        if isinstance(val, AwkArray):
            warn_once(
                "Support for Awkward Arrays is currently experimental. "
                "Behavior may change in the future. Please report any issues you may encounter!",
                ExperimentalFeatureWarning,
                # stacklevel=3,
            )
        elif isinstance(val, np.ndarray | CupyArray) and len(val.shape) == 1:
            val = val.reshape((val.shape[0], 1))
        elif isinstance(val, XDataset):
            val = Dataset2D(data_vars=val.data_vars, coords=val.coords, attrs=val.attrs)
        for i, axis in enumerate(self.axes):
            if self.parent.shape[axis] == axis_len(val, i):
                continue
            right_shape = tuple(self.parent.shape[a] for a in self.axes)
            actual_shape = tuple(axis_len(val, a) for a, _ in enumerate(self.axes))
            if actual_shape[i] is None and isinstance(val, AwkArray):
                dim = ("obs", "var")[i]
                msg = (
                    f"The AwkwardArray is of variable length in dimension {dim}.",
                    f"Try ak.to_regular(array, {i}) before including the array in AnnData",
                )
            else:
                dims = tuple(("obs", "var")[ax] for ax in self.axes)
                msg = (
                    f"Value passed for key {key!r} is of incorrect shape. "
                    f"Values of {self.attrname} must match dimensions {dims} of parent. "
                    f"Value had shape {actual_shape} while it should have had {right_shape}."
                )
            raise ValueError(msg)
        name = f"{self.attrname.title().rstrip('s')} {key!r}"
        return coerce_array(val, name=name, allow_df=self._allow_df)

    @property
    @abstractmethod
    def attrname(self) -> str:
        """What attr for the AnnData is this?"""

    @property
    @abstractmethod
    def axes(self) -> tuple[Literal[0, 1], ...]:
        """Which axes of the parent is this aligned to?"""

    @property
    @abstractmethod
    def is_view(self) -> bool: ...

    @property
    def parent(self) -> AnnData | Raw:
        return self._parent

    def copy(self) -> dict[str, Value]:
        # Shallow copy for awkward array since their buffers are immutable
        return {
            k: copy(v) if isinstance(v, AwkArray) else v.copy() for k, v in self.items()
        }

    def _view(self, parent: AnnData, subset_idx: I) -> AlignedView[Self, I]:
        """Returns a subset copy-on-write view of the object."""
        return self._view_class(self, parent, subset_idx)

    @deprecated("dict(obj)")
    def as_dict(self) -> dict:
        return dict(self)


class AlignedView(AlignedMappingBase, Generic[P, I]):
    is_view: ClassVar[Literal[True]] = True

    # override docstring
    parent: AnnData
    """Reference to parent AnnData view"""

    attrname: str
    """What attribute in the parent is this?"""

    parent_mapping: P
    """The object this is a view of."""

    subset_idx: I
    """The subset of the parent to view."""

    def __init__(self, parent_mapping: P, parent_view: AnnData, subset_idx: I):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
        if hasattr(parent_mapping, "_axis"):
            # LayersBase has no _axis, the rest does
            self._axis = parent_mapping._axis  # type: ignore

    def __getitem__(self, key: str) -> Value:
        return as_view(
            _subset(self.parent_mapping[key], self.subset_idx),
            ElementRef(self.parent, self.attrname, (key,)),
        )

    def __setitem__(self, key: str, value: Value) -> None:
        value = self._validate_value(value, key)  # Validate before mutating
        warnings.warn(
            f"Setting element `.{self.attrname}['{key}']` of view, "
            "initializing view as actual.",
            ImplicitModificationWarning,
            stacklevel=2,
        )
        with view_update(self.parent, self.attrname, ()) as new_mapping:
            new_mapping[key] = value

    def __delitem__(self, key: str) -> None:
        if key not in self:
            msg = f"{key!r} not found in view of {self.attrname}"
            raise KeyError(msg)  # Make sure it exists before bothering with a copy
        warnings.warn(
            f"Removing element `.{self.attrname}['{key}']` of view, "
            "initializing view as actual.",
            ImplicitModificationWarning,
            stacklevel=2,
        )
        with view_update(self.parent, self.attrname, ()) as new_mapping:
            del new_mapping[key]

    def __contains__(self, key: str) -> bool:
        return key in self.parent_mapping

    def __iter__(self) -> Iterator[str]:
        return iter(self.parent_mapping)

    def __len__(self) -> int:
        return len(self.parent_mapping)


class AlignedActual(AlignedMappingBase):
    is_view: ClassVar[Literal[False]] = False

    _data: MutableMapping[str, Value]
    """Underlying mapping to the data"""

    def __init__(self, parent: AnnData | Raw, *, store: MutableMapping[str, Value]):
        self._parent = parent
        self._data = store
        for k, v in self._data.items():
            self._data[k] = self._validate_value(v, k)

    def __getitem__(self, key: str) -> Value:
        return self._data[key]

    def __setitem__(self, key: str, value: Value):
        value = self._validate_value(value, key)
        self._data[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self._data

    def __delitem__(self, key: str):
        del self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)


class AxisArraysBase(AlignedMappingBase):
    """\
    Mapping of key→array-like,
    where array-like is aligned to an axis of parent AnnData.
    """

    _allow_df: ClassVar = True
    _dimnames: ClassVar = ("obs", "var")

    _axis: Literal[0, 1]

    @property
    def attrname(self) -> str:
        return f"{self.dim}m"

    @property
    def axes(self) -> tuple[Literal[0, 1]]:
        """Axes of the parent this is aligned to"""
        return (self._axis,)

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=self.dim_names)
        for key in self.keys():
            value = self[key]
            for icolumn, column in enumerate(value.T):
                df[f"{key}{icolumn + 1}"] = column
        return df

    def _validate_value(self, val: Value, key: str) -> Value:
        if isinstance(val, pd.DataFrame):
            raise_value_error_if_multiindex_columns(val, f"{self.attrname}[{key!r}]")
            if not val.index.equals(self.dim_names):
                # Could probably also re-order index if it’s contained
                try:
                    pd.testing.assert_index_equal(val.index, self.dim_names)
                except AssertionError as e:
                    msg = f"value.index does not match parent’s {self.dim} names:\n{e}"
                    raise ValueError(msg) from None
                else:
                    msg = "Index.equals and pd.testing.assert_index_equal disagree"
                    raise AssertionError(msg)
            val.index.name = (
                self.dim_names.name
            )  # this is consistent with AnnData.obsm.setter and AnnData.varm.setter
        return super()._validate_value(val, key)

    @property
    def dim_names(self) -> pd.Index:
        return (self.parent.obs_names, self.parent.var_names)[self._axis]


class AxisArrays(AlignedActual, AxisArraysBase):
    def __init__(
        self,
        parent: AnnData | Raw,
        *,
        axis: Literal[0, 1],
        store: MutableMapping[str, Value] | AxisArraysBase,
    ):
        if axis not in {0, 1}:
            raise ValueError()
        self._axis = axis
        super().__init__(parent, store=store)


class AxisArraysView(AlignedView[AxisArraysBase, OneDIdx], AxisArraysBase):
    pass


AxisArraysBase._view_class = AxisArraysView
AxisArraysBase._actual_class = AxisArrays


class LayersBase(AlignedMappingBase):
    """\
    Mapping of key: array-like, where array-like is aligned to both axes of the
    parent anndata.
    """

    _allow_df: ClassVar = False
    attrname: ClassVar[Literal["layers"]] = "layers"
    axes: ClassVar[tuple[Literal[0], Literal[1]]] = (0, 1)


class Layers(AlignedActual, LayersBase):
    pass


class LayersView(AlignedView[LayersBase, TwoDIdx], LayersBase):
    pass


LayersBase._view_class = LayersView
LayersBase._actual_class = Layers


class PairwiseArraysBase(AlignedMappingBase):
    """\
    Mapping of key: array-like, where both axes of array-like are aligned to
    one axis of the parent anndata.
    """

    _allow_df: ClassVar = False
    _dimnames: ClassVar = ("obs", "var")

    _axis: Literal[0, 1]

    @property
    def attrname(self) -> str:
        return f"{self.dim}p"

    @property
    def axes(self) -> tuple[Literal[0], Literal[0]] | tuple[Literal[1], Literal[1]]:
        """Axes of the parent this is aligned to"""
        return self._axis, self._axis  # type: ignore

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]


class PairwiseArrays(AlignedActual, PairwiseArraysBase):
    def __init__(
        self,
        parent: AnnData,
        *,
        axis: Literal[0, 1],
        store: MutableMapping[str, Value],
    ):
        if axis not in {0, 1}:
            raise ValueError()
        self._axis = axis
        super().__init__(parent, store=store)


class PairwiseArraysView(AlignedView[PairwiseArraysBase, OneDIdx], PairwiseArraysBase):
    pass


PairwiseArraysBase._view_class = PairwiseArraysView
PairwiseArraysBase._actual_class = PairwiseArrays


AlignedMapping = (
    AxisArrays
    | AxisArraysView
    | Layers
    | LayersView
    | PairwiseArrays
    | PairwiseArraysView
)
T = TypeVar("T", bound=AlignedMapping)
"""Pair of types to be aligned."""


@dataclass
class AlignedMappingProperty(property, Generic[T]):
    """A :class:`property` that creates an ephemeral AlignedMapping.

    The actual data is stored as `f'_{self.name}'` in the parent object.
    """

    name: str
    """Name of the attribute in the parent object."""
    cls: type[T]
    """Concrete type that will be constructed."""
    axis: Literal[0, 1] | None = None
    """Axis of the parent to align to."""

    def construct(self, obj: AnnData, *, store: MutableMapping[str, Value]) -> T:
        if self.axis is None:
            return self.cls(obj, store=store)
        return self.cls(obj, axis=self.axis, store=store)

    @property
    def fget(self) -> Callable[[], None]:
        """Fake fget for sphinx-autodoc-typehints."""

        def fake(): ...

        fake.__annotations__ = {"return": self.cls._actual_class | self.cls._view_class}
        return fake

    def __get__(self, obj: None | AnnData, objtype: type | None = None) -> T:
        if obj is None:
            # When accessed from the class, e.g. via `AnnData.obs`,
            # this needs to return a `property` instance, e.g. for Sphinx
            return self  # type: ignore
        if not obj.is_view:
            return self.construct(obj, store=getattr(obj, f"_{self.name}"))
        parent_anndata = obj._adata_ref
        idxs = (obj._oidx, obj._vidx)
        parent: AlignedMapping = getattr(parent_anndata, self.name)
        return parent._view(obj, tuple(idxs[ax] for ax in parent.axes))

    def __set__(
        self, obj: AnnData, value: Mapping[str, Value] | Iterable[tuple[str, Value]]
    ) -> None:
        value = convert_to_dict(value)
        _ = self.construct(obj, store=value)  # Validate
        if obj.is_view:
            obj._init_as_actual(obj.copy())
        setattr(obj, f"_{self.name}", value)

    def __delete__(self, obj) -> None:
        setattr(obj, self.name, dict())
