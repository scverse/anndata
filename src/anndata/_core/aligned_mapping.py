from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import MutableMapping
from copy import copy
from dataclasses import dataclass, field
from itertools import chain
from typing import TYPE_CHECKING, cast, overload

import h5py
import numpy as np
import pandas as pd
import zarr
from scverse_misc import Deprecation, deprecated

from .._warnings import ExperimentalFeatureWarning, ImplicitModificationWarning
from ..abc import CSCDataset, CSRDataset
from ..compat import AwkArray, CupyArray, XDataArray, XDataset
from ..types import SupportsArrayApiBase
from ..utils import (
    asarray,
    axis_len,
    convert_to_dict,
    deprecation_msg,
    raise_value_error_if_multiindex_columns,
    warn,
    warn_once,
)
from .access import ElementRef
from .file_backing import to_memory
from .index import Idx1D, _subset
from .storage import _non_2d_message, coerce_array
from .views import as_view, view_update
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping
    from typing import ClassVar, Literal, Self, TypeAlias

    from ..compat import ZappyArray
    from ..typing import AlignedArray, InMemoryArray
    from .anndata import AnnData
    from .raw import Raw
    from .sparse_dataset import BaseCompressedSparseDataset

    type _AlignedAny = AlignedArray | pd.DataFrame


OneDIdx: TypeAlias = tuple[Idx1D]  # noqa: UP040
TwoDIdx: TypeAlias = tuple[Idx1D, Idx1D]  # noqa: UP040

OnDiskArray = h5py.Dataset | zarr.Array | CSRDataset | CSCDataset
"""Element types that live on disk and so cannot be copied in memory."""


def _copy_value[V: _AlignedAny](v: V) -> V:
    """Copy an array element in memory without copying on-disk elements."""
    # awkward arrays have immutable buffers, and on-disk elements aren’t copyable
    if isinstance(v, AwkArray | OnDiskArray | None):
        return copy(v)
    if isinstance(v, pd.DataFrame | Dataset2D | XDataArray):
        return cast("V", v.copy())
    return _copy_array(v)


@overload
def _copy_array[V: InMemoryArray](a: V) -> V: ...
@overload
def _copy_array(a: OnDiskArray | ZappyArray) -> InMemoryArray: ...
def _copy_array(
    a: InMemoryArray | OnDiskArray | ZappyArray,
) -> InMemoryArray:
    """An in-memory copy of an array element."""
    if isinstance(a, OnDiskArray):
        return to_memory(a)
    if isinstance(a, np.ndarray):
        return a.copy()
    if isinstance(a, SupportsArrayApiBase):
        # the array API standard has no `.copy()`, e.g. `torch.Tensor` lacks it
        return a.__array_namespace__().asarray(a, copy=True)
    return a.copy()


def _on_disk_x(
    elem: h5py.Group | h5py.Dataset | BaseCompressedSparseDataset,
) -> h5py.Dataset | CSRDataset | CSCDataset:
    """Interpret an on-disk `.X` element, wrapping sparse groups."""
    if isinstance(elem, h5py.Group):
        from .sparse_dataset import sparse_dataset

        return sparse_dataset(elem)
    if not isinstance(elem, h5py.Dataset):
        msg = f"Unexpected on-disk element of type {type(elem).__name__}"
        raise AssertionError(msg)
    return elem


class AlignedMappingBase[I: (OneDIdx, TwoDIdx), K: (str, str | None), V: _AlignedAny](
    MutableMapping[K, V], ABC
):
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
    """The parent object that this mapping is aligned to (`Raw` only for `.varm`)."""

    def __repr__(self) -> str:
        return f"{type(self).__name__} with keys: {', '.join(map(repr, self.keys()))}"

    def _ipython_key_completions_(self) -> list[K]:
        return list(self.keys())

    def _validate_value(self, val: _AlignedAny, key: K) -> V:
        """Raises an error if value is invalid"""
        if isinstance(val, AwkArray):
            msg = (
                "Support for Awkward Arrays is currently experimental. "
                "Behavior may change in the future. Please report any issues you may encounter!"
            )
            warn_once(msg, ExperimentalFeatureWarning)
        elif isinstance(val, np.ndarray | CupyArray) and len(val.shape) == 1:
            val = val.reshape((val.shape[0], 1))
        elif isinstance(val, XDataset):
            val = Dataset2D(val)
        own_axes: tuple[Literal[0, 1], ...] = (0, 1) if len(self.axes) == 2 else (0,)
        for i, axis in zip(own_axes, self.axes, strict=True):
            if self.parent.shape[axis] == axis_len(val, i):
                continue
            right_shape = tuple(self.parent.shape[a] for a in self.axes)
            actual_shape = tuple(axis_len(val, a) for a in own_axes)
            if actual_shape[i] is None and isinstance(val, AwkArray):
                dim = ("obs", "var")[i]
                msg = (
                    f"The AwkwardArray is of variable length in dimension {dim}. "
                    f"Try ak.to_regular(array, {i}) before including the array in AnnData"
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

    @property
    def _anndata_parent(self) -> AnnData:
        """The parent, which is only a `Raw` for `Raw.varm`."""
        from .anndata import AnnData

        if not isinstance(parent := self._parent, AnnData):
            msg = f"{self.attrname} is parented by {type(parent).__name__}, not AnnData"
            raise AssertionError(msg)
        return parent

    def copy(self) -> dict[K, V]:
        return {k: _copy_value(v) for k, v in self.items()}

    def _view(self, parent: AnnData, subset_idx: I) -> AlignedView[Self, I, K, V]:
        """Returns a subset copy-on-write view of the object."""
        return self._view_class(self, parent, subset_idx)

    @deprecated(Deprecation("0.10.2", deprecation_msg("as_dict", "dict(obj)")))
    def as_dict(self) -> dict:
        return dict(self)


class AlignedView[
    P: AlignedMappingBase,
    I: (OneDIdx, TwoDIdx),
    K: (str, str | None),
    V: _AlignedAny,
](AlignedMappingBase[I, K, V]):
    is_view: ClassVar[Literal[True]] = True

    # override docstring
    parent: AnnData
    """Reference to parent AnnData view"""

    parent_mapping: P
    """The object this is a view of."""

    subset_idx: I
    """The subset of the parent to view."""

    def __init__(self, parent_mapping: P, parent_view: AnnData, subset_idx: I) -> None:
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
        # LayersBase has no _axis, the rest does
        if (axis := getattr(parent_mapping, "_axis", None)) is not None:
            self._axis = axis

    def __getitem__(self, key: K) -> V:
        # Only `LayersBase()[None]` is ever `None`, so we override the return type there
        if self.parent_mapping[key] is None:
            return cast("V", None)
        return as_view(
            _subset(self.parent_mapping[key], self.subset_idx),
            ElementRef(self.parent, self.attrname, (key,)),
        )

    def __setitem__(self, key: K, value: _AlignedAny) -> None:
        value = self._validate_value(value, key)  # Validate before mutating
        msg = (
            f"Setting element `.{self.attrname}['{key}']` of view, "
            "initializing view as actual."
        )
        warn(msg, ImplicitModificationWarning)
        with view_update(self.parent, self.attrname, ()) as new_mapping:
            if value is None:
                del new_mapping[key]
            else:
                new_mapping[key] = value

    def __delitem__(self, key: K) -> None:
        if key not in self:
            msg = f"{key!r} not found in view of {self.attrname}"
            raise KeyError(msg)  # Make sure it exists before bothering with a copy
        if not (key is None and self.attrname == "layers"):
            msg = (
                f"Removing element `.{self.attrname}['{key}']` of view, "
                "initializing view as actual."
            )
            warn(msg, ImplicitModificationWarning)
        with view_update(self.parent, self.attrname, ()) as new_mapping:
            del new_mapping[key]

    def __contains__(self, key: object) -> bool:
        return key in self.parent_mapping

    def __iter__(self) -> Iterator[K]:
        return iter(self.parent_mapping)

    def __len__(self) -> int:
        return len(self.parent_mapping)


class AlignedActual[I: (OneDIdx, TwoDIdx), K: (str, str | None), V: _AlignedAny](
    AlignedMappingBase[I, K, V]
):
    is_view: ClassVar[Literal[False]] = False

    _data: MutableMapping[K, V]
    """Underlying mapping to the data"""

    def __init__(
        self,
        parent: AnnData | Raw,
        *,
        store: MutableMapping[K, V],
        validate: bool = True,
    ):
        self._parent = parent
        self._data = store
        for k, v in self._data.items():
            if v is None:
                continue
            self._data[k] = self._validate_value(v, k) if validate else v

    def __getitem__(self, key: K) -> V:
        return self._data[key]

    def __setitem__(self, key: K, value: _AlignedAny):
        if value is not None:
            value = self._validate_value(value, key)
        if key is None and value is None:
            del self[key]
        else:
            # `value` is only ever `None` here for `LayersBase()[None]`
            self._data[key] = cast("V", value)

    def __contains__(self, key: object) -> bool:
        return key in self._data

    def __delitem__(self, key: K):
        if key is None:
            self._data.pop(key, None)
        else:
            del self._data[key]

    def __iter__(self) -> Iterator[K]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)


class AxisArraysBase[V: _AlignedAny](AlignedMappingBase[OneDIdx, str, V]):
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
            for icolumn, column in enumerate(asarray(self[key]).T):
                df[f"{key}{icolumn + 1}"] = column
        return df

    def _validate_value(self, val: _AlignedAny, key: str) -> V:
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


class AxisArrays(AlignedActual[OneDIdx, str, "_AlignedAny"], AxisArraysBase):
    def __init__(
        self,
        parent: AnnData,
        *,
        axis: Literal[0, 1],
        store: MutableMapping[str, _AlignedAny] | AxisArraysBase,
        validate: bool = True,
    ):
        if axis not in {0, 1}:
            raise ValueError()
        self._axis = axis
        super().__init__(parent, store=store, validate=validate)


class AxisArraysView(
    AlignedView[AxisArraysBase, OneDIdx, str, "_AlignedAny"],
    AxisArraysBase["_AlignedAny"],
):
    pass


AxisArraysBase._view_class = AxisArraysView
AxisArraysBase._actual_class = AxisArrays


class LayersBase(AlignedMappingBase[TwoDIdx, str | None, "AlignedArray"]):
    """\
    Mapping of key: array-like, where array-like is aligned to both axes of the
    parent anndata.
    """

    _allow_df: ClassVar = False
    attrname: ClassVar[Literal["layers"]] = "layers"
    axes: ClassVar[tuple[Literal[0], Literal[1]]] = (0, 1)

    isbacked: bool

    def __bool__(self) -> bool:
        return not self.keys() <= {None}

    def clear(self, *, keep_x: bool | None = None) -> None:
        """Remove all layers.

        Parameters
        ----------
        keep_x
            Whether to keep `.X` (stored under the `None` key).
            The default (`None`) keeps ``.X`` and warns that a future
            release may drop it by default - use this parameter explicitly to silence the warning.
        """
        if keep_x is None:
            warn(
                "`.layers.clear()` currently keeps `.X` (stored under the "
                "`None` key), but a future release may drop it. Pass "
                "`keep_x=True` to keep `.X` or `keep_x=False` to drop it "
                "and silence this warning.",
                FutureWarning,
            )
            keep_x = True
        for key in [k for k in self if not (keep_x and k is None)]:
            del self[key]

    def _validate_value(
        self, val: _AlignedAny | None, key: str | None
    ) -> AlignedArray | None:
        """Warn if storing ``val`` under ``key`` would violate the on-disk spec.

        Called from the explicit write paths (``__setitem__`` and the
        :class:`AlignedMappingProperty` full-reassign hook) only; we do
        *not* warn from :meth:`_validate_value` because that runs on
        every property access via :meth:`AlignedActual.__init__`.

        The ``key=None`` slot mirrors ``AnnData.X`` and is reported by
        the ``X`` setter itself (with the better name "X").
        """
        if key is not None and (msg := _non_2d_message(val, name=f"Layer {key!r}")):
            warn(msg, UserWarning)
        return super()._validate_value(val, key)


class Layers(AlignedActual[TwoDIdx, str | None, "AlignedArray"], LayersBase):
    def __init__(
        self,
        parent: AnnData,
        *,
        store: MutableMapping[str | None, AlignedArray | None],
        validate: bool = True,
    ):
        super().__init__(parent, store=store, validate=validate)
        self.isbacked = (
            None not in self._data and self._anndata_parent.filename is not None
        )

    def __getitem__(self, key: str | None) -> AlignedArray | None:
        if key is None and self.isbacked:
            parent = self._anndata_parent
            if not parent.file.is_open:
                parent.file.open()
            return _on_disk_x(parent.file["X"])
        return super().__getitem__(key)

    def __iter__(self) -> Iterator[str | None]:
        keys_iter = super().__iter__()
        if self.isbacked:
            yield from chain([None], keys_iter)
        yield from keys_iter

    def __len__(self) -> int:
        data_length = super().__len__()
        if self.isbacked:
            return data_length + 1
        return data_length

    def __contains__(self, key: object) -> bool:
        if key is None and self.isbacked:
            return True
        return super().__contains__(key)


class LayersView(
    AlignedView[LayersBase, TwoDIdx, str | None, "AlignedArray"], LayersBase
):
    def __init__(
        self, parent_mapping: LayersBase, parent_view: AnnData, subset_idx: TwoDIdx
    ) -> None:
        super().__init__(parent_mapping, parent_view, subset_idx)
        self.isbacked = parent_mapping.isbacked

    def __getitem__(self, key: str | None) -> AlignedArray | None:
        return super().__getitem__(key)


LayersBase._view_class = LayersView
LayersBase._actual_class = Layers


# OneDIdx because this is aligned to one axis
class PairwiseArraysBase(AlignedMappingBase[OneDIdx, str, "AlignedArray"]):
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
        return (0, 0) if self._axis == 0 else (1, 1)

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]


class PairwiseArrays(AlignedActual[OneDIdx, str, "AlignedArray"], PairwiseArraysBase):
    def __init__(
        self,
        parent: AnnData,
        *,
        axis: Literal[0, 1],
        store: MutableMapping[str, AlignedArray],
        validate: bool = True,
    ):
        if axis not in {0, 1}:
            raise ValueError()
        self._axis = axis
        super().__init__(parent, store=store, validate=validate)


class PairwiseArraysView(
    AlignedView[PairwiseArraysBase, OneDIdx, str, "AlignedArray"], PairwiseArraysBase
):
    pass


PairwiseArraysBase._view_class = PairwiseArraysView
PairwiseArraysBase._actual_class = PairwiseArrays


type AlignedMapping = (
    AxisArrays
    | AxisArraysView
    | Layers
    | LayersView
    | PairwiseArrays
    | PairwiseArraysView
)
"""Pair of types to be aligned."""


@dataclass
class AlignedMappingProperty[A: AlignedActual, V: AlignedView, K: (str, str | None)](
    property
):
    """A :class:`property` that creates an ephemeral AlignedMapping.

    The actual data is stored as `f'_{self.name}'` in the parent object.
    """

    cls: type[A]
    """Concrete type that will be constructed."""
    axis: Literal[0, 1] | None = None
    """Axis of the parent to align to."""
    name: str = field(init=False)
    """Name of the attribute in the parent object, set by `__set_name__`."""

    def construct(
        self,
        obj: AnnData | Raw,
        *,
        store: MutableMapping[K, AlignedArray],
        validate: bool = True,
    ) -> A:
        # only the 1D mappings (`AxisArrays`/`PairwiseArrays`) take an `axis`
        axis_kw: dict[str, Literal[0, 1]] = (
            {} if self.axis is None else {"axis": self.axis}
        )
        return self.cls(obj, store=store, validate=validate, **axis_kw)

    def __post_init__(self) -> None:
        """Install a fake `fget` for sphinx-autodoc-typehints."""

        def fake(obj): ...

        fake.__annotations__ = {"return": self.cls._actual_class | self.cls._view_class}
        property.__init__(self, fake)

    def __set_name__(self, owner: AnnData, name: str):
        self.name = name

    @overload
    def __get__(self, obj: None, objtype: type | None = None, /) -> Self: ...
    @overload
    def __get__(self, obj: AnnData | Raw, objtype: type | None = None, /) -> A | V: ...
    def __get__(
        self, obj: AnnData | Raw | None, objtype: type | None = None, /
    ) -> Self | A | V:
        if obj is None:
            # When accessed from the class, e.g. via `AnnData.obs`,
            # this needs to return a `property` instance, e.g. for Sphinx
            return self
        from .anndata import AnnData

        # only `AnnData` can be a view
        if not isinstance(obj, AnnData) or (parent_anndata := obj._adata_ref) is None:
            # all stores are created through `.__set__`, so no need to double-validate.
            return self.construct(
                obj, store=getattr(obj, f"_{self.name}"), validate=False
            )
        idxs = (obj._oidx, obj._vidx)
        parent: A = getattr(parent_anndata, self.name)
        axes = parent.axes
        subset_idx = (
            (idxs[axes[0]],) if len(axes) == 1 else (idxs[axes[0]], idxs[axes[1]])
        )
        return cast("V", parent._view(obj, subset_idx))

    def __set__(
        self,
        obj: AnnData | Raw,
        value: Mapping[K, AlignedArray] | Iterable[tuple[K, AlignedArray]] | None,
    ) -> None:
        from .anndata import AnnData

        value = convert_to_dict(value)
        _ = self.construct(obj, store=value)  # Validate and convert arrays in `value`
        if isinstance(obj, AnnData) and obj.is_view:  # only `AnnData` can be a view
            obj._init_as_actual(obj.copy())
        prev = getattr(obj, f"_{self.name}", None)
        if issubclass(self.cls, LayersBase):
            if None in value:
                if value[None] is None:
                    value = {k: v for k, v in value.items() if k is not None}
            elif prev is not None and (x := prev.get(None)) is not None:
                warn(
                    "Assigning to `.layers` without a `None` key currently "
                    "preserves `.X` (stored under the `None` key), but a "
                    "future release may drop it. Include `.X` explicitly as "
                    "`{None: adata.X, ...}` to keep it, or `{None: None, ...}` "
                    "to drop it and silence this warning.",
                    FutureWarning,
                )
                value = {**value, None: x}
        setattr(obj, f"_{self.name}", value)

    def __delete__(self, obj: AnnData) -> None:
        if issubclass(self.cls, LayersBase) and (
            (x := getattr(obj, self.name).get(None)) is not None
        ):
            warn(
                "`del adata.layers` currently keeps `.X` (stored under the "
                "`None` key), but a future release may drop it. Use "
                "`adata.layers.clear(keep_x=True)` to keep `.X` or "
                "`adata.layers.clear(keep_x=False)` to drop it and silence "
                "this warning.",
                FutureWarning,
            )
            # Reassign with `.X` prevents another, less-helpful warning
            setattr(obj, self.name, {None: x})
            return
        setattr(obj, self.name, {})
