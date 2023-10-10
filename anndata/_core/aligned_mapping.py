from abc import ABC, abstractmethod
from collections import abc as cabc
from copy import copy
from typing import Union, Optional, Type, ClassVar, TypeVar  # Special types
from typing import Iterator, Mapping, Sequence  # ABCs
from typing import Tuple, List, Dict  # Generic base types
import warnings

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix

from ..utils import deprecated, ensure_df_homogeneous, dim_len
from . import raw, anndata
from .views import as_view, view_update
from .access import ElementRef
from .index import _subset
from anndata.compat import AwkArray
from anndata._warnings import ExperimentalFeatureWarning, ImplicitModificationWarning


OneDIdx = Union[Sequence[int], Sequence[bool], slice]
TwoDIdx = Tuple[OneDIdx, OneDIdx]

I = TypeVar("I", OneDIdx, TwoDIdx, covariant=True)
# TODO: pd.DataFrame only allowed in AxisArrays?
V = Union[pd.DataFrame, spmatrix, np.ndarray]


class AlignedMapping(cabc.MutableMapping, ABC):
    """\
    An abstract base class for Mappings containing array-like values aligned
    to either one or both AnnData axes.
    """

    _allow_df: ClassVar[bool]
    """If this mapping supports heterogeneous DataFrames"""

    _view_class: ClassVar[Type["AlignedViewMixin"]]
    """The view class for this aligned mapping."""

    _actual_class: ClassVar[Type["AlignedActualMixin"]]
    """The actual class (which has it’s own data) for this aligned mapping."""

    def __repr__(self):
        return f"{type(self).__name__} with keys: {', '.join([k + f'[{str(self[k].dtype)}]' for k in self.keys()])}"

    def _ipython_key_completions_(self) -> List[str]:
        return list(self.keys())

    def _validate_value(self, val: V, key: str) -> V:
        """Raises an error if value is invalid"""
        if isinstance(val, AwkArray):
            warnings.warn(
                "Support for Awkward Arrays is currently experimental. "
                "Behavior may change in the future. Please report any issues you may encounter!",
                ExperimentalFeatureWarning,
                # stacklevel=3,
            )
            # Prevent from showing up every time an awkward array is used
            # You'd think `once` works, but it doesn't at the repl and in notebooks
            warnings.filterwarnings(
                "ignore",
                category=ExperimentalFeatureWarning,
                message="Support for Awkward Arrays is currently experimental.*",
            )
        for i, axis in enumerate(self.axes):
            if self.parent.shape[axis] != dim_len(val, i):
                right_shape = tuple(self.parent.shape[a] for a in self.axes)
                actual_shape = tuple(dim_len(val, a) for a, _ in enumerate(self.axes))
                if actual_shape[i] is None and isinstance(val, AwkArray):
                    raise ValueError(
                        f"The AwkwardArray is of variable length in dimension {i}.",
                        f"Try ak.to_regular(array, {i}) before including the array in AnnData",
                    )
                else:
                    raise ValueError(
                        f"Value passed for key {key!r} is of incorrect shape. "
                        f"Values of {self.attrname} must match dimensions "
                        f"{self.axes} of parent. Value had shape {actual_shape} while "
                        f"it should have had {right_shape}."
                    )

        if not self._allow_df and isinstance(val, pd.DataFrame):
            name = self.attrname.title().rstrip("s")
            val = ensure_df_homogeneous(val, f"{name} {key!r}")
        return val

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
    def is_view(self) -> bool:
        pass

    @property
    def parent(self) -> Union["anndata.AnnData", "raw.Raw"]:
        return self._parent

    def copy(self):
        d = self._actual_class(self.parent, self._axis)
        for k, v in self.items():
            if isinstance(v, AwkArray):
                # Shallow copy since awkward array buffers are immutable
                d[k] = copy(v)
            else:
                d[k] = v.copy()
        return d

    def _view(self, parent: "anndata.AnnData", subset_idx: I):
        """Returns a subset copy-on-write view of the object."""
        if parent.is_view or subset_idx is not None:  # and or or?
            return self._view_class(self, parent, subset_idx)
        return self

    @deprecated("dict(obj)")
    def as_dict(self) -> dict:
        return dict(self)


class AlignedViewMixin:
    parent: "anndata.AnnData"
    """Reference to parent AnnData view"""

    attrname: str
    """What attribute in the parent is this?"""

    parent_mapping: Mapping[str, V]
    """The object this is a view of."""

    is_view = True

    def __getitem__(self, key: str) -> V:
        return as_view(
            _subset(self.parent_mapping[key], self.subset_idx),
            ElementRef(self.parent, self.attrname, (key,)),
        )

    def __setitem__(self, key: str, value: V):
        value = self._validate_value(value, key)  # Validate before mutating
        warnings.warn(
            f"Setting element `.{self.attrname}['{key}']` of view, "
            "initializing view as actual.",
            ImplicitModificationWarning,
            stacklevel=2,
        )
        with view_update(self.parent, self.attrname, ()) as new_mapping:
            new_mapping[key] = value

    def __delitem__(self, key: str):
        _ = key in self  # Make sure it exists before bothering with a copy
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


class AlignedActualMixin:
    _data: Dict[str, V]
    """Underlying mapping to the data"""

    is_view = False

    def __getitem__(self, key: str) -> V:
        return self._data[key]

    def __setitem__(self, key: str, value: V):
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


class AxisArraysBase(AlignedMapping):
    """\
    Mapping of key→array-like,
    where array-like is aligned to an axis of parent AnnData.
    """

    _allow_df = True
    _dimnames = ("obs", "var")

    @property
    def attrname(self) -> str:
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
                df[f"{key}{icolumn + 1}"] = column
        return df

    def _validate_value(self, val: V, key: str) -> V:
        if (
            hasattr(val, "index")
            and isinstance(val.index, cabc.Collection)
            and not val.index.equals(self.dim_names)
        ):
            # Could probably also re-order index if it’s contained
            try:
                pd.testing.assert_index_equal(val.index, self.dim_names)
            except AssertionError as e:
                msg = f"value.index does not match parent’s axis {self.axes[0]} names:\n{e}"
                raise ValueError(msg) from None
            else:
                msg = "Index.equals and pd.testing.assert_index_equal disagree"
                raise AssertionError(msg)
        return super()._validate_value(val, key)

    @property
    def dim_names(self) -> pd.Index:
        return (self.parent.obs_names, self.parent.var_names)[self._axis]


class AxisArrays(AlignedActualMixin, AxisArraysBase):
    def __init__(
        self,
        parent: Union["anndata.AnnData", "raw.Raw"],
        axis: int,
        vals: Union[Mapping, AxisArraysBase, None] = None,
    ):
        self._parent = parent
        if axis not in (0, 1):
            raise ValueError()
        self._axis = axis
        self._data = dict()
        if vals is not None:
            self.update(vals)


class AxisArraysView(AlignedViewMixin, AxisArraysBase):
    def __init__(
        self,
        parent_mapping: AxisArraysBase,
        parent_view: "anndata.AnnData",
        subset_idx: OneDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx
        self._axis = parent_mapping._axis


AxisArraysBase._view_class = AxisArraysView
AxisArraysBase._actual_class = AxisArrays


class LayersBase(AlignedMapping):
    """\
    Mapping of key: array-like, where array-like is aligned to both axes of the
    parent anndata.
    """

    _allow_df = False
    attrname = "layers"
    axes = (0, 1)

    # TODO: I thought I had a more elegant solution to overiding this...
    def copy(self) -> "Layers":
        d = self._actual_class(self.parent)
        for k, v in self.items():
            d[k] = v.copy()
        return d


class Layers(AlignedActualMixin, LayersBase):
    def __init__(self, parent: "anndata.AnnData", vals: Optional[Mapping] = None):
        self._parent = parent
        self._data = dict()
        if vals is not None:
            self.update(vals)


class LayersView(AlignedViewMixin, LayersBase):
    def __init__(
        self,
        parent_mapping: LayersBase,
        parent_view: "anndata.AnnData",
        subset_idx: TwoDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = subset_idx


LayersBase._view_class = LayersView
LayersBase._actual_class = Layers


class PairwiseArraysBase(AlignedMapping):
    """\
    Mapping of key: array-like, where both axes of array-like are aligned to
    one axis of the parent anndata.
    """

    _allow_df = False
    _dimnames = ("obs", "var")

    @property
    def attrname(self) -> str:
        return f"{self.dim}p"

    @property
    def axes(self) -> Tuple[int, int]:
        """Axes of the parent this is aligned to"""
        return self._axis, self._axis

    @property
    def dim(self) -> str:
        """Name of the dimension this aligned to."""
        return self._dimnames[self._axis]


class PairwiseArrays(AlignedActualMixin, PairwiseArraysBase):
    def __init__(
        self,
        parent: "anndata.AnnData",
        axis: int,
        vals: Optional[Mapping] = None,
    ):
        self._parent = parent
        if axis not in (0, 1):
            raise ValueError()
        self._axis = axis
        self._data = dict()
        if vals is not None:
            self.update(vals)


class PairwiseArraysView(AlignedViewMixin, PairwiseArraysBase):
    def __init__(
        self,
        parent_mapping: PairwiseArraysBase,
        parent_view: "anndata.AnnData",
        subset_idx: OneDIdx,
    ):
        self.parent_mapping = parent_mapping
        self._parent = parent_view
        self.subset_idx = (subset_idx, subset_idx)
        self._axis = parent_mapping._axis


PairwiseArraysBase._view_class = PairwiseArraysView
PairwiseArraysBase._actual_class = PairwiseArrays
