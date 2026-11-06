from __future__ import annotations

from collections.abc import Sequence
from types import EllipsisType
from typing import TYPE_CHECKING, Never, TypeVar

import h5py
import numpy as np
import pandas as pd
import zarr
from numpy import ma
from numpy.typing import NDArray

from anndata.types import SupportsArrayApiBase

from . import abc
from ._core.anndata import AnnData
from ._core.xarray import Dataset2D
from .compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    CupyArray,
    CupyCSMatrix,
    DaskArray,
    IndexManager,
    XDataArray,
    ZappyArray,
)

__all__ = ["AlignedArray", "Index", "Index1D", "RWAble", "Storable"]


if TYPE_CHECKING:
    from typing import TypeAlias, TypeIs

# private aliases are no `type`s, so subscripting substitutes into the union and Sphinx expands them in the public ones
_M = TypeVar("_M", bound=IndexManager, default=Never)
_I = TypeVar("_I")
_Index1DNormBase: TypeAlias = (  # noqa: UP040
    slice | NDArray[np.bool_] | NDArray[np.integer] | SupportsArrayApiBase
)
_Index1DNorm: TypeAlias = _Index1DNormBase | _M  # noqa: PYI047
# TODO: pd.Index[???]
type Index1D = (
    # 0D index
    int
    | str
    | np.int64
    # normalized 1D idex
    | _Index1DNormBase
    # different containers for mask, obs/varnames, or numerical index,
    | Sequence[int]
    | Sequence[str]
    | Sequence[bool]
    | pd.Series  # bool, int, str
    | pd.Index
    | pd.api.extensions.ExtensionArray  # bool | int | str
    | NDArray[np.str_]
    | np.matrix  # bool
    | CSMatrix  # bool
    | CSArray  # bool
)
"""Index each :class:`~anndata.AnnData` object’s axis can be sliced with."""

_IndexOf: TypeAlias = (  # noqa: UP040
    _I
    | EllipsisType
    | tuple[_I | EllipsisType]
    | tuple[_I | EllipsisType, _I | EllipsisType]
    | tuple[_I, _I, EllipsisType]
    | tuple[EllipsisType, _I, _I]
    | tuple[_I, EllipsisType, _I]
    | CSMatrix
    | CSArray
)
type Index = _IndexOf[Index1D]
"""Index an :class:`~anndata.AnnData` object can be sliced with."""

# internal variants that also accept the `IndexManager`s stored by views
_Index1D: TypeAlias = Index1D | _M
_Index: TypeAlias = _IndexOf[_Index1D[_M]]  # noqa: PYI047

type InMemoryArray = (
    np.ndarray
    | ma.MaskedArray
    | CSMatrix
    | CSArray
    | DaskArray
    | CupyArray
    | CupyCSMatrix
    | SupportsArrayApiBase
)
"""An Array that is possibly stored in Memory (Dask Arrays are possibly stored on disk)."""


def _is_in_memory(a: object, /) -> TypeIs[InMemoryArray]:
    return isinstance(a, InMemoryArray.__value__)


type AlignedArray = (
    InMemoryArray
    | h5py.Dataset
    | zarr.Array
    | ZappyArray
    | abc.CSRDataset
    | abc.CSCDataset
    | AwkArray
    | XDataArray
    | Dataset2D
)
"""Every array-like data structure anndata accepts, i.e. the full set of array types that can be stored in `X`, `layers`, and `{obs,var}{m,p}`."""

_StorableArrayOrScalar: TypeAlias = AlignedArray | pd.DataFrame | np.number | str  # noqa: UP040
type Storable = _StorableArrayOrScalar | dict[str, "Storable"] | list["Storable"]
"""A serializable object, excluding :class:`anndata.AnnData` objects i.e., something that can be stored in `uns` or `obsm`."""

type RWAble = Storable | AnnData | pd.api.extensions.ExtensionArray
"""A superset of :type:`anndata.typing.Storable` (i.e., including :class:`anndata.AnnData`) which is everything can be read/written by :func:`anndata.io.read_elem` and :func:`anndata.io.write_elem`."""


def __getattr__(name: str) -> object:
    from ._warnings import warn

    if name == "AxisStorable":
        msg = "AxisStorable is deprecated, use Storable instead"
        warn(msg, FutureWarning)

    raise AttributeError(name)
