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
    CupySparseMatrix,
    DaskArray,
    IndexManager,
    XDataArray,
    ZappyArray,
)

__all__ = ["AlignedArray", "Index", "Index1D", "RWAble", "Storable"]


if TYPE_CHECKING:
    from typing import TypeAlias

# private aliases are no `type`s, so Sphinx expands them in the public ones
_M = TypeVar("_M", bound=IndexManager, default=Never)
_Index1DNorm: TypeAlias = (
    slice | NDArray[np.bool_] | NDArray[np.integer] | SupportsArrayApiBase | _M
)
# TODO: pd.Index[???]
type Index1D[M: IndexManager = Never] = (
    # 0D index
    int
    | str
    | np.int64
    # normalized 1D idex
    | _Index1DNorm[M]
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
    | M
)
"""Index each :class:`~anndata.AnnData` object’s axis can be sliced with."""

type Index[M: IndexManager = Never] = (
    Index1D[M]
    | EllipsisType
    | tuple[Index1D[M] | EllipsisType]
    | tuple[Index1D[M] | EllipsisType, Index1D[M] | EllipsisType]
    | tuple[Index1D[M], Index1D[M], EllipsisType]
    | tuple[EllipsisType, Index1D[M], Index1D[M]]
    | tuple[Index1D[M], EllipsisType, Index1D[M]]
    | CSMatrix
    | CSArray
    | M
)
"""Index an :class:`~anndata.AnnData` object can be sliced with."""

type InMemoryArray = (
    np.ndarray
    | ma.MaskedArray
    | CSMatrix
    | CSArray
    | DaskArray
    | CupyArray
    | CupySparseMatrix
    | SupportsArrayApiBase
)
"""An Array that is possibly stored in Memory (Dask Arrays are possibly stored on disk)."""

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
