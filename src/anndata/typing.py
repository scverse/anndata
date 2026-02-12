from __future__ import annotations

import sys
from collections.abc import Sequence
from types import EllipsisType
from typing import TYPE_CHECKING, Never, TypeVar

import numpy as np
import pandas as pd
from numpy import ma
from numpy.typing import NDArray

from anndata.types import SupportsArrayApi

from . import abc
from ._core.anndata import AnnData
from .compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    IndexManager,
    XDataArray,
    ZappyArray,
    ZarrArray,
)

if TYPE_CHECKING:
    from typing import TypeAlias


__all__ = ["AxisStorable", "Index", "Index1D", "RWAble"]

if TYPE_CHECKING or sys.version_info >= (3, 13):
    _M = TypeVar("_M", IndexManager, Never, default=Never)
else:
    _M = TypeVar("_M")

_Index1DNorm: TypeAlias = (  # noqa: UP040
    slice | NDArray[np.bool_] | NDArray[np.integer] | SupportsArrayApi | _M
)
# TODO: pd.Index[???]
type Index1D[_M: IndexManager] = (
    # 0D index
    int
    | str
    | np.int64
    # normalized 1D idex
    | _Index1DNorm[_M]
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
    | _M
)
"""Index each :class:`~anndata.AnnData` object’s axis can be sliced with."""

type Index[_M: IndexManager] = (
    Index1D[_M]
    | EllipsisType
    | tuple[Index1D[_M] | EllipsisType, Index1D[_M] | EllipsisType]
    | tuple[Index1D, Index1D[_M], EllipsisType]
    | tuple[EllipsisType, Index1D[_M], Index1D[_M]]
    | tuple[Index1D[_M], EllipsisType, Index1D[_M]]
    | CSMatrix
    | CSArray
    | _M
)
"""Index an :class:`~anndata.AnnData` object can be sliced with."""

_XDataType: TypeAlias = (  # noqa: UP040
    np.ndarray
    | ma.MaskedArray
    | CSMatrix
    | CSArray
    | H5Array
    | ZarrArray
    | ZappyArray
    | abc.CSRDataset
    | abc.CSCDataset
    | DaskArray
    | CupyArray
    | CupySparseMatrix
)
_ArrayDataStructureTypes: TypeAlias = _XDataType | AwkArray | XDataArray  # noqa: UP040
_InMemoryArrayOrScalarType: TypeAlias = (  # noqa: UP040
    pd.DataFrame | np.number | str | _ArrayDataStructureTypes
)
type AxisStorable = (
    _InMemoryArrayOrScalarType | dict[str, "AxisStorable"] | list["AxisStorable"]
)
"""A serializable object, excluding :class:`anndata.AnnData` objects i.e., something that can be stored in `uns` or `obsm`."""

type RWAble = AxisStorable | AnnData | pd.Categorical | pd.api.extensions.ExtensionArray
"""A superset of :type:`anndata.typing.AxisStorable` (i.e., including :class:`anndata.AnnData`) which is everything can be read/written by :func:`anndata.io.read_elem` and :func:`anndata.io.write_elem`."""
