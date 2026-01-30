from __future__ import annotations

from collections.abc import Sequence
from types import EllipsisType

import numpy as np
import pandas as pd
from numpy import ma
from numpy.typing import NDArray

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
    XDataArray,
    ZappyArray,
    ZarrArray,
)

__all__ = ["AxisStorable", "Index", "Index1D", "RWAble"]


_Index1DNorm = slice | NDArray[np.bool_] | NDArray[np.integer]
# TODO: pd.Index[???]
type Index1D = (
    # 0D index
    int
    | str
    | np.int64
    # normalized 1D idex
    | _Index1DNorm
    # different containers for mask, obs/varnames, or numerical index
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

type Index = (
    Index1D
    | EllipsisType
    | tuple[Index1D | EllipsisType, Index1D | EllipsisType]
    | tuple[Index1D, Index1D, EllipsisType]
    | tuple[EllipsisType, Index1D, Index1D]
    | tuple[Index1D, EllipsisType, Index1D]
    | CSMatrix
    | CSArray
)
"""Index an :class:`~anndata.AnnData` object can be sliced with."""

_XDataType = (
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
_ArrayDataStructureTypes = _XDataType | AwkArray | XDataArray
_InMemoryArrayOrScalarType = pd.DataFrame | np.number | str | _ArrayDataStructureTypes
type AxisStorable = (
    _InMemoryArrayOrScalarType | dict[str, "AxisStorable"] | list["AxisStorable"]
)
"""A serializable object, excluding :class:`anndata.AnnData` objects i.e., something that can be stored in `uns` or `obsm`."""

type RWAble = AxisStorable | AnnData | pd.Categorical | pd.api.extensions.ExtensionArray
"""A superset of :type:`anndata.typing.AxisStorable` (i.e., including :class:`anndata.AnnData`) which is everything can be read/written by :func:`anndata.io.read_elem` and :func:`anndata.io.write_elem`."""
