from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy import ma

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
    ZappyArray,
    ZarrArray,
)
from .compat import Index as _Index

if TYPE_CHECKING:
    from typing import TypeAlias


__all__ = ["Index", "RWAble", "AxisStorable"]


Index = _Index
"""1D or 2D index an :class:`~anndata.AnnData` object can be sliced with."""

ArrayDataStructureType: TypeAlias = (
    np.ndarray
    | ma.MaskedArray
    | CSMatrix
    | CSArray
    | AwkArray
    | H5Array
    | ZarrArray
    | ZappyArray
    | abc.CSRDataset
    | abc.CSCDataset
    | DaskArray
    | CupyArray
    | CupySparseMatrix
)


InMemoryArrayOrScalarType: TypeAlias = (
    pd.DataFrame | np.number | str | ArrayDataStructureType
)


AxisStorable: TypeAlias = (
    InMemoryArrayOrScalarType | dict[str, "AxisStorable"] | list["AxisStorable"]
)
"""A serializable object, excluding :class:`anndata.AnnData` objects i.e., something that can be stored in `uns` or `obsm`."""

RWAble: TypeAlias = (
    AxisStorable | AnnData | pd.Categorical | pd.api.extensions.ExtensionArray
)
"""A superset of :type:`anndata.typing.AxisStorable` (i.e., including :class:`anndata.AnnData`) which is everything can be read/written by :func:`anndata.io.read_elem` and :func:`anndata.io.write_elem`."""
