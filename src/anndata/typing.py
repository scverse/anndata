from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np
import pandas as pd
from numpy import ma
from scipy import sparse

from . import abc
from ._core.anndata import AnnData
from .compat import (
    AwkArray,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    SpArray,
    ZappyArray,
    ZarrArray,
)

if TYPE_CHECKING:
    from typing import TypeAlias


__all__ = ["RWAble", "AxisStorable"]


ArrayDataStructureType: TypeAlias = Union[
    np.ndarray,
    ma.MaskedArray,
    sparse.csr_matrix,
    sparse.csc_matrix,
    SpArray,
    AwkArray,
    H5Array,
    ZarrArray,
    ZappyArray,
    abc.CSRDataset,
    abc.CSCDataset,
    DaskArray,
    CupyArray,
    CupySparseMatrix,
]


InMemoryArrayOrScalarType: TypeAlias = Union[
    pd.DataFrame, np.number, str, ArrayDataStructureType
]


AxisStorable: TypeAlias = Union[
    InMemoryArrayOrScalarType, dict[str, "AxisStorable"], list["AxisStorable"]
]
"""A serializable object, excluding :class:`anndata.AnnData` objects i.e., something that can be stored in `uns` or `obsm`."""

RWAble: TypeAlias = Union[
    AxisStorable,
    AnnData,
    pd.Categorical,
    pd.api.extensions.ExtensionArray,
]
"""A superset of :type:`anndata.typing.AxisStorable` (i.e., including :class:`anndata.AnnData`) which is everything can be read/written by :func:`anndata.read_elem` and :func:`anndata.write_elem`."""
