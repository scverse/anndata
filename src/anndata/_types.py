"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse

from anndata._core.sparse_dataset import BaseCompressedSparseDataset
from anndata.compat import (
    AwkArray,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    H5Group,
    SpArray,
    ZappyArray,
    ZarrArray,
    ZarrGroup,
)

if TYPE_CHECKING:
    from typing import TypeAlias

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

InMemoryArrayOrScalarType: TypeAlias = Union[
    NDArray,
    np.ma.MaskedArray,
    sparse.spmatrix,
    SpArray,
    H5Array,
    ZarrArray,
    ZappyArray,
    BaseCompressedSparseDataset,
    DaskArray,
    CupyArray,
    CupySparseMatrix,
    AwkArray,
    pd.DataFrame,
    np.number,
    str,
]

ArrayStorageType = Union[ZarrArray, H5Array]
GroupStorageType = Union[ZarrGroup, H5Group]
StorageType = Union[ArrayStorageType, GroupStorageType]
