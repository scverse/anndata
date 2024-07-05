"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd
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

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

DictElemType = (
    np.ndarray
    | np.ma.MaskedArray
    | sparse.spmatrix
    | SpArray
    | H5Array
    | ZarrArray
    | ZappyArray
    | BaseCompressedSparseDataset
    | DaskArray
    | CupyArray
    | CupySparseMatrix
    | AwkArray
    | pd.DataFrame
    | np.number
    | str
)

ArrayStorageType = Union[ZarrArray, H5Array]
GroupStorageType = Union[ZarrGroup, H5Group]
StorageType = Union[ArrayStorageType, GroupStorageType]
