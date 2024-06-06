from __future__ import annotations

from enum import Enum
from typing import Any

import numpy as np
import pandas as pd
from numpy import ma
from scipy import sparse

from ..compat import (
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    SpArray,
    ZappyArray,
    ZarrArray,
)
from ..utils import ensure_df_homogeneous
from .sparse_dataset import BaseCompressedSparseDataset


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sparse.spmatrix
    ZarrArray = ZarrArray
    ZappyArray = ZappyArray
    DaskArray = DaskArray
    CupyArray = CupyArray
    CupySparseMatrix = CupySparseMatrix
    BackedSparseMatrix = BaseCompressedSparseDataset
    SparseArray = SpArray

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


def coerce_array(value: Any, *, name: str, allow_df: bool = True):
    """Coerce arrays stored in layers/X, and aligned arrays ({obs,var}{m,p})."""
    if (
        isinstance(value, StorageType.classes()) or np.isscalar(value)
    ) and not isinstance(value, np.matrix):
        return value
    if not allow_df and isinstance(value, pd.DataFrame):
        return ensure_df_homogeneous(value, name)
    # TODO: asarray? asanyarray?
    return np.array(value)
