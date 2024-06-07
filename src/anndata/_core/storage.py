from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from numpy import ma
from scipy import sparse

from ..compat import (
    AwkArray,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    SpArray,
    ZappyArray,
    ZarrArray,
)
from ..utils import ensure_df_homogeneous, join_english
from .sparse_dataset import BaseCompressedSparseDataset

if TYPE_CHECKING:
    from collections.abc import Generator


class StorageType(Enum):
    Array = (np.ndarray, "np.ndarray")
    Masked = (ma.MaskedArray, "numpy.ma.core.MaskedArray")
    Sparse = (sparse.spmatrix, "scipy.sparse.spmatrix")
    ZarrArray = (ZarrArray, "zarr.Array")
    ZappyArray = (ZappyArray, "zappy.base.ZappyArray")
    DaskArray = (DaskArray, "dask.array.Array")
    CupyArray = (CupyArray, "cupy.ndarray")
    CupySparseMatrix = (CupySparseMatrix, "cupyx.scipy.sparse.spmatrix")
    BackedSparseMatrix = (
        BaseCompressedSparseDataset,
        "anndata.experimental.[CSC,CSR]Dataset",
    )
    SparseArray = (SpArray, "scipy.sparse.sparray")
    AwkArray = (AwkArray, "awkward.Array")

    @property
    def cls(self):
        return self.value[0]

    @property
    def qualname(self):
        return self.value[1]

    @classmethod
    def classes(cls) -> tuple[type, ...]:
        return tuple(v.cls for v in cls)

    @classmethod
    def qualnames(cls) -> Generator[str, None, None]:
        yield from (v.qualname for v in cls)


def coerce_array(
    value: Any,
    *,
    name: str,
    allow_df: bool = False,
    allow_array_like: bool = False,
):
    """Coerce arrays stored in layers/X, and aligned arrays ({obs,var}{m,p})."""
    # If value is a scalar and we allow that, return it
    if allow_array_like and np.isscalar(value):
        return value
    # If value is one of the allowed types, return it
    if isinstance(value, StorageType.classes()):
        return value.A if isinstance(value, np.matrix) else value
    if isinstance(value, pd.DataFrame):
        return value if allow_df else ensure_df_homogeneous(value, name)
    # if value is an array-like object, try to convert it
    if allow_array_like:
        try:
            # TODO: asarray? asanyarray?
            return np.array(value)
        except (ValueError, TypeError):
            pass
    # if value isn’t the right type or convertible, raise an error
    msg = f"X needs to be of one of {join_english(StorageType.qualnames())}, not {type(value)}."
    raise ValueError(msg)
