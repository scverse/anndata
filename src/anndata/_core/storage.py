from __future__ import annotations

import warnings
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy import ma
from scipy import sparse

from .._warnings import ImplicitModificationWarning
from ..compat import (
    AwkArray,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    ZappyArray,
    ZarrArray,
)
from ..utils import ensure_df_homogeneous, join_english
from .sparse_dataset import BaseCompressedSparseDataset

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Any


class StorageType(Enum):
    # Memory
    Array = (np.ndarray, "np.ndarray")
    Masked = (ma.MaskedArray, "numpy.ma.core.MaskedArray")
    Sparse = (sparse.spmatrix, "scipy.sparse.spmatrix")
    AwkArray = (AwkArray, "awkward.Array")
    # Backed
    HDF5Dataset = (H5Array, "h5py.Dataset")
    ZarrArray = (ZarrArray, "zarr.Array")
    ZappyArray = (ZappyArray, "zappy.base.ZappyArray")
    BackedSparseMatrix = (
        BaseCompressedSparseDataset,
        "anndata.experimental.[CSC,CSR]Dataset",
    )
    # Distributed
    DaskArray = (DaskArray, "dask.array.Array")
    CupyArray = (CupyArray, "cupy.ndarray")
    CupySparseMatrix = (CupySparseMatrix, "cupyx.scipy.sparse.spmatrix")

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
        if isinstance(value, np.matrix):
            msg = f"{name} should not be a np.matrix, use np.ndarray instead."
            warnings.warn(msg, ImplicitModificationWarning)
            value = value.A
        return value
    if isinstance(value, pd.DataFrame):
        return value if allow_df else ensure_df_homogeneous(value, name)
    # if value is an array-like object, try to convert it
    e = None
    if allow_array_like:
        try:
            # TODO: asarray? asanyarray?
            return np.array(value)
        except (ValueError, TypeError) as _e:
            e = _e
    # if value isn’t the right type or convertible, raise an error
    msg = f"{name} needs to be of one of {join_english(StorageType.qualnames())}, not {type(value)}."
    if e is not None:
        msg += " (Failed to convert it to an array, see above for details.)"
    raise ValueError(msg) from e
