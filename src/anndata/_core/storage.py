from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Union, get_args

import numpy as np
import pandas as pd
from numpy import ma
from scipy import sparse

from .._warnings import ImplicitModificationWarning
from ..compat import (
    AwkArray,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    SpArray,
    ZappyArray,
    ZarrArray,
)
from ..utils import (
    ensure_df_homogeneous,
    join_english,
    raise_value_error_if_multiindex_columns,
)
from .sparse_dataset import CSCDataset, CSRDataset

if TYPE_CHECKING:
    from typing import Any, TypeAlias

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
    CSRDataset,
    CSCDataset,
    DaskArray,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    CupySparseMatrix,
]


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
    array_data_structure_types = get_args(ArrayDataStructureType)
    if isinstance(value, array_data_structure_types):
        if isinstance(value, np.matrix):
            msg = f"{name} should not be a np.matrix, use np.ndarray instead."
            warnings.warn(msg, ImplicitModificationWarning)
            value = value.A
        return value
    elif isinstance(value, sparse.spmatrix):
        msg = (
            f"AnnData previously had undefined behavior around matrices of type {type(value)}."
            "In 0.12, passing in this type will throw an error. Please convert to a supported type."
            "Continue using for this minor version at your own risk."
        )
        warnings.warn(msg, FutureWarning)
        return value
    if isinstance(value, pd.DataFrame):
        if allow_df:
            raise_value_error_if_multiindex_columns(value, name)
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
    msg = f"{name} needs to be of one of {join_english(map(str, array_data_structure_types))}, not {type(value)}."
    if e is not None:
        msg += " (Failed to convert it to an array, see above for details.)"
    raise ValueError(msg) from e
