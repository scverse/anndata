from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, get_args

import numpy as np
import pandas as pd
from scipy import sparse

from anndata.compat import CSArray, CSMatrix

from .._warnings import ImplicitModificationWarning
from ..compat import XDataset
from ..utils import (
    ensure_df_homogeneous,
    join_english,
    raise_value_error_if_multiindex_columns,
)
from .xarray import Dataset2D

if TYPE_CHECKING:
    from typing import Any


def coerce_array(
    value: Any,
    *,
    name: str,
    allow_df: bool = False,
    allow_array_like: bool = False,
):
    """Coerce arrays stored in layers/X, and aligned arrays ({obs,var}{m,p})."""
    from ..typing import ArrayDataStructureTypes

    # If value is a scalar and we allow that, return it
    if allow_array_like and np.isscalar(value):
        return value
    # If value is one of the allowed types, return it
    array_data_structure_types = get_args(ArrayDataStructureTypes)
    if isinstance(value, XDataset) and not isinstance(value, Dataset2D):
        value = Dataset2D(value.data_vars, value.coords, value.attrs)
    if isinstance(value, (*array_data_structure_types, Dataset2D)):
        if isinstance(value, np.matrix):
            msg = f"{name} should not be a np.matrix, use np.ndarray instead."
            warnings.warn(msg, ImplicitModificationWarning, stacklevel=3)
            value = value.A
        return value
    is_non_csc_r_array_or_matrix = (
        (isinstance(value, base) and not isinstance(value, csr_c_format))
        for base, csr_c_format in [
            (sparse.spmatrix, CSMatrix),
            (sparse.sparray, CSArray),
        ]
    )
    if any(is_non_csc_r_array_or_matrix):
        msg = f"Only CSR and CSC {'matrices' if isinstance(value, sparse.spmatrix) else 'arrays'} are supported."
        raise ValueError(msg)
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
