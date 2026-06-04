from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import sparse

from anndata.compat import CSArray, CSMatrix

from .._warnings import ImplicitModificationWarning
from ..compat import XDataset, has_xp
from ..utils import (
    ensure_df_homogeneous,
    get_union_members,
    join_english,
    raise_value_error_if_multiindex_columns,
    warn,
)
from .xarray import Dataset2D

if TYPE_CHECKING:
    from typing import Any

    from .._core.anndata import AnnData


def _spec_violation_message(value: Any, *, name: str) -> str | None:
    """Return a spec-violation message for higher-than-2D ``X``/``layers`` values.

    AnnData's on-disk specification requires ``X`` and every entry of
    ``layers`` to be two-dimensional. Returns ``None`` for conforming values
    (no ``shape`` attribute, ``None``, or ``ndim <= 2``).
    """
    if value is None:
        return None
    shape = getattr(value, "shape", None)
    if shape is None:
        return None
    ndim = len(shape)
    if ndim <= 2:
        return None
    return (
        f"{name} must be 2-dimensional, but got an array with shape "
        f"{tuple(shape)} (ndim={ndim}). Storing higher-dimensional arrays "
        f"in `X` or `layers` violates the AnnData specification, and cannot be written to disk."
    )


def _check_x_and_layers_are_2d_on_write(adata: AnnData) -> None:
    """Reject writing AnnData objects whose ``X`` or ``layers`` are not 2D.

    AnnData's spec requires ``X`` and ``layers`` entries on disk to be
    2-dimensional. In-memory we do not block higher-dimensional values, but
    a 3D+ array would propagate the spec violation onto disk, so we hard-fail
    at the IO boundary.
    """
    if (X := adata.layers.get(None)) is not None:
        msg = _spec_violation_message(X, name="X")
        if msg is not None:
            raise ValueError(msg)
    for key, value in adata.layers.items():
        if key is None:
            continue
        msg = _spec_violation_message(value, name=f"Layer {key!r}")
        if msg is not None:
            raise ValueError(msg)


def coerce_array(
    value: Any,
    *,
    name: str,
    allow_df: bool = False,
    allow_array_like: bool = False,
):
    """Coerce arrays stored in layers/X, and aligned arrays ({obs,var}{m,p})."""
    from ..typing import _ArrayDataStructureTypes

    # If value is a scalar and we allow that, return it
    if allow_array_like and np.isscalar(value):
        return value
    # If value is one of the allowed types, return it
    array_data_structure_types = get_union_members(_ArrayDataStructureTypes)
    if isinstance(value, XDataset):
        value = Dataset2D(value)
    if isinstance(value, (*array_data_structure_types, Dataset2D)):
        if isinstance(value, np.matrix):
            msg = f"{name} should not be a np.matrix, use np.ndarray instead."
            warn(msg, ImplicitModificationWarning)
            value = value.A
        return value
    if has_xp(value):
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
