# SPDX-License-Identifier: MPL-2.0
from typing import Any, Literal

import h5py
import zarr
from numpy.typing import ArrayLike, DTypeLike

from .._core import ndarray

def asarray(
    a: ArrayLike | h5py.Dataset | zarr.Array[Any],
    dtype: DTypeLike | None = None,
    order: Literal["C", "F", "A", "K"] | None = None,
    *,
    blocking: bool = False,
) -> ndarray: ...
def array(
    obj: ArrayLike | ndarray,
    dtype: DTypeLike | None = None,
    copy: bool = True,
    order: Literal["C", "F", "A", "K"] = "K",
    subok: bool = False,
    ndmin: int = 0,
    *,
    blocking: bool = False,
) -> ndarray: ...
