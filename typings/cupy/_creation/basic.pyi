# SPDX-License-Identifier: MPL-2.0
from typing import Literal

from numpy.typing import ArrayLike, DTypeLike

from .._core import ndarray

def full(
    shape: int | tuple[int, ...],
    fill_value: ArrayLike | ndarray,
    dtype: DTypeLike | None = None,
    order: Literal["C", "F"] = "C",
) -> ndarray: ...
