# SPDX-License-Identifier: MPL-2.0
from collections.abc import Sequence
from typing import Literal

from numpy.typing import ArrayLike, DTypeLike

from .._core import ndarray

def full(
    shape: int | Sequence[int],
    fill_value: ArrayLike | ndarray,
    dtype: DTypeLike | None = None,
    order: Literal["C", "F"] = "C",
) -> ndarray: ...
