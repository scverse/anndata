# SPDX-License-Identifier: MPL-2.0
from collections.abc import Sequence
from typing import Literal

from numpy.typing import DTypeLike

from ._core import ndarray

def broadcast_to(array: ndarray, shape: int | tuple[int, ...]) -> ndarray: ...
def concatenate(
    tup: Sequence[ndarray],
    axis: int | None = 0,
    out: ndarray | None = None,
    *,
    dtype: DTypeLike | None = None,
    casting: Literal["no", "equiv", "safe", "same_kind", "unsafe"] = "same_kind",
) -> ndarray: ...
