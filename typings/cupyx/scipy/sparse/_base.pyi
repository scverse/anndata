# SPDX-License-Identifier: MPL-2.0
from typing import Any, Literal, Self

import cupy.cuda
import numpy as np
import scipy.sparse as sps
from numpy.typing import NDArray

class spmatrix:
    dtype: np.dtype[Any]
    shape: tuple[int, int]
    ndim: int
    def toarray(
        self, order: Literal["C", "F"] | None = None, out: None = None
    ) -> cupy.ndarray: ...
    def copy(self) -> spmatrix: ...
    def __power__(self, other: int) -> Self: ...
    def __array__(self) -> NDArray[Any]: ...
    def get(self, stream: cupy.cuda.Stream | None = None) -> sps.spmatrix: ...
