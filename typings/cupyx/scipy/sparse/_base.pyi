# SPDX-License-Identifier: MPL-2.0
from typing import Any, Literal, Self

import cupy.cuda
import numpy as np
import scipy.sparse as sps
from cupy import ndarray
from numpy.typing import NDArray

from ._csc import csc_matrix
from ._csr import csr_matrix

# constructor inputs all sparse formats accept
type _ToSparse = (  # noqa: PYI047
    ndarray
    | spmatrix
    | sps.spmatrix
    | sps.sparray
    | tuple[int, int]
    | tuple[ndarray, tuple[ndarray, ndarray]]
)

class spmatrix:
    dtype: np.dtype[Any]
    shape: tuple[int, int]
    ndim: int
    def toarray(
        self, order: Literal["C", "F"] | None = None, out: None = None
    ) -> cupy.ndarray: ...
    def copy(self) -> spmatrix: ...
    def tocsr(self, copy: bool = False) -> csr_matrix: ...
    def tocsc(self, copy: bool = False) -> csc_matrix: ...
    def __matmul__(
        self, other: spmatrix
    ) -> csr_matrix: ...  # sparse @ sparse is always CSR
    def __power__(self, other: int) -> Self: ...
    def __array__(self) -> NDArray[Any]: ...
    def get(self, stream: cupy.cuda.Stream | None = None) -> sps.spmatrix: ...
