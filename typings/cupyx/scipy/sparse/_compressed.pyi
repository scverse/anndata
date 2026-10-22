# SPDX-License-Identifier: MPL-2.0
from typing import Literal, Self, overload

import numpy as np
import scipy.sparse as sps
from cupy import ndarray
from numpy.typing import DTypeLike, NDArray

from ._base import spmatrix

# `list`, not `Sequence`: a `tuple` means one index per axis
type _NonIntIndex = (
    slice | list[int] | list[bool] | ndarray | NDArray[np.integer | np.bool]
)
type _Index = int | _NonIntIndex

class _compressed_sparse_matrix(spmatrix):
    format: Literal["csr", "csc"]
    data: ndarray
    indices: ndarray
    indptr: ndarray

    @overload
    def __init__(
        self, arg1: ndarray | spmatrix | sps.spmatrix | sps.sparray
    ) -> None: ...
    @overload
    def __init__(
        self, arg1: tuple[int, int], *, dtype: DTypeLike | None = None
    ) -> None: ...
    @overload
    def __init__(self, arg1: tuple[ndarray, tuple[ndarray, ndarray]]) -> None: ...
    @overload
    def __init__(
        self,
        arg1: tuple[ndarray, ndarray, ndarray],
        shape: tuple[int, int] | None = None,
    ) -> None: ...
    @overload
    def __getitem__(self, key: tuple[int, int]) -> ndarray: ...  # 0-dimensional
    @overload
    def __getitem__(
        self, key: _Index | tuple[_NonIntIndex, _Index] | tuple[int, _NonIntIndex]
    ) -> Self: ...
    def __setitem__(
        self, key: _Index | tuple[_Index, _Index], value: object
    ) -> None: ...

    # methods
    def astype(self, dtype: DTypeLike | None) -> Self: ...
    def power(self, n: int, dtype: DTypeLike | None = None) -> Self: ...
    def sum(
        self,
        axis: Literal[0, 1] | None = None,
        dtype: DTypeLike | None = None,
        out: Self | None = None,
    ) -> ndarray: ...
