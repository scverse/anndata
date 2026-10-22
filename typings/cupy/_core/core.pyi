# SPDX-License-Identifier: MPL-2.0
from collections.abc import Iterator
from types import EllipsisType
from typing import Any, Literal, Self, overload

import numpy as np
from cupy.cuda import MemoryPointer, Stream
from numpy._core.multiarray import flagsobj
from numpy.typing import DTypeLike, NDArray

class ndarray:
    dtype: np.dtype[Any]
    shape: tuple[int, ...]
    size: int
    ndim: int
    flags: flagsobj
    nbytes: int

    def __init__(
        self,
        shape: tuple[int, ...],
        dtype: DTypeLike | None = ...,
        memptr: MemoryPointer | None = None,
        strides: tuple[int, ...] | None = None,
        order: Literal["C", "F"] = "C",
    ) -> None: ...

    # cupy-specific
    def get(
        self,
        stream: Stream | None = None,
        order: Literal["C", "F", "A"] = "C",
        out: NDArray[Any] | None = None,
        blocking: bool = True,
    ) -> NDArray[Any]: ...

    # operators
    def __array__(self) -> NDArray[Any]: ...
    def __len__(self) -> int: ...
    def __getitem__(  # never returns scalars
        self,
        index: int
        | slice
        | EllipsisType
        | NDArray[np.integer | np.bool]
        | ndarray
        | tuple[int | slice | EllipsisType | None, ...],
    ) -> Self: ...
    def __iter__(self) -> Iterator[Self]: ...
    def __add__(self, other: object) -> Self: ...
    def __sub__(self, other: object) -> Self: ...
    def __eq__(self, value: object) -> ndarray: ...  # type: ignore[override]
    def __power__(self, other: int) -> Self: ...

    # methods
    def copy(self, order: Literal["C", "F", "A", "K"] = "C") -> ndarray: ...
    def item(self) -> bool | int | float | complex: ...
    def view(
        self, dtype: DTypeLike | None = None, type: type[ndarray] | None = None
    ) -> ndarray: ...
    def astype(
        self,
        dtype: DTypeLike | None,
        order: Literal["C", "F", "A", "K"] = "K",
        casting: None = None,
        subok: None = None,
        copy: bool = True,
    ) -> Self: ...
    @property
    def T(self) -> Self: ...
    @overload
    def all(self, axis: None = None) -> np.bool: ...
    @overload
    def all(self, axis: int) -> ndarray: ...
    def reshape(self, shape: tuple[int, ...] | int) -> ndarray: ...
    def squeeze(self, axis: int | None = None) -> Self: ...
    def ravel(self, order: Literal["C", "F", "A", "K"] = "C") -> Self: ...
    def flatten(self, order: Literal["C", "F", "A", "K"] = "C") -> Self: ...
    @property
    def flat(self) -> _FlatIter: ...
    def sum(
        self,
        axis: int | None = None,
        dtype: DTypeLike | None = None,
        out: ndarray | None = None,
        keepdims: bool = False,
    ) -> ndarray: ...

class _FlatIter:
    def __next__(self) -> np.float32 | np.float64: ...
    def __iter__(self) -> _FlatIter: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> np.float32 | np.float64: ...
