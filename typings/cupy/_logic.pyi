# SPDX-License-Identifier: MPL-2.0
from ._core import ndarray

def array_equal(a1: ndarray, a2: ndarray, equal_nan: bool = False) -> ndarray: ...
