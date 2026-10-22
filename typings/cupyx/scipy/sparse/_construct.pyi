# SPDX-License-Identifier: MPL-2.0
from collections.abc import Sequence
from typing import Literal

from numpy.typing import DTypeLike

from ._base import spmatrix

def hstack(
    blocks: Sequence[spmatrix],
    format: Literal["csr", "csc", "coo"] | None = None,
    dtype: DTypeLike | None = None,
) -> spmatrix: ...
def vstack(
    blocks: Sequence[spmatrix],
    format: Literal["csr", "csc", "coo"] | None = None,
    dtype: DTypeLike | None = None,
) -> spmatrix: ...
