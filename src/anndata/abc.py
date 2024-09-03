from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    import numpy as np


class _AbstractCSDataset(ABC):
    """Base for the public API for CSRDataset/CSCDataset."""

    format: ClassVar[Literal["csr", "csc"]]
    """The format of the sparse matrix."""

    shape: tuple[int, int]
    """Shape of the matrix."""

    dtype: np.dtype
    """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""

    # TODO: __getitem__ and __setitem__?


_sparse_dataset_doc = """\
On disk {format} sparse matrix.

Analogous to :class:`h5py.Dataset` or :class:`zarr.core.Array`, but for sparse matrices.
"""


class CSRDataset(_AbstractCSDataset, ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSR")
    format = "csr"


class CSCDataset(_AbstractCSDataset, ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSC")
    format = "csc"


__all__ = [
    "CSRDataset",
    "CSCDataset",
]
