from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    import numpy as np
    from scipy.sparse import csc_matrix, csr_matrix

    from .compat import Index, SpArray


__all__ = ["CSRDataset", "CSCDataset"]


class _AbstractCSDataset(ABC):
    """Base for the public API for CSRDataset/CSCDataset."""

    format: ClassVar[Literal["csr", "csc"]]
    """The format of the sparse matrix."""

    shape: tuple[int, int]
    """Shape of the matrix."""

    dtype: np.dtype
    """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""

    backend: Literal["zarr", "hdf5"]
    """Which file type is used on-disk."""

    @abstractmethod
    def __getitem__(self, index: Index) -> float | csr_matrix | csc_matrix | SpArray:
        """Load a slice or an element from the sparse dataset into memory.

        Parameters
        ----------
        index
            Index to load.

        Returns
        -------
        The desired data read off disk.
        """

    @abstractmethod
    def to_memory(self) -> csr_matrix | csc_matrix | SpArray:
        """Load the sparse dataset into memory.

        Returns
        -------
        The in-memory representation of the sparse dataset.
        """


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
