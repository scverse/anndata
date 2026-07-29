from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    import numpy as np

    from ._types import _GroupStorageType
    from .compat import CSArray, CSMatrix
    from .typing import Index


__all__ = ["CSCDataset", "CSRDataset"]


class _AbstractCSDataset(ABC):
    """Base for the public API for CSRDataset/CSCDataset."""

    format: ClassVar[Literal["csr", "csc"]]
    """The format of the sparse matrix."""

    @property
    @abstractmethod
    def shape(self) -> tuple[int, int]:
        """Shape of the matrix."""

    @property
    @abstractmethod
    def dtype(self) -> np.dtype:
        """The :class:`numpy.dtype` of the `data` attribute of the sparse matrix."""

    @property
    @abstractmethod
    def backend(self) -> Literal["zarr", "hdf5"]:
        """Which file type is used on-disk."""

    @property
    @abstractmethod
    def group(self) -> _GroupStorageType:
        """The group underlying the backed matrix."""

    @abstractmethod
    def __getitem__(self, index: Index) -> float | CSMatrix | CSArray:
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
    def to_memory(self) -> CSMatrix | CSArray:
        """Load the sparse dataset into memory.

        Returns
        -------
        The in-memory representation of the sparse dataset.
        """

    @abstractmethod
    def append(self, sparse_matrix: CSMatrix | CSArray | _AbstractCSDataset) -> None:
        """Append an in-memory or on-disk sparse matrix to the current object's store.

        Parameters
        ----------
        sparse_matrix
            The matrix to append.

        Raises
        ------
        NotImplementedError
            If the matrix to append is not one of :class:`~scipy.sparse.csr_array`, :class:`~scipy.sparse.csc_array`, :class:`~scipy.sparse.csr_matrix`, or :class:`~scipy.sparse.csc_matrix`.
        ValueError
            If both the on-disk and to-append matrices are not of the same format i.e., `csr` or `csc`.
        OverflowError
            If the underlying data store has a 32 bit indptr, and the new matrix is too large to fit in it i.e., would cause a 64 bit `indptr` to be written.
        AssertionError
            If the on-disk data does not have `csc` or `csr` format.
        """


_sparse_dataset_doc = """\
On disk {format} sparse matrix.

Analogous to :class:`h5py.Dataset` or :class:`zarr.Array`, but for sparse matrices.
"""


def _redeclare_abstract_methods[T: type[_AbstractCSDataset]](cls: T) -> T:
    """Rebind the abstract methods so Sphinx doesn’t interpret them as inherited."""
    for name in ("__getitem__", "to_memory", "append"):
        setattr(cls, name, getattr(_AbstractCSDataset, name))
    return cls


@_redeclare_abstract_methods
class CSRDataset(_AbstractCSDataset, ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSR")
    format = "csr"


@_redeclare_abstract_methods
class CSCDataset(_AbstractCSDataset, ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSC")
    format = "csc"
