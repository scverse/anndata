from __future__ import annotations

from abc import ABC

_sparse_dataset_doc = """\
    On disk {format} sparse matrix.

    Analogous to :class:`h5py.Dataset` or :class:`zarr.core.Array`, but for sparse matrices.

    Parameters
    ----------
    group
        The backing group store.
"""


class CSRDataset(ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSR")
    format = "csr"


class CSCDataset(ABC):
    __doc__ = _sparse_dataset_doc.format(format="CSC")
    format = "csc"


__all__ = [
    "CSRDataset",
    "CSCDataset",
]
