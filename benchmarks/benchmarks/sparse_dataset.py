from __future__ import annotations

import numpy as np
import zarr
from scipy import sparse

from anndata import AnnData
from anndata.experimental import sparse_dataset, write_elem


def make_alternating_mask(n):
    mask_alternating = np.ones(10_000, dtype=bool)
    for i in range(0, 10_000, n):
        mask_alternating[i] = False
    return mask_alternating


class SparseCSRContiguousSlice:
    params = (
        [
            (10_000, 10_000),
            # (10_000, 500)
        ],
        [
            slice(0, 1000),
            slice(0, 9000),
            slice(None, 9000, -1),
            slice(None, None, 2),
            make_alternating_mask(10),
        ],
    )
    param_names = ["shape", "slice"]

    def setup(self, shape, slice):
        X = sparse.random(
            *shape, density=0.01, format="csr", random_state=np.random.default_rng(42)
        )
        self.slice = slice
        g = zarr.group()
        write_elem(g, "X", X)
        self.x = sparse_dataset(g["X"])
        self.adata = AnnData(self.x)

    def time_getitem(self, shape, slice):
        self.x[self.slice]

    def peakmem_getitem(self, shape, slice):
        self.x[self.slice]

    def time_getitem_adata(self, shape, slice):
        self.adata[self.slice]

    def peakmem_getitem_adata(self, shape, slice):
        self.adata[self.slice]
