from __future__ import annotations

from types import MappingProxyType

import numpy as np
import zarr
from scipy import sparse

from anndata import AnnData
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import write_elem


def make_alternating_mask(n):
    mask_alternating = np.ones(10_000, dtype=bool)
    for i in range(0, 10_000, n):
        mask_alternating[i] = False
    return mask_alternating


class SparseCSRContiguousSlice:
    _slices = MappingProxyType(
        {
            "0:1000": slice(0, 1000),
            "0:9000": slice(0, 9000),
            ":9000:-1": slice(None, 9000, -1),
            "::-2": slice(None, None, 2),
            "array": np.array([0, 5000, 9999]),
            "arange": np.arange(0, 1000),
            "first": 0,
            "alternating": make_alternating_mask(10),
        }
    )
    params = (
        [
            (10_000, 10_000),
            # (10_000, 500)
        ],
        _slices.keys(),
    )
    param_names = ("shape", "slice")

    def setup(self, shape: tuple[int, int], slice: str):
        X = sparse.random(
            *shape, density=0.01, format="csr", random_state=np.random.default_rng(42)
        )
        self.slice = self._slices[slice]
        g = zarr.group()
        write_elem(g, "X", X)
        self.x = sparse_dataset(g["X"])
        self.adata = AnnData(self.x)

    def time_getitem(self, *_):
        self.x[self.slice]

    def peakmem_getitem(self, *_):
        self.x[self.slice]

    def time_getitem_adata(self, *_):
        self.adata[self.slice]

    def peakmem_getitem_adata(self, *_):
        self.adata[self.slice]
