from __future__ import annotations

from types import MappingProxyType

import numpy as np
import zarr
from dask.array.core import Array as DaskArray
from scipy import sparse

from anndata import AnnData, concat
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import write_elem
from anndata.experimental import read_elem_lazy


def make_alternating_mask(n):
    mask_alternating = np.ones(10_000, dtype=bool)
    for i in range(0, 10_000, n):
        mask_alternating[i] = False
    return mask_alternating


class SparseCSRContiguousSlice:
    _indexers = MappingProxyType({
        "0:1000": slice(0, 1000),
        "0:9000": slice(0, 9000),
        ":9000:-1": slice(None, 9000, -1),
        "::-2": slice(None, None, 2),
        "array": np.array([0, 5000, 9999]),
        "arange": np.arange(0, 1000),
        "first": 0,
        "alternating": make_alternating_mask(10),
    })
    filepath = "data.zarr"
    params = (
        list(_indexers.keys()),
        [True, False],
    )
    param_names = (
        "index",
        "use_dask",
    )

    def setup_cache(self) -> None:
        x = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        g = zarr.group(self.filepath)
        write_elem(g, "X", x)

    def setup(self, index: str, use_dask: bool):  # noqa: FBT001
        g = zarr.open(self.filepath)
        self.x = read_elem_lazy(g["X"]) if use_dask else sparse_dataset(g["X"])
        self.adata = AnnData(self.x)
        self.index = self._indexers[index]

    def time_getitem(self, *_):
        res = self.x[self.index]
        if isinstance(res, DaskArray):
            res.compute()

    def peakmem_getitem(self, *_):
        res = self.x[self.index]
        if isinstance(res, DaskArray):
            res.compute()

    def time_getitem_adata(self, *_):
        res = self.adata[self.index]
        if isinstance(res, DaskArray):
            res.compute()

    def peakmem_getitem_adata(self, *_):
        res = self.adata[self.index]
        if isinstance(res, DaskArray):
            res.compute()


class SparseCSRDask:
    filepath = "data.zarr"

    def setup_cache(self) -> None:
        x = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        g = zarr.group(self.filepath)
        write_elem(g, "X", x)

    def setup(self):
        self.group = zarr.group(self.filepath)
        self.adata = AnnData(X=read_elem_lazy(self.group["X"]))

    def time_concat(self):
        concat([self.adata for i in range(100)])

    def peakmem_concat(self):
        concat([self.adata for i in range(100)])

    def time_read(self):
        AnnData(X=read_elem_lazy(self.group["X"]))

    def peakmem_read(self):
        AnnData(X=read_elem_lazy(self.group["X"]))
