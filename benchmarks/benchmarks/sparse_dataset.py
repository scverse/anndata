from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import zarr
from dask.array.core import Array as DaskArray
from scipy import sparse

from anndata import AnnData, concat
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import write_elem
from anndata.experimental import read_elem_lazy

if TYPE_CHECKING:
    from typing import Literal


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

    def setup_cache(self):
        X = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        g = zarr.group(self.filepath)
        write_elem(g, "X", X)

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


class SparseCSRDaskConcat:
    filepath = "data.zarr"

    params = (["inner", "outer"], [0, -1])
    param_names = ("join", "fill_value")

    def setup_cache(self):
        X = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        g = zarr.group(self.filepath)
        write_elem(g, "X", X)

    def setup(self, *_):
        self.group = zarr.group(self.filepath)
        self.adatas = [
            AnnData(
                var=pd.DataFrame(
                    index=[
                        f"gene_{j}{f'_{i}' if (j % 100 == 0) else ''}"
                        for j in range(10_000)
                    ]
                ),
                X=read_elem_lazy(self.group["X"]),
            )
            for i in range(10)
        ]

    def time_concat(self, join: Literal["inner", "outer"], fill_value: Literal[0, -1]):
        concat(self.adatas, join=join, fill_value=fill_value)

    def peakmem_concat(
        self, join: Literal["inner", "outer"], fill_value: Literal[0, -1]
    ):
        concat(self.adatas, join=join, fill_value=fill_value)

    def time_concat_with_mem(
        self, join: Literal["inner", "outer"], fill_value: Literal[0, -1]
    ):
        concat(self.adatas, join=join, fill_value=fill_value).to_memory()

    def peakmem_concat_with_mem(
        self, join: Literal["inner", "outer"], fill_value: Literal[0, -1]
    ):
        concat(self.adatas, join=join, fill_value=fill_value).to_memory()


class SparseCSRDask:
    filepath = "data.zarr"

    def setup_cache(self):
        X = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        g = zarr.group(self.filepath)
        write_elem(g, "X", X)

    def setup(self, *_):
        self.group = zarr.group(self.filepath)

    def time_read(self, *_):
        AnnData(X=read_elem_lazy(self.group["X"]))

    def peakmem_read(self, *_):
        AnnData(X=read_elem_lazy(self.group["X"]))
