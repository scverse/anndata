from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import h5py
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


def make_integer_indexers() -> MappingProxyType:
    rng = np.random.default_rng(42)
    fragmented_sample = rng.choice(10_000, size=2_048, replace=False)
    return MappingProxyType({
        "single_run": np.arange(2_048),
        "multiple_runs": np.concatenate([
            np.arange(0, 512),
            np.arange(1_500, 2_012),
            np.arange(3_000, 3_512),
            np.arange(7_000, 7_512),
        ]),
        "fragmented_sorted": np.sort(fragmented_sample),
        "fragmented_unsorted": rng.permutation(fragmented_sample),
        "clustered_shuffled": rng.permutation(
            np.concatenate([
                np.arange(100, 612),
                np.arange(1_500, 2_012),
                np.arange(3_000, 3_512),
                np.arange(7_000, 7_512),
            ])
        ),
        "clustered_duplicates": rng.permutation(
            np.concatenate([
                np.repeat(np.arange(100, 356), 2),
                np.repeat(np.arange(1_500, 1_756), 2),
                np.repeat(np.arange(3_000, 3_256), 2),
                np.repeat(np.arange(7_000, 7_256), 2),
            ])
        ),
    })


INT_INDEXERS = make_integer_indexers()


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
        g = zarr.open(self.filepath, mode="w")
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


class SparseBackedIntegerIndexing:
    filepath = "data"
    params = (
        ("h5ad", "zarr"),
        ("csr", "csc"),
        list(INT_INDEXERS.keys()),
    )
    param_names = ("store_type", "sparse_format", "index_case")

    def setup_cache(self):
        rng = np.random.default_rng(42)
        csr = sparse.random(
            10_000,
            10_000,
            density=0.01,
            format="csr",
            random_state=rng,
        )
        csc = csr.tocsc()
        for store_type in ["h5ad", "zarr"]:
            path = f"{self.filepath}.{store_type}"
            if store_type == "h5ad":
                with h5py.File(path, mode="w") as f:
                    write_elem(f, "csr", csr)
                    write_elem(f, "csc", csc)
            else:
                g = zarr.open(path, mode="w")
                write_elem(g, "csr", csr)
                write_elem(g, "csc", csc)

    def setup(self, store_type: str, sparse_format: str, index_case: str):
        self._h5_file = None
        if store_type == "h5ad":
            self._h5_file = h5py.File(f"{self.filepath}.h5ad", mode="r")
            self.group = self._h5_file
        else:
            self.group = zarr.open(f"{self.filepath}.zarr", mode="r")
        self.x = sparse_dataset(self.group[sparse_format])
        self.index = INT_INDEXERS[index_case]
        self.is_csr = sparse_format == "csr"

    def teardown(self, *_):
        if self._h5_file is not None:
            self._h5_file.close()

    def time_getitem(self, *_):
        if self.is_csr:
            self.x[self.index, :]
        else:
            self.x[:, self.index]

    def peakmem_getitem(self, *_):
        if self.is_csr:
            self.x[self.index, :]
        else:
            self.x[:, self.index]


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
        g = zarr.open(self.filepath, mode="w")
        write_elem(g, "X", X)

    def setup(self, *_):
        self.group = zarr.group(self.filepath)
        self.adatas = [
            AnnData(
                var=pd.DataFrame(
                    index=[
                        f"gene_{j}{f'_{i}' if (j % 500 == 0) else ''}"
                        for j in range(10_000)
                    ]
                ),
                X=read_elem_lazy(self.group["X"]),
            )
            for i in range(5)
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
        g = zarr.open(self.filepath, mode="w")
        write_elem(g, "X", X)

    def setup(self, *_):
        self.group = zarr.group(self.filepath)

    def time_read(self, *_):
        AnnData(X=read_elem_lazy(self.group["X"]))

    def peakmem_read(self, *_):
        AnnData(X=read_elem_lazy(self.group["X"]))
