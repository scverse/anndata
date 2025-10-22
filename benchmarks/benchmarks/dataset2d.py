from __future__ import annotations

from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr

import anndata as ad

if TYPE_CHECKING:
    from typing import Literal


class Dataset2D:
    param_names = ("store_type", "chunks")
    params = (
        ("zarr", "h5ad"),
        ((-1,), None),
    )

    def setup_cache(self):
        n_obs = 100000
        df = pd.DataFrame(
            {
                "a": pd.Categorical(np.array(["a"] * n_obs)),
                "b": np.arange(n_obs),
            },
            index=[f"cell{i}" for i in range(n_obs)],
        )
        for store in [
            h5py.File("data.h5ad", mode="w"),
            zarr.open("data.zarr", mode="w", zarr_version=2),
        ]:
            ad.io.write_elem(store, "obs", df)

    def setup(self, store_type: Literal["zarr", "h5ad"], chunks: None | tuple[int]):
        store = (
            h5py.File("data.h5ad", mode="r")
            if store_type == "h5ad"
            else zarr.open("data.zarr")
        )
        self.ds = ad.experimental.read_elem_lazy(store["obs"], chunks=chunks)
        self.n_obs = self.ds.shape[0]

    def time_getitem_slice(self, *_):
        self.ds.iloc[0 : (self.n_obs // 2)].to_memory()

    def peakmem_getitem_slice(self, *_):
        self.ds.iloc[0 : (self.n_obs // 2)].to_memory()

    def time_full_to_memory(self, *_):
        self.ds.to_memory()

    def peakmem_full_to_memory(self, *_):
        self.ds.to_memory()

    def time_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)].to_memory()

    def peakmem_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)].to_memory()
