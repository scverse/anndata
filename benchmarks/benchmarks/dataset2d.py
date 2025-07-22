from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr

import anndata as ad

if TYPE_CHECKING:
    from collections.abc import Callable


class Dataset2D:
    param_names = ("gen_store",)
    params = (
        (
            lambda: h5py.File(Path(tempfile.mkdtemp()) / "data.h5ad", mode="w"),
            lambda: zarr.open(
                Path(tempfile.mkdtemp()) / "data.zarr", mode="w", zarr_version=2
            ),
        ),
    )

    def setup(self, gen_store: Callable[[], zarr.Group | h5py.File]):
        self.n_obs = 10000
        df = pd.DataFrame(
            {"a": pd.Categorical(np.array(["a"] * self.n_obs))},
            index=[f"cell{i}" for i in range(self.n_obs)],
        )
        store = gen_store()
        ad.io.write_elem(store, "obs", df)
        self.ds = ad.experimental.read_elem_lazy(store["obs"])

    def time_getitem_slice(self, *_):
        self.ds.iloc[0 : (self.n_obs // 2)].to_memory()

    def peakmem_getitem_slivce(self, *_):
        self.ds.iloc[0 : (self.n_obs // 2)].to_memory()

    def time_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)].to_memory()

    def peakmem_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)].to_memory()
