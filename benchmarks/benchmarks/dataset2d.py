from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr

import anndata as ad

if TYPE_CHECKING:
    from collections.abc import Callable


def make_alternating_mask(n):
    mask_alternating = np.ones(10_000, dtype=bool)
    for i in range(0, 10_000, n):
        mask_alternating[i] = False
    return mask_alternating


class Dataset2D:
    param_names = ("gen_store",)
    params = (
        (
            lambda: h5py.File(Path(__file__).parent / "data/df.h5", mode="w"),
            lambda: zarr.open(Path(__file__).parent / "data/df.h5", mode="w"),
        ),
    )

    def setup(self, gen_store: Callable[[], zarr.Group | h5py.File]):
        self.n_obs = 500
        df = pd.DataFrame(
            {"a": np.array(["a"] * self.n_obs)},
            index=[f"cell{i}" for i in range(self.n_obs)],
        )
        store = gen_store()
        ad.io.write_elem(store, "obs", df)
        self.ds = ad.experimental.read_elem_lazy(store["obs"])

    def time_getitem_slice(self, *_):
        self.ds[0 : (self.n_obs // 2)]

    def peakmem_getitem_slivce(self, *_):
        self.ds[0 : (self.n_obs // 2)]

    def time_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)]

    def peakmem_getitem_bool_mask(self, *_):
        self.ds.iloc[np.random.randint(0, self.n_obs, self.n_obs // 2)]
