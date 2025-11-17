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
    param_names = ("store_type", "chunks", "array_type")
    params = (
        ("zarr", "h5ad"),
        ((-1,), None),
        ("cat", "numeric", "string-array", "nullable-string-array"),
    )

    def setup_cache(self):
        n_obs = 1_000_000
        array_types = {
            "numeric": np.arange(n_obs),
            "string-array": np.array(["a"] * n_obs),
            "nullable-string-array": pd.array(
                ["a", pd.NA] * (n_obs // 2), dtype="string"
            ),
            "cat": pd.Categorical(np.array(["a"] * n_obs)),
        }
        for k, v in array_types.items():
            for store in [
                h5py.File(f"data_{k}.h5ad", mode="w"),
                zarr.open(f"data_{k}.zarr", mode="w", zarr_version=2),
            ]:
                df = pd.DataFrame({"a": v}, index=[f"cell{i}" for i in range(n_obs)])
                if writing_string_array_on_disk := (
                    isinstance(v, np.ndarray) and df["a"].dtype == "string"
                ):
                    df["a"] = df["a"].to_numpy()
                with ad.settings.override(allow_write_nullable_strings=True):
                    ad.io.write_elem(store, "df", df)
                if writing_string_array_on_disk:
                    assert store["df"]["a"].attrs["encoding-type"] == "string-array"

    def setup(
        self,
        store_type: Literal["zarr", "h5ad"],
        chunks: None | tuple[int],
        array_type: Literal["cat", "numeric", "string-array", "nullable-string-array"],
    ):
        self.store = (
            h5py.File(f"data_{array_type}.h5ad", mode="r")
            if store_type == "h5ad"
            else zarr.open(f"data_{array_type}.zarr")
        )
        self.ds = ad.experimental.read_elem_lazy(self.store["df"], chunks=chunks)
        self.n_obs = self.ds.shape[0]

    def time_read_lazy_default(self, *_):
        ad.experimental.read_elem_lazy(self.store["df"])

    def peakmem_read_lazy_default(self, *_):
        ad.experimental.read_elem_lazy(self.store["df"])

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

    def time_concat(self, *_):
        adatas = [ad.AnnData(obs=self.ds)] * 50
        ad.concat(adatas, join="outer")
