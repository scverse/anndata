from __future__ import annotations

import tempfile
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr

import anndata as ad

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal


class Dataset2D:
    param_names = ("gen_store", "chunks", "gen_array_key")
    arr_generator = MappingProxyType({
        "cat": lambda n_obs: pd.Categorical(np.array(["a"] * n_obs)),
        "numeric": lambda n_obs: np.arange(n_obs),
        "string-array": lambda n_obs: np.array(["a"] * n_obs),
        "nullable-string-array": lambda n_obs: pd.array(
            ["a", pd.NA] * (n_obs // 2), dtype="string"
        ),
    })
    params = (
        (
            lambda: h5py.File(Path(tempfile.mkdtemp()) / "data.h5ad", mode="w"),
            lambda: zarr.open(
                Path(tempfile.mkdtemp()) / "data.zarr", mode="w", zarr_version=2
            ),
        ),
        ((-1,), None),
        ("cat", "numeric", "string-array", "nullable-string-array"),
    )

    def setup(
        self,
        gen_store: Callable[[], zarr.Group | h5py.File],
        chunks: None | tuple[int],
        gen_array_key: Literal[
            "cat", "numeric", "string-array", "nullable-string-array"
        ],
    ):
        self.n_obs = 100000
        df = pd.DataFrame(
            {"a": self.arr_generator[gen_array_key](self.n_obs)},
            index=[f"cell{i}" for i in range(self.n_obs)],
        )
        # Force writing string-array format to check read performance
        if writing_string_array_on_disk := (
            isinstance(self.arr_generator[gen_array_key](1), np.ndarray)
            and df["a"].dtype == "string"
        ):
            df["a"] = df["a"].to_numpy()
        self.store = gen_store()
        with ad.settings.override(allow_write_nullable_strings=True):
            ad.io.write_elem(self.store, "obs", df)
            if writing_string_array_on_disk:
                assert self.store["obs"]["a"].attrs["encoding-type"] == "string-array"
        self.ds = ad.experimental.read_elem_lazy(self.store["obs"], chunks=chunks)

    def time_read_lazy_default(self, *_):
        ad.experimental.read_elem_lazy(self.store["obs"])

    def peakmem_read_lazy_default(self, *_):
        ad.experimental.read_elem_lazy(self.store["obs"])

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
