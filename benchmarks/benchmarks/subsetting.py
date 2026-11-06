from __future__ import annotations

from typing import TYPE_CHECKING

import dask.array as da
import numpy as np
import pandas as pd
import zarr
from scipy import sparse

import anndata as ad
from anndata.experimental import read_elem_lazy

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray

    type Index = tuple[NDArray[np.bool_] | NDArray[np.intp] | slice, ...]

N_OBS, N_VAR = 20_000, 2_000
CHUNKS = {"rows": (2_500, -1), "cols": (-1, 250), "both": (2_500, 250)}


def make_index(name: str, n_obs: int, n_var: int) -> Index:
    rng = np.random.default_rng(0)
    obs_mask = rng.random(n_obs) < 0.5
    var_mask = rng.random(n_var) < 0.1
    indexers: dict[str, Index] = {
        "obs_mask": (obs_mask, slice(None)),
        "obs_sorted": (np.flatnonzero(rng.random(n_obs) < 0.1), slice(None)),
        "var_mask": (slice(None), var_mask),
        "both_masks": (obs_mask, var_mask),
    }
    return indexers[name]


class DaskSubset:
    """Subsetting dask arrays held in memory, chunked along either or both axes."""

    params = (
        ("dense", "csr", "csc"),
        tuple(CHUNKS),
        ("obs_mask", "obs_sorted", "var_mask", "both_masks"),
    )
    param_names = ("array_type", "chunks", "index")

    def setup(
        self,
        array_type: Literal["dense", "csr", "csc"],
        chunks: Literal["rows", "cols", "both"],
        index: str,
    ):
        # sparse matrices are chunked along their major axis
        if (array_type == "csr" and chunks != "rows") or (
            array_type == "csc" and chunks != "cols"
        ):
            raise NotImplementedError
        x = np.random.default_rng(0).random((N_OBS, N_VAR), dtype=np.float32)
        if array_type != "dense":
            x[x < 0.9] = 0
        X = da.from_array(x, chunks=CHUNKS[chunks])
        if array_type != "dense":
            cls = getattr(sparse, f"{array_type}_matrix")
            meta = cls((0, 0), dtype=X.dtype)
            X = X.map_blocks(cls, dtype=X.dtype, meta=meta).persist()
        self.adata = ad.AnnData(X)
        self.index = make_index(index, N_OBS, N_VAR)

    def time_subset(self, *_):
        self.adata[self.index].X.compute()

    def peakmem_subset(self, *_):
        self.adata[self.index].X.compute()


class DaskSubsetLazyDense:
    """Subsetting a dense matrix read lazily from zarr."""

    filepath = "dense.zarr"
    params = (("obs_mask", "var_mask", "both_masks"),)
    param_names = ("index",)

    def setup_cache(self):
        g = zarr.open_group(self.filepath, mode="w")
        x = np.random.default_rng(0).random((N_OBS, N_VAR), dtype=np.float32)
        ad.io.write_elem(g, "X", x)

    def setup(self, index: str):
        g = zarr.open_group(self.filepath, mode="r")
        self.adata = ad.AnnData(read_elem_lazy(g["X"], chunks=CHUNKS["rows"]))
        self.index = make_index(index, N_OBS, N_VAR)

    def time_subset(self, *_):
        self.adata[self.index].X.compute()

    def peakmem_subset(self, *_):
        self.adata[self.index].X.compute()


class SubsetCategorical:
    """Subsetting with categorical annotations, whose unused categories are removed."""

    params = (("obs_mask", "var_mask"),)
    param_names = ("index",)

    def setup(self, index: str):
        n_obs = 500_000
        rng = np.random.default_rng(0)
        obs = pd.DataFrame(
            {
                # every category is used, so subsetting only genes keeps all of them
                f"cat_{n}": pd.Categorical.from_codes(
                    rng.permutation(np.arange(n_obs) % n),
                    categories=pd.Index([f"c{i}" for i in range(n)]),
                )
                for n in (10, 1_000, 100_000)
            },
            index=pd.RangeIndex(n_obs).astype(str),
        )
        var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VAR)])
        self.adata = ad.AnnData(obs=obs, var=var)
        self.index = make_index(index, n_obs, N_VAR)

    def time_view(self, *_):
        self.adata[self.index]

    def time_copy(self, *_):
        self.adata[self.index].copy()

    def peakmem_copy(self, *_):
        self.adata[self.index].copy()
