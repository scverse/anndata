from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

import anndata as ad

if TYPE_CHECKING:
    from numpy.typing import NDArray

    type Index = tuple[NDArray[np.bool_] | slice, ...]

N_VAR = 2_000


def make_index(name: str, n_obs: int, n_var: int) -> Index:
    rng = np.random.default_rng(0)
    obs_mask = rng.random(n_obs) < 0.5
    var_mask = rng.random(n_var) < 0.1
    indexers: dict[str, Index] = {
        "obs_mask": (obs_mask, slice(None)),
        "var_mask": (slice(None), var_mask),
    }
    return indexers[name]


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
