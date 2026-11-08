from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata.acc import A
from anndata.tests.helpers import DEFAULT_KEY_TYPES, AwkArray, assert_equal, gen_adata

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from zarr.storage import MemoryStore

    from anndata.acc import AdRef


@pytest.fixture(scope="module")
def adata() -> ad.AnnData:
    types = (*DEFAULT_KEY_TYPES, AwkArray)
    adata = gen_adata((20, 15), obsm_types=types, varm_types=types)
    rng = np.random.default_rng(0)
    for dim in ("obs", "var"):
        df = getattr(adata, dim)
        df["keep"] = rng.random(len(df)) > 0.3
        df["keep2"] = rng.random(len(df)) > 0.3
        df["none"] = False
        strs = [f"s{i}" if i % 3 else None for i in range(len(df))]
        df["nstr"] = pd.array(strs, dtype=pd.StringDtype())
        df["nstr_nan"] = pd.array(strs, dtype=pd.StringDtype(na_value=np.nan))
    adata.raw = adata
    return adata


@pytest.mark.parametrize(
    "filter",
    [
        pytest.param(A.obs["keep"], id="obs"),
        pytest.param(A.var["keep"], id="var"),
        pytest.param([A.obs["keep"], A.var["keep"]], id="both"),
        pytest.param([A.obs["keep"], A.obs["keep2"]], id="obs-and"),
        pytest.param(A.obs["none"], id="empty"),
    ],
)
def test_read_filter(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    filter: AdRef | list[AdRef],
) -> None:
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    refs = filter if isinstance(filter, list) else [filter]
    masks = {"obs": np.ones(adata.n_obs, bool), "var": np.ones(adata.n_vars, bool)}
    for ref in refs:
        (dim,) = ref.dims
        masks[dim] &= ref.acc.get(adata)[ref.idx].to_numpy()

    read = getattr(ad.io, f"read_{diskfmt}")(diskfmt_store, filter=filter)
    assert_equal(read, adata[masks["obs"], masks["var"]])


@pytest.mark.parametrize(
    ("filter", "match"),
    [
        pytest.param(A.obs["obs_cat"], r"not a boolean column", id="not-bool"),
        pytest.param(A.obs.index, r"only `obs`/`var` columns", id="index"),
        pytest.param(A.obsm["array"][:, 0], r"only `obs`/`var` columns", id="obsm"),
    ],
)
def test_read_filter_invalid(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    filter: AdRef,
    match: str,
) -> None:
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    with pytest.raises(ValueError, match=match):
        getattr(ad.io, f"read_{diskfmt}")(diskfmt_store, filter=filter)
