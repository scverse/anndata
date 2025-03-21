from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pytest

from anndata import AnnData
from anndata.experimental.backed._io import read_lazy

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.mark.parametrize("fmt", ["zarr", "h5ad", "loom", "csvs"])
@pytest.mark.parametrize("key", ["obs", "var", "obsm", "varm"])
def test_write_error(
    tmp_path: Path,
    fmt: Literal["zarr", "h5ad", "loom", "csvs"],
    key: Literal["obs", "var", "obsm", "varm"],
):
    path = tmp_path / "adata.h5ad"
    X = np.random.random((4, 4))
    adata = AnnData(X=X)
    if key.endswith("m"):
        elem = {"df": getattr(adata, key[:-1])}
        setattr(adata, key, elem)
    adata.write_h5ad(path)
    adata_lazy = read_lazy(path)
    if key.endswith("m"):
        adata_lazy.obs = adata_lazy.obs.to_memory()
        adata_lazy.obs = adata_lazy.var.to_memory()
    noop_path = tmp_path / f"adata_noop.{fmt}"
    with pytest.raises(
        NotImplementedError,
        match=r"Writing AnnData objects with a Dataset2D not supported yet. Please use `ds.to_memory`",
    ):
        getattr(adata_lazy, f"write_{fmt}")(noop_path)
    assert not noop_path.exists(), (
        "Found a directory at the path at which no data should have been written"
    )
