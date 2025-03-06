from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import pandas as pd
import pytest

from anndata import AnnData
from anndata.experimental.backed._io import read_lazy

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.mark.parametrize("fmt", ["zarr", "h5ad", "loom", "csvs"])
def test_write_error(
    tmp_path: Path,
    fmt: Literal["zarr", "h5ad", "loom", "csvs"],
):
    path = tmp_path / "adata.h5ad"
    AnnData(obs=pd.DataFrame({"column": list("abcd")})).write_h5ad(path)
    adata_remote = read_lazy(path)
    noop_path = tmp_path / f"adata_noop.{fmt}"
    with pytest.raises(
        NotImplementedError,
        match=r"Writing AnnData objects with a Dataset2D not supported. Please use `ds.to_memory`",
    ):
        getattr(adata_remote, f"write_{fmt}")(noop_path)
    assert not noop_path.exists(), (
        "Found a directory at the path at which no data should have been written"
    )
