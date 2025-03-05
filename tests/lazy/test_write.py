from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from anndata import AnnData

pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.mark.parametrize("fmt", ["zarr", "h5ad", "loom", "csvs"])
def test_write_error(
    adata_remote: AnnData,
    tmp_path: Path,
    fmt: Literal["zarr", "h5ad", "loom", "csvs"],
):
    path = tmp_path / f"adata.{fmt}"
    with pytest.raises(
        NotImplementedError,
        match=r"Writing AnnData objects with a Dataset2D not supported. Please use `ds.to_memory`",
    ):
        getattr(adata_remote, f"write_{fmt}")(path)
