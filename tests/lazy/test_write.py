from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from anndata import AnnData


@pytest.mark.parametrize("fmt", ["zarr", "h5ad", "loom", "csvs"])
def test_write_error(
    adata_remote: AnnData,
    tmp_path: Path,
    fmt: Literal["zarr", "h5ad", "loom", "csvs"],
):
    path = tmp_path / f"adata.{fmt}"
    with pytest.raises(
        NotImplementedError,
        match=r"Writing AnnData objects with a Dataset2D not supported",
    ):
        getattr(adata_remote, f"write_{fmt}")(path)
