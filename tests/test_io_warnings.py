from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pytest
from zarr.storage import MemoryStore

import anndata as ad
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, gen_adata

HERE = Path(__file__).parent
DATA_DIR = HERE / "data"


def test_old_format_warning_thrown() -> None:
    def msg_re(entry: str) -> str:
        return re.escape(
            f"Moving element from .uns['neighbors'][{entry!r}] to .obsp[{entry!r}]."
        )

    pth = DATA_DIR / "archives/v0.5.0/adata.h5ad"
    warnings.simplefilter("default", FutureWarning)
    with (
        pytest.warns(FutureWarning, match=msg_re("distances")),
        pytest.warns(FutureWarning, match=msg_re("connectivities")),
        pytest.warns(ad.OldFormatWarning),
    ):
        ad.read_h5ad(pth)


def test_old_format_warning_not_thrown(tmp_path: Path) -> None:
    pth = tmp_path / "current.h5ad"
    adata = gen_adata((20, 10), **GEN_ADATA_NO_XARRAY_ARGS)
    adata.write_h5ad(pth)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always", ad.OldFormatWarning)

        ad.read_h5ad(pth)

    if len(record) != 0:
        msg_content = "\n".join([
            f"\t{w.category.__name__}('{w.message}')" for w in record
        ])
        pytest.fail(
            f"Warnings were thrown when they shouldn't be. Got:\n\n{msg_content}"
        )


def test_zarr_write_format_deprecated():
    store = MemoryStore()
    # use a context so we don't trigger the warning when the setting resets
    with (
        pytest.warns(DeprecationWarning, match="zarr v3 will be the default"),
        ad.settings.override(zarr_write_format=2),
        pytest.warns(DeprecationWarning, match="zarr v3 will become the only option"),
    ):
        ad.AnnData(np.ones((10, 10))).write_zarr(store)
