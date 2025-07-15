from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import assert_equal

ARCHIVE_PTH = Path(__file__).parent / "data/archives"


@pytest.fixture(params=list(ARCHIVE_PTH.glob("v*")), ids=lambda x: x.name)
def archive_dir(request):
    return request.param


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


def test_backwards_compat_files(archive_dir) -> None:
    with pytest.warns(ad.OldFormatWarning):
        from_h5ad = ad.read_h5ad(archive_dir / "adata.h5ad")
    with pytest.warns(ad.OldFormatWarning):
        from_zarr = ad.read_zarr(archive_dir / "adata.zarr.zip")

    assert_equal(from_h5ad, from_zarr, exact=True)


def test_clean_uns_backwards_compat(tmp_path, diskfmt):
    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)

    orig = ad.AnnData(
        sparse.csr_matrix((3, 5), dtype="float32"),
        obs=pd.DataFrame(
            {"a": pd.Categorical(list("aab")), "b": [1, 2, 3]},
            index=[f"cell_{i}" for i in range(3)],
        ),
        uns={
            "a_categories": "some string",
            "b_categories": "another string",
        },
    )

    write(orig, pth)
    from_disk = read(pth)

    assert_equal(orig, from_disk)
