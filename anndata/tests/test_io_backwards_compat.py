from pathlib import Path

import pytest

import anndata as ad
from anndata.tests.helpers import assert_equal

ARCHIVE_PTH = Path(__file__).parent / "data/archives"


@pytest.fixture(params=list(ARCHIVE_PTH.glob("v*")), ids=lambda x: x.name)
def archive_dir(request):
    return request.param


def test_backwards_compat_files(archive_dir):
    with pytest.warns(ad.OldFormatWarning):
        from_h5ad = ad.read_h5ad(archive_dir / "adata.h5ad")
    with pytest.warns(ad.OldFormatWarning):
        from_zarr = ad.read_zarr(archive_dir / "adata.zarr.zip")

    assert_equal(from_h5ad, from_zarr, exact=True)
