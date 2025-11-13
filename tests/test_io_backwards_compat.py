from __future__ import annotations

from pathlib import Path
from subprocess import run
from typing import TYPE_CHECKING

import pandas as pd
import pytest
import zarr
import zarr.storage
from scipy import sparse

import anndata as ad
from anndata.compat import is_zarr_v2
from anndata.tests.helpers import assert_equal

if TYPE_CHECKING:
    from typing import Literal

ARCHIVE_PTH = Path(__file__).parent / "data/archives"


@pytest.fixture(params=list(ARCHIVE_PTH.glob("v*")), ids=lambda x: x.name)
def archive_dir(request: pytest.FixtureRequest) -> Path:
    assert isinstance(request.param, Path)
    if request.param.name == "v0.5.0":
        reason = "v0.5.0 has no zarr and very different format"
        request.applymarker(pytest.mark.xfail(reason=reason))
    return request.param


def read_archive(
    archive_dir: Path, format: Literal["h5ad", "zarr"]
) -> tuple[ad.AnnData, Path]:
    if format == "h5ad":
        path = archive_dir / "adata.h5ad"
        return ad.read_h5ad(path), path
    if format == "zarr":
        path = archive_dir / "adata.zarr.zip"
        store = path if is_zarr_v2() else zarr.storage.ZipStore(path)
        return ad.read_zarr(store), path
    pytest.fail(f"Unknown format: {format}")


@pytest.mark.filterwarnings("ignore::anndata.OldFormatWarning")
def test_backwards_compat_files(archive_dir: Path) -> None:
    from_h5ad, _ = read_archive(archive_dir, "h5ad")
    from_zarr, _ = read_archive(archive_dir, "zarr")
    assert_equal(from_h5ad, from_zarr, exact=True)


def test_no_diff(
    request: pytest.FixtureRequest, tmp_path: Path, archive_dir: Path
) -> None:
    if pd.options.future.infer_string:
        reason = "even old files get read as `string` dtype if this is active"
        request.applymarker(pytest.mark.xfail(reason=reason))
    if archive_dir.name in {"v0.7.8", "v0.7.0"}:
        pytest.skip("DataFrame encoding changed between 0.7 and now")
    adata, in_path = read_archive(archive_dir, "h5ad")
    adata.write_h5ad(out_path := tmp_path / "adata.h5ad")
    diff_proc = run(
        ["h5diff", "-c", in_path, out_path], check=False, capture_output=True, text=True
    )
    assert diff_proc.returncode == 0, diff_proc.stdout


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
