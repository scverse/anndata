"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from typing import TYPE_CHECKING, get_args

import h5py
import numpy as np
import pytest
import zarr
from packaging.version import Version
from zarr.storage import MemoryStore

import anndata as ad
from anndata.io import write_elem

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


def _shape_on_disk(
    store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"], key: str
) -> tuple[int, ...]:
    if diskfmt == "zarr":
        elem = zarr.open_group(store, mode="r")[key]
        assert isinstance(elem, zarr.Array)
        return elem.shape
    with h5py.File(store, "r") as f:
        return f[key].shape


def test_write_compat_compares_as_versions() -> None:
    assert Version("0.9") < ad.WriteCompat.V0_12
    assert Version("0.12") == ad.WriteCompat.V0_12


def test_literal_matches_enum() -> None:
    assert get_args(ad.WriteCompatStr.__value__) == tuple(
        str(c) for c in ad.WriteCompat
    )


def test_default_is_alias() -> None:
    assert ad.WriteCompat.DEFAULT in ad.WriteCompat
    assert "DEFAULT" not in {c.name for c in ad.WriteCompat}


def test_takes_str() -> None:
    assert ad.WriteCompat("0.12") is ad.WriteCompat.V0_12
    # this anndata version can’t write anything older than 0.12
    with pytest.raises(ValueError, match=r"0\.11"):
        ad.AnnData().unwriteable(compat="0.11")  # type: ignore[arg-type]


def test_non_2d_layer(
    *,
    diskfmt_store: Path | MemoryStore,
    diskfmt: Literal["h5ad", "zarr"],
    arr2d: np.ndarray,
    arr3d: np.ndarray,
) -> None:
    """`compat="0.12"` writes 3-D layers, which 0.13 banned."""
    adata = ad.AnnData(X=arr2d)
    with pytest.warns(UserWarning, match=r"must be 2-dimensional"):
        adata.layers["L"] = arr3d

    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.12")

    assert _shape_on_disk(diskfmt_store, diskfmt, "layers/L") == arr3d.shape


def test_non_2d_layer_hint(
    diskfmt_store: Path | MemoryStore,
    diskfmt: Literal["h5ad", "zarr"],
    arr2d: np.ndarray,
    arr3d: np.ndarray,
) -> None:
    """The default profile refuses and points at the escape hatch."""
    adata = ad.AnnData(X=arr2d)
    with pytest.warns(UserWarning, match=r"must be 2-dimensional"):
        adata.layers["L"] = arr3d
    with pytest.raises(ValueError, match=r'Pass `compat="0.12"`'):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)


def test_unwriteable(arr2d: np.ndarray, arr3d: np.ndarray) -> None:
    adata = ad.AnnData(X=arr2d)
    with pytest.warns(UserWarning, match=r"must be 2-dimensional"):
        adata.layers["L"] = arr3d
    assert adata.unwriteable()
    assert not adata.unwriteable(compat="0.12")


def _write_h5(
    writer: Literal["write_h5ad", "write_elem"], path: Path, **kwargs
) -> None:
    match writer:
        case "write_h5ad":
            adata = ad.AnnData(shape=(3, 4), uns={"bad/key": np.ones(3)})
            adata.write_h5ad(path, **kwargs)
        case "write_elem":
            with h5py.File(path, "w") as f:
                write_elem(f, "uns/bad/key", np.ones(3), **kwargs)
        case _:
            pytest.fail(f"Unknown writer: {writer}")


@pytest.mark.parametrize("writer", ["write_h5ad", "write_elem"])
def test_slash_key_h5_compat(
    tmp_path: Path, writer: Literal["write_h5ad", "write_elem"]
) -> None:
    with pytest.warns(FutureWarning, match=r"Forward slashes"):
        _write_h5(writer, tmp_path / "allowed.h5ad", compat="0.12")
    with h5py.File(tmp_path / "allowed.h5ad") as f:
        assert f["uns/bad/key"].shape == (3,)


@pytest.mark.parametrize("writer", ["write_h5ad", "write_elem"])
def test_slash_key_h5_error(
    tmp_path: Path, writer: Literal["write_h5ad", "write_elem"]
) -> None:
    with pytest.raises(ValueError, match=r'Pass `compat="0.12"`'):
        _write_h5(writer, tmp_path / "banned.h5ad")


def test_slash_key_zarr_error() -> None:
    adata = ad.AnnData(shape=(3, 4), uns={"bad/key": np.ones(3)})
    with pytest.raises(ValueError, match=r"Forward slashes"):
        adata.write_zarr(MemoryStore(), compat="0.12")


@pytest.mark.parametrize(
    ("compat", "disallow"), [("0.12", True), ("0.13", False)], ids=["disallow", "allow"]
)
@pytest.mark.filterwarnings(
    "ignore:This will be removed in 0.15.*pass `compat=`:DeprecationWarning"
)
def test_deprecated_setting_wins(
    *, tmp_path: Path, disallow: bool, compat: ad.WriteCompatStr
) -> None:
    adata = ad.AnnData(shape=(3, 4), uns={"bad/key": np.ones(3)})
    ctx = (
        pytest.raises(ValueError, match=r"Forward slashes")
        if disallow
        else pytest.warns(FutureWarning, match=r"Forward slashes")
    )
    with ad.settings.override(disallow_forward_slash_in_h5ad=disallow), ctx:
        adata.write_h5ad(tmp_path / "test.h5ad", compat=compat)
