"""Tests that AnnData enforces a 2-D ``X`` and ``layers`` at the IO boundary."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr

import anndata as ad
from anndata.io import write_elem

if TYPE_CHECKING:
    from pathlib import Path


MSG_PATTERN = r"must be 2-dimensional"

# Where to put the non-2-D array on disk for the read tests, and which attribute
# to make non-2-D in memory for the write tests.
WhichAttr = Literal["X", "layers"]

DISK_FORMATS = [
    pytest.param("h5ad", id="h5ad"),
    pytest.param("zarr", id="zarr"),
]
WHICH_ATTRS = [
    pytest.param("X", id="X"),
    pytest.param("layers", id="layers"),
]


@pytest.fixture
def arr2d() -> np.ndarray:
    return np.zeros((3, 4))


@pytest.fixture
def arr3d() -> np.ndarray:
    return np.arange(3 * 4 * 5).reshape((3, 4, 5))


def test_in_memory_non_2d_x_and_layers_allowed(
    arr2d: np.ndarray, arr3d: np.ndarray
) -> None:
    """All in-memory mutation paths accept non-2-D arrays for ``X`` / ``layers``.

    Only the IO layer enforces the 2-D spec, so users can freely hold
    higher-dimensional intermediates in memory.
    """
    # via constructor
    adata = ad.AnnData(X=arr3d)
    assert adata.X.shape == arr3d.shape

    adata = ad.AnnData(X=arr2d, layers={"L": arr3d})
    assert adata.layers["L"].shape == arr3d.shape

    # via attribute setter
    adata = ad.AnnData(X=arr2d)
    adata.X = arr3d
    assert adata.X.shape == arr3d.shape

    # via item assignment
    adata = ad.AnnData(X=arr2d)
    adata.layers["L"] = arr3d
    assert adata.layers["L"].shape == arr3d.shape

    # via full mapping reassignment
    adata.layers = {"M": arr3d}
    assert adata.layers["M"].shape == arr3d.shape


def _make_non_conforming(
    tmp_path: Path,
    arr2d: np.ndarray,
    arr3d: np.ndarray,
    *,
    diskfmt: Literal["h5ad", "zarr"],
    which: WhichAttr,
) -> Path:
    """Build a syntactically-valid AnnData file where exactly one of ``X`` /
    ``layers["L"]`` is non-2-D.

    We can't go through ``adata.write_*`` because the writer rejects non-2-D
    payloads; instead we lay out the on-disk structure directly via the
    low-level ``write_elem`` API.
    """
    n_obs, n_var = arr3d.shape[:2]
    obs = pd.DataFrame(index=[f"o{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"v{i}" for i in range(n_var)])
    if diskfmt == "h5ad":
        pth = tmp_path / "non_conforming.h5ad"
        # Start from a conforming file so obs/var/encoding are well-formed,
        # then surgically replace just the part we want to be 3-D.
        ad.AnnData(X=arr2d, obs=obs, var=var).write_h5ad(pth)
        with h5py.File(pth, "a") as f:
            if which == "X":
                del f["X"]
                write_elem(f, "X", arr3d)
            else:  # "layers"
                write_elem(f, "layers", {"L": arr3d})
        return pth
    if diskfmt == "zarr":
        pth = tmp_path / "non_conforming.zarr"
        f = zarr.open_group(pth, mode="w")
        f.attrs.setdefault("encoding-type", "anndata")
        f.attrs.setdefault("encoding-version", "0.1.0")
        write_elem(f, "X", arr3d if which == "X" else arr2d)
        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        if which == "layers":
            write_elem(f, "layers", {"L": arr3d})
        # Consolidated metadata is not yet part of the Zarr v3 spec; zarr-python
        # warns about that on every call but it is orthogonal to what we test
        # here, so we silence just this one helper warning.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Consolidated metadata is",
                category=UserWarning,
            )
            zarr.consolidate_metadata(f.store)
        return pth
    msg = f"Unknown diskfmt: {diskfmt!r}"
    raise ValueError(msg)


def _read(pth: Path, diskfmt: Literal["h5ad", "zarr"]) -> ad.AnnData:
    return ad.read_h5ad(pth) if diskfmt == "h5ad" else ad.read_zarr(pth)


def _write(adata: ad.AnnData, pth: Path, diskfmt: Literal["h5ad", "zarr"]) -> None:
    if diskfmt == "h5ad":
        adata.write_h5ad(pth)
    else:
        adata.write_zarr(pth)


@pytest.mark.parametrize("diskfmt", DISK_FORMATS)
@pytest.mark.parametrize("which", WHICH_ATTRS)
def test_read_with_non_2d_warns(
    tmp_path: Path,
    arr2d: np.ndarray,
    arr3d: np.ndarray,
    diskfmt: Literal["h5ad", "zarr"],
    which: WhichAttr,
) -> None:
    """Reading a file with a non-2-D ``X`` *or* a non-2-D layer warns but succeeds.

    The 2-D attribute on the same file must still load with its
    original (2-D) shape and must not trigger the spec warning.
    """
    pth = _make_non_conforming(tmp_path, arr2d, arr3d, diskfmt=diskfmt, which=which)
    with pytest.warns(UserWarning, match=MSG_PATTERN) as record:
        adata = _read(pth, diskfmt)
    if which == "X":
        assert adata.X.shape == arr3d.shape
        # The 2-D layers entry doesn't exist on this file, so nothing more to check.
        # Make sure the warning was specifically about X.
        assert any("X must be" in str(w.message) for w in record)
        assert not any("Layer" in str(w.message) for w in record)
    else:  # "layers"
        assert adata.X.shape == arr2d.shape  # untouched
        assert adata.layers["L"].shape == arr3d.shape
        assert any("Layer 'L'" in str(w.message) for w in record)
        assert not any("X must be" in str(w.message) for w in record)


@pytest.mark.parametrize("diskfmt", DISK_FORMATS)
@pytest.mark.parametrize("which", WHICH_ATTRS)
def test_write_with_non_2d_raises(
    tmp_path: Path,
    arr2d: np.ndarray,
    arr3d: np.ndarray,
    diskfmt: Literal["h5ad", "zarr"],
    which: WhichAttr,
) -> None:
    """Writing an AnnData with a non-2-D ``X`` or layer raises a ``ValueError``."""
    adata = ad.AnnData(X=arr2d)
    if which == "X":
        adata.X = arr3d  # allowed in memory
    else:
        adata.layers["L"] = arr3d  # allowed in memory
    out = tmp_path / f"out.{diskfmt}"
    with pytest.raises(ValueError, match=MSG_PATTERN):
        _write(adata, out, diskfmt)
