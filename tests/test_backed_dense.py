"""Tests for backing by just sticking zarr/h5py objects into AnnData."""

from __future__ import annotations

from typing import TYPE_CHECKING

import h5py
import numpy as np
import pytest
import zarr

from anndata import AnnData
from anndata._io.zarr import open_write_group
from anndata.io import write_elem
from anndata.tests.helpers import assert_equal

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


@pytest.fixture
def file(tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]) -> h5py.File | zarr.Group:
    path = tmp_path / f"test.{diskfmt}"
    if diskfmt == "zarr":
        return open_write_group(path, mode="a")
    if diskfmt == "h5ad":
        return h5py.File(path, "a")
    pytest.fail(f"Unknown diskfmt: {diskfmt}")


@pytest.mark.parametrize("assign", ["init", "assign"])
@pytest.mark.parametrize("attr", ["X", "obsm", "varm", "layers"])
def test_create_delete(
    diskfmt: Literal["h5ad", "zarr"],
    file: h5py.File | zarr.Group,
    assign: Literal["init", "assign"],
    attr: Literal["X", "obsm", "varm", "layers"],
):
    x = np.random.randn(10, 10)
    write_elem(file, "a", x)

    # initialize (and if applicable, assign)
    if assign == "init":
        kw = (
            dict(X=file["a"])
            if attr == "X"
            else {attr: dict(a=file["a"]), "shape": x.shape}
        )
        adata = AnnData(**kw)
    elif assign == "assign":
        adata = AnnData(shape=x.shape)
        if attr == "X":
            adata.X = file["a"]
        else:
            getattr(adata, attr)["a"] = file["a"]
    else:
        pytest.fail(f"Unexpected assign: {assign}")

    # check equality
    if attr == "X":
        # TODO: should that be inverted, e.g. when the Dataset’s path matches the backed mode path?
        assert not adata.isbacked
        backed_array = adata.X
    else:
        backed_array = getattr(adata, attr)["a"]
    assert isinstance(backed_array, zarr.Array if diskfmt == "zarr" else h5py.Dataset)
    assert_equal(backed_array, x)

    # check that there’s no error deleting it either
    if attr == "X":
        del adata.X
    else:
        del getattr(adata, attr)["a"]


def test_assign_x_subset(file: h5py.File | zarr.Group):
    x = np.ones((10, 10))
    write_elem(file, "a", x)

    adata = AnnData(file["a"])

    view = adata[3:7, 6:8]
    view.X = np.zeros((4, 2))

    expected = x.copy()
    expected[3:7, 6:8] = np.zeros((4, 2))
    assert_equal(adata.X, expected)
