"""Tests for backing by just sticking zarr/h5py objects into AnnData."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import h5py
import numpy as np
import pytest
import zarr

from anndata import AnnData
from anndata._io.specs import write_elem
from anndata.tests.helpers import assert_equal

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.mark.parametrize("assign", ["init", "assign"])
@pytest.mark.parametrize("attr", ["X", "obsm", "varm", "layers"])
def test_backed_init(
    tmp_path: Path,
    diskfmt: Literal["h5ad", "zarr"],
    assign: Literal["init", "assign"],
    attr: Literal["X", "obsm", "varm", "layers"],
):
    path = tmp_path / f"test.{diskfmt}"
    if diskfmt == "zarr":
        f = zarr.open_group(path, mode="a")
    elif diskfmt == "h5ad":
        f = h5py.File(path, mode="a")
    else:
        pytest.fail(f"Unexpected diskfmt: {diskfmt}")

    x = np.random.randn(10, 10)
    write_elem(f, "a", x)

    # initialize (and if applicable, assign)
    if assign == "init":
        kw = dict(X=f["a"]) if attr == "X" else {attr: dict(a=f["a"]), "shape": x.shape}
        adata_backed = AnnData(**kw)
    elif assign == "assign":
        adata_backed = AnnData(shape=x.shape)
        if attr == "X":
            adata_backed.X = f["a"]
        else:
            getattr(adata_backed, attr)["a"] = f["a"]
    else:
        pytest.fail(f"Unexpected assign: {assign}")

    # check equality
    if attr == "X":
        # TODO: should that be inverted, e.g. when the Dataset’s path matches the backed mode path?
        assert not adata_backed.isbacked
        backed_array = adata_backed.X
    else:
        backed_array = getattr(adata_backed, attr)["a"]
    assert isinstance(backed_array, zarr.Array if diskfmt == "zarr" else h5py.Dataset)
    assert_equal(backed_array, x)

    # check that there’s no error deleting it either
    if attr == "X":
        del adata_backed.X
    else:
        del getattr(adata_backed, attr)["a"]
