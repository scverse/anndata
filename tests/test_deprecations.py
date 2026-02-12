"""\
This file contains tests for deprecated functions.

This includes correct behaviour as well as throwing warnings.
"""

from __future__ import annotations

import warnings

import h5py
import numpy as np
import pytest
from scipy import sparse

import anndata.experimental
from anndata import AnnData
from anndata.tests.helpers import assert_equal


@pytest.fixture
def adata() -> AnnData:
    adata = AnnData(
        X=sparse.csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
        obs=dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        var=dict(var_names=["a", "b", "c"]),
    )
    adata.raw = adata.copy()
    adata.layers["x2"] = adata.X * 2
    adata.var["anno2"] = ["p1", "p2", "p3"]
    adata.X = adata.X / 2
    return adata


def test_get_obsvar_array_warn(adata: AnnData) -> None:
    with pytest.warns(FutureWarning):
        adata.obs_vector("a")
    with pytest.warns(FutureWarning):
        adata.var_vector("s1")


def test_deprecated_write_attribute(tmp_path):
    pth = tmp_path / "file.h5"
    A = np.random.randn(20, 10)
    from anndata._io.utils import read_attribute, write_attribute
    from anndata.io import read_elem

    with h5py.File(pth, "w") as f, pytest.warns(FutureWarning, match=r"write_elem"):
        write_attribute(f, "written_attribute", A)

    with h5py.File(pth, "r") as f:
        elem_A = read_elem(f["written_attribute"])
        with pytest.warns(FutureWarning, match=r"read_elem"):
            attribute_A = read_attribute(f["written_attribute"])

        assert_equal(elem_A, attribute_A)
        assert_equal(A, attribute_A)


@pytest.mark.parametrize(
    ("old_name", "new_name", "module"),
    (
        (old_name, new_name, module)
        for module in [anndata, anndata.experimental]
        for (old_name, new_name) in module._DEPRECATED.items()
    ),
)
def test_warn_on_import_with_redirect(old_name: str, new_name: str, module):
    with pytest.warns(FutureWarning, match=rf"Importing {old_name}.*is deprecated"):
        getattr(module, old_name)


def test_warn_on_deprecated__io_module():
    with pytest.warns(
        FutureWarning, match=r"Importing read_h5ad from `anndata._io` is deprecated"
    ):
        from anndata._io import read_h5ad  # noqa


@pytest.mark.parametrize("name", ["obs", "var", "obsm", "varm", "uns"])
def test_keys_function_warns(adata: AnnData, name) -> None:
    with pytest.warns(FutureWarning, match=rf"{name}_keys is deprecated"):
        getattr(adata, f"{name}_keys")()
