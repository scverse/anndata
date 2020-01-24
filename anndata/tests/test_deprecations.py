"""\
This file contains tests for deprecated functions.

This includes correct behaviour as well as throwing warnings.
"""
import h5py
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData

from anndata.tests.helpers import assert_equal


@pytest.fixture
def adata():
    adata = AnnData(
        X=sparse.csr_matrix([[0, 2, 3], [0, 5, 6]]),
        obs=dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        var=dict(var_names=["a", "b", "c"]),
    )
    adata.raw = adata
    adata.layers["x2"] = adata.X * 2
    adata.var["anno2"] = ["p1", "p2", "p3"]
    adata.X = adata.X / 2
    return adata


def test_get_obsvar_array_warn(adata):
    with pytest.warns(DeprecationWarning):
        adata._get_obs_array("a")
    with pytest.warns(DeprecationWarning):
        adata._get_var_array("s1")


# TODO: Why doesn’t this mark work?
# @pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_obsvar_array(adata):
    with pytest.warns(DeprecationWarning):  # Just to hide warnings
        assert np.allclose(adata._get_obs_array("a"), adata.obs_vector("a"))
        assert np.allclose(
            adata._get_obs_array("a", layer="x2"), adata.obs_vector("a", layer="x2"),
        )
        assert np.allclose(
            adata._get_obs_array("a", use_raw=True), adata.raw.obs_vector("a")
        )
        assert np.allclose(adata._get_var_array("s1"), adata.var_vector("s1"))
        assert np.allclose(
            adata._get_var_array("s1", layer="x2"), adata.var_vector("s1", layer="x2"),
        )
        assert np.allclose(
            adata._get_var_array("s1", use_raw=True), adata.raw.var_vector("s1")
        )


def test_obsvar_vector_Xlayer(adata):
    with pytest.warns(FutureWarning):
        adata.var_vector("s1", layer="X")
    with pytest.warns(FutureWarning):
        adata.obs_vector("a", layer="X")

    adata = adata.copy()
    adata.layers["X"] = adata.X * 3

    with pytest.warns(None) as records:
        adata.var_vector("s1", layer="X")
        adata.obs_vector("a", layer="X")

    for r in records:
        # This time it shouldn’t throw a warning
        if "anndata" in r.filename:
            assert r.category is not FutureWarning


def test_force_dense_deprecated(tmp_path):
    dense_pth = tmp_path / "dense.h5ad"
    adata = AnnData(X=sparse.random(10, 10, format="csr"))
    adata.raw = adata

    with pytest.warns(FutureWarning):
        adata.write_h5ad(dense_pth, force_dense=True)
    with h5py.File(dense_pth, "r") as f:
        assert isinstance(f["X"], h5py.Dataset)
        assert isinstance(f["raw/X"], h5py.Dataset)

    dense = ad.read_h5ad(dense_pth)

    assert isinstance(dense.X, np.ndarray)
    assert isinstance(dense.raw.X, np.ndarray)
    assert_equal(adata, dense)
