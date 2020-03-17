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


#######################################
# Dealing with uns adj matrices
#######################################


def test_get_uns_neighbors_deprecated(adata):
    n = adata.shape[0]
    mtx = sparse.random(n, n, density=0.3, format="csr")
    adata.obsp["connectivities"] = mtx
    adata.uns["neighbors"] = {}

    with pytest.warns(FutureWarning):
        from_uns = adata.uns["neighbors"]["connectivities"]

    assert_equal(from_uns, mtx)

    with pytest.warns(None) as rec:
        v = adata[: n // 2]
        assert not rec

    with pytest.warns(FutureWarning):
        from_uns_v = v.uns["neighbors"]["connectivities"]

    assert_equal(from_uns_v, v.obsp["connectivities"])


def test_set_uns_neighbors_deprecated(adata):
    n = adata.shape[0]
    mtx = sparse.random(n, n, format="csr")
    adata.uns["neighbors"] = {}

    with pytest.warns(FutureWarning):
        adata.uns["neighbors"]["connectivities"] = sparse.random(n, n, format="csr")

    assert_equal(adata.obsp["connectivities"], mtx)
    with pytest.warns(FutureWarning):
        assert_equal(adata.uns["neighbors"]["connectivities"], mtx)

    # Make sure that we can write to uns normally:
    adata.uns["new_key"] = 100
    assert adata.uns["new_key"] == 100


def test_slice_uns_sparse_deprecated(adata):
    n = adata.shape[0]
    mtx = sparse.random(n, n, density=0.2, format="csr")
    adata.uns["sparse_mtx"] = mtx

    with pytest.warns(FutureWarning):
        v = adata[: n // 2]

    assert_equal(adata.uns["sparse_mtx"], mtx)
    assert_equal(v.uns["sparse_mtx"], mtx[: n // 2, :][:, n // 2])


@pytest.fixture
def adata_neighbors():
    return ad.AnnData(
        X=sparse.random(100, 200, format="csr"),
        obsp=dict(
            distances=sparse.random(100, 100, format="csr"),
            connectivities=sparse.random(100, 100, format="csr"),
        ),
        uns={"neighbors": {"params": {"method": "umap", "n_neighbors": 10}}},
    )


def test_deprecated_neighbors_get_mtx(adata_neighbors):
    """Test getting neighbor matrices from adata.uns"""
    adata = adata_neighbors

    with pytest.warns(FutureWarning):
        assert_equal(adata.obsp["distances"], adata.uns["neighbors"]["distances"])
    with pytest.warns(FutureWarning):
        assert_equal(
            adata.obsp["connectivities"], adata.uns["neighbors"]["connectivities"]
        )


def test_deprecated_neighbors_get_other(adata_neighbors):
    """Test getting other fields from adata.uns"""
    adata = adata_neighbors

    # This shouldn't throw a warning
    with pytest.warns(None) as rec:
        assert adata.uns["neighbors"]["params"] == {"method": "umap", "n_neighbors": 10}
        assert not rec


def test_deprecated_neighbors_set_other(adata_neighbors):
    adata = adata_neighbors

    # This shouldn't throw a warning
    with pytest.warns(None) as rec:
        adata.uns["neighbors"]["new_key"] = 10
        assert adata.uns["neighbors"]["new_key"] == 10
        # Test nested
        adata.uns["neighbors"]["params"]["new_param"] = 100
        assert adata.uns["neighbors"]["params"]["new_param"] == 100
        assert adata.uns["neighbors"]["params"] == {
            "method": "umap",
            "n_neighbors": 10,
            "new_param": 100,
        }

        assert not rec
