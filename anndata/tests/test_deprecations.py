"""This file contains tests for deprecated functions.

This includes correct behaviour as well as throwing warnings.
"""
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from anndata import AnnData


@pytest.fixture
def adata():
    adata = AnnData(
        X=csr_matrix([[0, 2, 3], [0, 5, 6]]),
        obs=dict(obs_names=['s1', 's2'], anno1=['c1', 'c2']),
        var=dict(var_names=['a', 'b', 'c']),
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


# TODO: Why doesn't this mark work?
# @pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_obsvar_array(adata):
    with pytest.warns(DeprecationWarning):  # Just to hide warnings
        assert np.allclose(adata._get_obs_array("a"), adata.obs_vector("a"))
        assert np.allclose(
            adata._get_obs_array("a", layer="x2"),
            adata.obs_vector("a", layer="x2"),
        )
        assert np.allclose(
            adata._get_obs_array("a", use_raw=True), adata.raw.obs_vector("a")
        )
        assert np.allclose(adata._get_var_array("s1"), adata.var_vector("s1"))
        assert np.allclose(
            adata._get_var_array("s1", layer="x2"),
            adata.var_vector("s1", layer="x2"),
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
        # This time it shouldn't throw a warning
        if "anndata" in r.filename:
            assert r.category is not FutureWarning
