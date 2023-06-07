"""\
This file contains tests for deprecated functions.

This includes correct behaviour as well as throwing warnings.
"""
from pathlib import Path
import warnings

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
        X=sparse.csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
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

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        adata.var_vector("s1", layer="X")
        adata.obs_vector("a", layer="X")


# This should break in 0.9
def test_dtype_warning():
    # Tests a warning is thrown
    with pytest.warns(PendingDeprecationWarning):
        a = AnnData(np.ones((3, 3)), dtype=np.float32)
    assert a.X.dtype == np.float32

    # This shouldn't warn, shouldn't copy
    with warnings.catch_warnings(record=True) as record:
        b_X = np.ones((3, 3), dtype=np.float64)
        b = AnnData(b_X)
        assert not record
    assert b_X is b.X
    assert b.X.dtype == np.float64

    # Should warn, should copy
    with pytest.warns(PendingDeprecationWarning):
        c_X = np.ones((3, 3), dtype=np.float32)
        c = AnnData(c_X, dtype=np.float64)
        assert not record
    assert c_X is not c.X
    assert c.X.dtype == np.float64


def test_deprecated_write_attribute(tmp_path):
    pth = tmp_path / "file.h5"
    A = np.random.randn(20, 10)
    from anndata._io.utils import read_attribute, write_attribute
    from anndata._io.specs import read_elem

    with h5py.File(pth, "w") as f:
        with pytest.warns(DeprecationWarning, match="write_elem"):
            write_attribute(f, "written_attribute", A)

    with h5py.File(pth, "r") as f:
        elem_A = read_elem(f["written_attribute"])
        with pytest.warns(DeprecationWarning, match="read_elem"):
            attribute_A = read_attribute(f["written_attribute"])

        assert_equal(elem_A, attribute_A)
        assert_equal(A, attribute_A)
