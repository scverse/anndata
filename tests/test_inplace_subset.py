from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from anndata.tests.helpers import (
    as_dense_dask_array,
    assert_equal,
    gen_adata,
)
from anndata.utils import asarray


@pytest.fixture(
    params=[
        np.array,
        sparse.csr_matrix,
        sparse.csc_matrix,
        sparse.csr_array,
        sparse.csc_array,
        as_dense_dask_array,
    ],
    ids=[
        "np_array",
        "scipy_csr",
        "scipy_csc",
        "scipy_csr_array",
        "scipy_csc_array",
        "dask_array",
    ],
)
def matrix_type(request):
    return request.param


def subset_dim(adata, *, obs=slice(None), var=slice(None)):
    # Should probably get used for test_inplace_subset_var and test_inplace_subset_obs
    from anndata._core.index import _subset

    return _subset(adata, (obs, var))


# TODO: Test values of .uns
def test_inplace_subset_var(matrix_type, subset_func):
    orig = gen_adata((30, 30), X_type=matrix_type)
    subset_idx = subset_func(orig.var_names)

    modified = orig.copy()
    from_view = orig[:, subset_idx].copy()
    modified._inplace_subset_var(subset_idx)

    assert_equal(asarray(from_view.X), asarray(modified.X), exact=True)
    assert_equal(from_view.obs, modified.obs, exact=True)
    assert_equal(from_view.var, modified.var, exact=True)
    for k in from_view.obsm:
        assert_equal(from_view.obsm[k], modified.obsm[k], exact=True)
        assert_equal(orig.obsm[k], modified.obsm[k], exact=True)
    for k in from_view.varm:
        assert_equal(from_view.varm[k], modified.varm[k], exact=True)
    for k in from_view.layers:
        assert_equal(from_view.layers[k], modified.layers[k], exact=True)


def test_inplace_subset_obs(matrix_type, subset_func):
    orig = gen_adata((30, 30), X_type=matrix_type)
    subset_idx = subset_func(orig.obs_names)

    modified = orig.copy()
    from_view = orig[subset_idx, :].copy()
    modified._inplace_subset_obs(subset_idx)

    assert_equal(asarray(from_view.X), asarray(modified.X), exact=True)
    assert_equal(from_view.obs, modified.obs, exact=True)
    assert_equal(from_view.var, modified.var, exact=True)
    for k in from_view.obsm:
        assert_equal(from_view.obsm[k], modified.obsm[k], exact=True)
    for k in from_view.varm:
        assert_equal(from_view.varm[k], modified.varm[k], exact=True)
        assert_equal(orig.varm[k], modified.varm[k], exact=True)
    for k in from_view.layers:
        assert_equal(from_view.layers[k], modified.layers[k], exact=True)


@pytest.mark.parametrize("dim", ["obs", "var"])
def test_inplace_subset_no_X(subset_func, dim):
    orig = gen_adata((30, 30))
    del orig.X

    subset_idx = subset_func(getattr(orig, f"{dim}_names"))

    modified = orig.copy()
    from_view = subset_dim(orig, **{dim: subset_idx}).copy()
    getattr(modified, f"_inplace_subset_{dim}")(subset_idx)

    assert_equal(modified, from_view, exact=True)
