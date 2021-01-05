import numpy as np
import pytest
from scipy import sparse

from anndata.tests.helpers import assert_equal, gen_adata, subset_func
from anndata.utils import asarray


@pytest.fixture(
    params=[np.array, sparse.csr_matrix, sparse.csc_matrix],
    ids=["np_array", "scipy_csr", "scipy_csc"],
)
def matrix_type(request):
    return request.param


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
        assert_equal(asarray(from_view.obsm[k]), asarray(modified.obsm[k]), exact=True)
        assert_equal(asarray(orig.obsm[k]), asarray(modified.obsm[k]), exact=True)
    for k in from_view.varm:
        assert_equal(asarray(from_view.varm[k]), asarray(modified.varm[k]), exact=True)
    for k in from_view.layers:
        assert_equal(
            asarray(from_view.layers[k]), asarray(modified.layers[k]), exact=True
        )


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
        assert_equal(asarray(from_view.obsm[k]), asarray(modified.obsm[k]), exact=True)
    for k in from_view.varm:
        assert_equal(asarray(from_view.varm[k]), asarray(modified.varm[k]), exact=True)
        assert_equal(asarray(orig.varm[k]), asarray(modified.varm[k]), exact=True)
    for k in from_view.layers:
        assert_equal(
            asarray(from_view.layers[k]), asarray(modified.layers[k]), exact=True
        )
