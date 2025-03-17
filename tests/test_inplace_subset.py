from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from fast_array_utils.conv import to_dense

from anndata.tests.helpers import assert_equal, gen_adata
from testing.fast_array_utils import SUPPORTED_TYPES, Flags

if TYPE_CHECKING:
    from testing.fast_array_utils import ArrayType


SPARSE_DASK = {
    at for at in SUPPORTED_TYPES if at.flags & Flags.Sparse and at.flags & Flags.Dask
}


def subset_dim(adata, *, obs=slice(None), var=slice(None)):
    # Should probably get used for test_inplace_subset_var and test_inplace_subset_obs
    from anndata._core.index import _subset

    return _subset(adata, (obs, var))


@pytest.mark.array_type(skip={Flags.Gpu | Flags.Disk, *SPARSE_DASK})
# TODO: Test values of .uns
def test_inplace_subset_var(array_type: ArrayType, subset_func) -> None:
    orig = gen_adata((30, 30), X_type=array_type)
    subset_idx = subset_func(orig.var_names)

    modified = orig.copy()
    from_view = orig[:, subset_idx].copy()
    modified._inplace_subset_var(subset_idx)

    assert_equal(to_dense(from_view.X), to_dense(modified.X), exact=True)
    assert_equal(from_view.obs, modified.obs, exact=True)
    assert_equal(from_view.var, modified.var, exact=True)
    for k in from_view.obsm:
        assert_equal(from_view.obsm[k], modified.obsm[k], exact=True)
        assert_equal(orig.obsm[k], modified.obsm[k], exact=True)
    for k in from_view.varm:
        assert_equal(from_view.varm[k], modified.varm[k], exact=True)
    for k in from_view.layers:
        assert_equal(from_view.layers[k], modified.layers[k], exact=True)


@pytest.mark.array_type(skip={Flags.Gpu | Flags.Disk, *SPARSE_DASK})
def test_inplace_subset_obs(array_type: ArrayType, subset_func) -> None:
    orig = gen_adata((30, 30), X_type=array_type)
    subset_idx = subset_func(orig.obs_names)

    modified = orig.copy()
    from_view = orig[subset_idx, :].copy()
    modified._inplace_subset_obs(subset_idx)

    assert_equal(to_dense(from_view.X), to_dense(modified.X), exact=True)
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
