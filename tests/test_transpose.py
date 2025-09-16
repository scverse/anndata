from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import assert_equal, gen_adata, shares_memory


def test_transpose_orig():
    """
    Original test for transpose, should be covered by more thorough tests below, but
    keeping around just in case.
    """
    adata = gen_adata((5, 3))
    adata.varp = {f"varp_{k}": v for k, v in adata.varp.items()}
    adata1 = adata.T
    adata1.uns["test123"] = 1
    assert "test123" in adata.uns
    assert_equal(adata1.X.shape, (3, 5))
    assert_equal(adata1.obsp.keys(), adata.varp.keys())


def _add_raw(adata, *, var_subset=slice(None)):
    new = adata[:, var_subset].copy()
    new.raw = adata.copy()
    return new


# TODO: Cases to add:
# * Views
# * X is None should have the xfail marker removed
# * Backed
@pytest.fixture(
    params=[
        pytest.param(gen_adata((50, 20)), id="csr_X"),
        pytest.param(gen_adata((50, 20), sparse.csc_matrix), id="csc_X"),
        pytest.param(_add_raw(gen_adata((50, 20))), id="with_raw"),
        pytest.param(gen_adata((20, 10), X_type=None), id="None_X"),
    ]
)
def adata(request):
    return request.param


def test_transpose_doesnt_copy():
    adata = ad.AnnData(
        sparse.random(50, 20, format="csr"),
        layers={
            "sparse": sparse.random(50, 20, format="csc"),
            "dense": np.random.rand(50, 20),
        },
        obsm={
            "sparse": sparse.random(50, 10, format="csc"),
            "dense": np.random.rand(50, 10),
        },
        obsp={
            "sparse": sparse.random(50, 50, format="csc"),
            "dense": np.random.rand(50, 50),
        },
    )

    t = adata.T

    assert shares_memory(adata.X, t.X)
    for k in adata.obsm:
        assert shares_memory(adata.obsm[k], t.varm[k])
    for k in adata.obsp:
        assert shares_memory(adata.obsp[k], t.varp[k])
    for k in adata.layers:
        assert shares_memory(adata.layers[k], t.layers[k])


def test_transpose_removes_raw(adata):
    """
    Since Raw must have the same `obs_names` as AnnData, but does not have the same
    `var_names`, transpose doesn't really make sense for Raw. So it should just get
    deleted.
    """
    assert adata.T.raw is None


def test_transposed_contents(adata):
    t = adata.T

    if adata.X is not None:
        assert_equal(adata.X.T, t.X)
    else:
        assert adata.X is t.X is None

    assert_equal({k: v.T for k, v in adata.layers.items()}, dict(t.layers))
    assert_equal(adata.obs, t.var)
    assert_equal(adata.var, t.obs)
    assert_equal(dict(adata.obsm), dict(t.varm))
    assert_equal(dict(adata.varm), dict(t.obsm))
    assert_equal(dict(adata.obsp), dict(t.varp))
    assert_equal(dict(adata.varp), dict(t.obsp))
    assert_equal(adata.uns, t.uns)


def test_transpose_roundtrip(adata):
    del adata.raw
    assert_equal(adata, adata.T.T)
