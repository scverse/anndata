import pytest
import anndata as ad
import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix, issparse
from anndata.experimental.pytorch import AnnLoader

_dense = lambda a: a.toarray() if issparse(a) else a
torch = pytest.importorskip("torch")


@pytest.fixture
def adatas(request):
    adata1 = ad.AnnData(X=request.param([[1, 2, 0, 0], [4, 5, 0, 0], [7, 8, 0, 0]]))

    adata1.obs.index = ["1", "2", "3"]
    adata1.obs["o_test"] = ["a", "b", "c"]

    adata1.var.index = ["1", "2", "3", "4"]
    adata1.var["v_test"] = ["x", "y", "z", "w"]

    adata2 = ad.AnnData(X=request.param([[1, 3, 1, 0], [9, 8, 0, 1]]))

    adata2.obs.index = ["4", "5"]
    adata2.obs["o_test"] = ["d", "e"]

    adata2.var.index = ["1", "2", "3", "4"]
    adata2.var["v_test"] = ["x", "y", "z", "w"]

    return adata1, adata2


@pytest.fixture
def adata(request):
    adata1 = ad.AnnData(
        X=request.param(
            [
                [1, 2, 0, 0, 1],
                [4, 5, 0, 0, 0],
                [7, 8, 0, 0, 0],
                [1, 3, 1, 0, 0],
                [9, 8, 0, 1, 0],
            ]
        )
    )

    adata1.obs.index = ["1", "2", "3", "4", "5"]
    adata1.obs["o_test"] = ["a", "b", "c", "d", "e"]

    adata1.var.index = ["1", "2", "3", "4", "5"]
    adata1.var["v_test"] = ["x", "y", "z", "w", "v"]

    return adata1


@pytest.mark.parametrize("adata", [np.array, csr_matrix], indirect=True)
def test_annloader_obs(adata):
    loader = AnnLoader(adata, batch_size=2, axis="obs")

    for i, x in enumerate(loader):
        np.testing.assert_allclose(x.X, _dense(adata.X[2 * i : 2 * (i + 1), :]))
        pd.testing.assert_frame_equal(x.obs.df, adata.obs[2 * i : 2 * (i + 1)])
        pd.testing.assert_frame_equal(x.var.df, adata.var)


@pytest.mark.parametrize("adata", [np.array, csr_matrix], indirect=True)
def test_annloader_vars(adata):
    loader = AnnLoader(adata, batch_size=2, axis="vars")

    for i, x in enumerate(loader):
        np.testing.assert_allclose(x.X, _dense(adata.X[:, 2 * i : 2 * (i + 1)]))
        pd.testing.assert_frame_equal(x.obs.df, adata.obs)
        pd.testing.assert_frame_equal(x.var.df, adata.var[2 * i : 2 * (i + 1)])
