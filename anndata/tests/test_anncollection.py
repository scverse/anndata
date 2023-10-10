from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix, issparse
from sklearn.preprocessing import LabelEncoder

import anndata as ad
from anndata.experimental.multi_files import AnnCollection

_dense = lambda a: a.toarray() if issparse(a) else a


@pytest.fixture
def adatas(request):
    adata1 = ad.AnnData(X=request.param([[1, 2, 0], [4, 5, 0], [7, 8, 0]]))
    adata1.obs["a_test"] = ["a", "a", "b"]
    adata1.obsm["o_test"] = np.ones((adata1.n_obs, 2))

    adata2 = ad.AnnData(X=request.param([[1, 3, 0], [9, 8, 0]]))
    adata2.obs["a_test"] = ["c", "c"]
    adata2.obsm["o_test"] = np.zeros((adata2.n_obs, 2))

    return adata1, adata2


@pytest.mark.parametrize("adatas", [np.array, csr_matrix], indirect=True)
def test_full_selection(adatas):
    dat = AnnCollection(adatas, index_unique="_")
    adt_concat = ad.concat(adatas, index_unique="_")

    # sorted selection from one adata
    dat_view = dat[:2, :2]
    for adata in (adatas[0], adt_concat):
        adt_view = adata[:2, :2]
        np.testing.assert_allclose(_dense(dat_view.X), _dense(adt_view.X))
        np.testing.assert_allclose(dat_view.obsm["o_test"], adt_view.obsm["o_test"])
        np.testing.assert_array_equal(dat_view.obs["a_test"], adt_view.obs["a_test"])

    # sorted and unsorted selection from 2 adatas
    rand_idxs = np.random.choice(dat.shape[0], 4, replace=False)
    for select in (slice(2, 5), [4, 2, 3], rand_idxs):
        dat_view = dat[select, :2]
        adt_view = adt_concat[select, :2]
        np.testing.assert_allclose(_dense(dat_view.X), _dense(adt_view.X))
        np.testing.assert_allclose(dat_view.obsm["o_test"], adt_view.obsm["o_test"])
        np.testing.assert_array_equal(dat_view.obs["a_test"], adt_view.obs["a_test"])

    # test duplicate selection
    idxs = [1, 2, 4, 4]
    dat_view = dat[idxs, :2]
    np.testing.assert_allclose(
        _dense(dat_view.X), np.array([[4, 5], [7, 8], [9, 8], [9, 8]])
    )


@pytest.mark.parametrize("adatas", [np.array, csr_matrix], indirect=True)
def test_creation(adatas):
    adatas_inner = [adatas[0], adatas[1][:, :2].copy()]

    dat = AnnCollection(adatas_inner, join_vars="inner", index_unique="_")
    adt_concat = ad.concat(adatas_inner, index_unique="_")
    np.testing.assert_array_equal(dat.var_names, adt_concat.var_names)


@pytest.mark.parametrize("adatas", [np.array], indirect=True)
def test_convert(adatas):
    dat = AnnCollection(adatas, index_unique="_")

    le = LabelEncoder()
    le.fit(dat[:].obs["a_test"])

    obs_no_convert = dat[:].obs["a_test"]
    convert = dict(obs={"a_test": lambda a: le.transform(a)})
    dat.convert = convert
    np.testing.assert_array_equal(dat[:].obs["a_test"], le.transform(obs_no_convert))
