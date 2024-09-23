from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
from anndata import ImplicitModificationWarning
from anndata.tests.helpers import GEN_ADATA_DASK_ARGS, assert_equal, gen_adata

# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

data = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
]  # data matrix of shape n_obs × n_vars

obs_dict = dict(  # annotation of observations / rows
    row_names=["name1", "name2", "name3"],  # row annotation
    oanno1=["cat1", "cat2", "cat2"],  # categorical annotation
    oanno2=["o1", "o2", "o3"],  # string annotation
    oanno3=[2.1, 2.2, 2.3],  # float annotation
)

var_dict = dict(  # annotation of variables / columns
    col_names=["var1", "var2", "var3"], vanno1=[3.1, 3.2, 3.3]
)

uns_dict = dict(  # unstructured annotation
    oanno1_colors=["#000000", "#FFFFFF"], uns2=["some annotation"]
)


@pytest.fixture
def adata_raw() -> ad.AnnData:
    adata = ad.AnnData(
        np.array(data, dtype="int32"), obs=obs_dict, var=var_dict, uns=uns_dict
    )
    adata.raw = adata.copy()
    # Make them different shapes
    adata = adata[:, [0, 1]].copy()
    return adata


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_raw_init(adata_raw: ad.AnnData):
    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_del(adata_raw: ad.AnnData):
    del adata_raw.raw
    assert adata_raw.raw is None


def test_raw_set_as_none(adata_raw: ad.AnnData):
    # Test for scverse/anndata#445
    a = adata_raw
    b = adata_raw.copy()

    del a.raw
    b.raw = None

    assert_equal(a, b)


def test_raw_of_view(adata_raw: ad.AnnData):
    adata_view = adata_raw[adata_raw.obs["oanno1"] == "cat2"]
    assert adata_view.raw.X.tolist() == [
        [4, 5, 6],
        [7, 8, 9],
    ]


def test_raw_rw(adata_raw: ad.AnnData, backing_h5ad):
    adata_raw.write(backing_h5ad)
    adata_read = ad.read_h5ad(backing_h5ad)

    assert_equal(adata_read, adata_raw, exact=True)

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_view_rw(adata_raw: ad.AnnData, backing_h5ad):
    # Make sure it still writes correctly if the object is a view
    adata_raw_view = adata_raw[:, adata_raw.var_names]
    assert_equal(adata_raw_view, adata_raw)
    with pytest.warns(
        ImplicitModificationWarning, match=r"initializing view as actual"
    ):
        adata_raw_view.write(backing_h5ad)
    adata_read = ad.read_h5ad(backing_h5ad)

    assert_equal(adata_read, adata_raw_view, exact=True)

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_backed(adata_raw: ad.AnnData, backing_h5ad):
    adata_raw.filename = backing_h5ad

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    if adata_raw.raw[:, 0].X.shape[1] != 1:
        pytest.xfail("Raw is broken for backed slices")
    assert adata_raw.raw[:, 0].X[:].tolist() == [[1], [4], [7]]


def test_raw_view_backed(adata_raw: ad.AnnData, backing_h5ad):
    adata_raw.filename = backing_h5ad

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    if adata_raw.raw[:, 0].X.shape[1] != 1:
        pytest.xfail("Raw is broken for backed slices")
    assert adata_raw.raw[:, 0].X[:].tolist() == [[1], [4], [7]]


def test_raw_as_parent_view():
    # https://github.com/scverse/anndata/issues/288
    a = ad.AnnData(np.ones((4, 3)))
    a.varm["PCs"] = np.ones((3, 3))
    a.raw = a.copy()
    # create a Raw containing views. This used to trigger #288.
    b = a.raw[:, "0"]
    # actualize
    b.varm["PCs"] = np.array([[1, 2, 3]])


def test_to_adata():
    # https://github.com/scverse/anndata/pull/404
    adata = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)

    with_raw = adata[:, ::2].copy()
    with_raw.raw = adata.copy()

    # Raw doesn't do layers or varp currently
    # Deleting after creation so we know to rewrite the test if they are supported
    del adata.layers, adata.varp

    assert_equal(adata, with_raw.raw.to_adata())


def test_to_adata_populates_obs():
    adata = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)

    del adata.layers, adata.uns, adata.varp
    adata_w_raw = adata.copy()

    raw = adata.copy()
    del raw.obs, raw.obsm, raw.obsp, raw.uns

    adata_w_raw.raw = raw
    from_raw = adata_w_raw.raw.to_adata()

    assert_equal(adata, from_raw)


def test_no_copy():
    adata = gen_adata((20, 10), X_type=np.asarray)
    adata.raw = adata  # no .copy() herer
    np.log1p(adata.X, out=adata.X)
    assert adata.X is adata.raw.X
