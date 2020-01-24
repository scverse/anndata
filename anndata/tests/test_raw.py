import numpy as np
import pytest

import anndata as ad
from anndata._core.anndata import ImplicitModificationWarning
from anndata.tests.helpers import assert_equal


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
def adata_raw():
    adata = ad.AnnData(
        np.array(data), obs=obs_dict, var=var_dict, uns=uns_dict, dtype="int32"
    )
    adata.raw = adata
    # Make them different shapes
    adata = adata[:, [0, 1]].copy()
    return adata


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_raw_init(adata_raw):
    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_del(adata_raw):
    del adata_raw.raw
    assert adata_raw.raw is None


def test_raw_of_view(adata_raw):
    adata_view = adata_raw[adata_raw.obs["oanno1"] == "cat2"]
    assert adata_view.raw.X.tolist() == [
        [4, 5, 6],
        [7, 8, 9],
    ]


def test_raw_rw(adata_raw, backing_h5ad):
    adata_raw.write(backing_h5ad)
    adata_read = ad.read(backing_h5ad)

    assert_equal(adata_read, adata_raw, exact=True)

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_view_rw(adata_raw, backing_h5ad):
    # Make sure it still writes correctly if the object is a view
    adata_raw_view = adata_raw[:, adata_raw.var_names]
    assert_equal(adata_raw_view, adata_raw)
    with pytest.warns(ImplicitModificationWarning, match="Initializing view as actual"):
        adata_raw_view.write(backing_h5ad)
    adata_read = ad.read(backing_h5ad)

    assert_equal(adata_read, adata_raw_view, exact=True)

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    assert adata_raw.raw[:, 0].X.tolist() == [[1], [4], [7]]


def test_raw_backed(adata_raw, backing_h5ad):
    adata_raw.filename = backing_h5ad

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    if adata_raw.raw[:, 0].X.shape[1] != 1:
        pytest.xfail("Raw is broken for backed slices")
    assert adata_raw.raw[:, 0].X[:].tolist() == [[1], [4], [7]]


def test_raw_view_backed(adata_raw, backing_h5ad):
    adata_raw.filename = backing_h5ad

    assert adata_raw.var_names.tolist() == ["var1", "var2"]
    assert adata_raw.raw.var_names.tolist() == ["var1", "var2", "var3"]
    if adata_raw.raw[:, 0].X.shape[1] != 1:
        pytest.xfail("Raw is broken for backed slices")
    assert adata_raw.raw[:, 0].X[:].tolist() == [[1], [4], [7]]


def test_raw_as_parent_view():
    # https://github.com/theislab/anndata/issues/288
    a = ad.AnnData(np.ones((4, 3)))
    a.varm["PCs"] = np.ones((3, 3))
    a.raw = a
    # create a Raw containing views. This used to trigger #288.
    b = a.raw[:, "0"]
    # actualize
    b.varm["PCs"] = np.array([[1, 2, 3]])
