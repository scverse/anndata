from __future__ import annotations

import re
import warnings
from functools import partial
from itertools import product
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from numpy import ma
from scipy import sparse as sp
from scipy.sparse import csr_matrix, issparse

import anndata as ad
from anndata import AnnData, ImplicitModificationWarning
from anndata._core.raw import Raw
from anndata._settings import settings
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    assert_equal,
    gen_adata,
    get_multiindex_columns_df,
)

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

# some test objects that we use below
adata_dense = AnnData(np.array([[1, 2], [3, 4]]))
adata_dense.layers["test"] = adata_dense.X
adata_sparse = AnnData(
    csr_matrix([[0, 2, 3], [0, 5, 6]]),
    dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
    dict(var_names=["a", "b", "c"]),
)


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(np.array([[1, 2], [3, 4]]), {}, {})
    AnnData(ma.array([[1, 2], [3, 4]]), uns=dict(mask=[0, 1, 1, 0]))
    AnnData(sp.eye(2, format="csr"))
    AnnData(sp.csr_array([[1, 0], [0, 1]]))
    X = np.array([[1, 2, 3], [4, 5, 6]])
    adata = AnnData(
        X=X,
        obs=dict(Obs=["A", "B"]),
        var=dict(Feat=["a", "b", "c"]),
        obsm=dict(X_pca=np.array([[1, 2], [3, 4]])),
        raw=dict(X=X, var=dict(var_names=["a", "b", "c"])),
    )

    assert adata.raw.X.tolist() == X.tolist()
    assert adata.raw.var_names.tolist() == ["a", "b", "c"]

    # init with empty data matrix
    shape = (3, 5)
    adata = AnnData(None, uns=dict(test=np.array((3, 3))), shape=shape)
    assert adata.X is None
    assert adata.shape == shape
    assert "test" in adata.uns


@pytest.mark.parametrize(
    ("src", "src_arg", "dim_msg"),
    [
        pytest.param(
            "X",
            adata_dense.X,
            "`{dim}` must have as many rows as `X` has {mat_dim}s",
            id="x",
        ),
        pytest.param(
            "shape", (2, 2), "`shape` is inconsistent with `{dim}`", id="shape"
        ),
    ],
)
@pytest.mark.parametrize("dim", ["obs", "var"])
@pytest.mark.parametrize(
    ("dim_arg", "msg"),
    [
        pytest.param(
            lambda _: dict(TooLong=[1, 2, 3, 4]),
            "Length of values (4) does not match length of index (2)",
            id="too_long_col",
        ),
        pytest.param(
            lambda dim: {f"{dim}_names": ["a", "b", "c"]}, None, id="too_many_names"
        ),
        pytest.param(
            lambda _: pd.DataFrame(index=["a", "b", "c"]), None, id="too_long_df"
        ),
    ],
)
def test_creation_error(src, src_arg, dim_msg, dim, dim_arg, msg: str | None):
    if msg is None:
        mat_dim = "row" if dim == "obs" else "column"
        msg = dim_msg.format(dim=dim, mat_dim=mat_dim)
    with pytest.raises(ValueError, match=re.escape(msg)):
        AnnData(**{src: src_arg, dim: dim_arg(dim)})


def test_invalid_X():
    with pytest.raises(
        ValueError,
        match=r"X needs to be of one of <class 'numpy.ndarray'>.*not <class 'str'>\.",
    ):
        AnnData("string is not a valid X")


def test_create_with_dfs():
    X = np.ones((6, 3))
    obs = pd.DataFrame(dict(cat_anno=pd.Categorical(["a", "a", "a", "a", "b", "a"])))
    obs_copy = obs.copy()
    adata = AnnData(X=X, obs=obs)
    assert obs.index.equals(obs_copy.index)
    assert obs.index.astype(str).equals(adata.obs.index)


def test_create_from_df():
    df = pd.DataFrame(np.ones((3, 2)), index=["a", "b", "c"], columns=["A", "B"])
    ad = AnnData(df)
    assert df.values.tolist() == ad.X.tolist()
    assert df.columns.tolist() == ad.var_names.tolist()
    assert df.index.tolist() == ad.obs_names.tolist()


@pytest.mark.parametrize("attr", ["X", "obs", "obsm"])
def test_error_create_from_multiindex_df(attr):
    df = get_multiindex_columns_df((100, 20))
    val = df if attr != "obsm" else {"df": df}
    with pytest.raises(ValueError, match=r"MultiIndex columns are not supported"):
        AnnData(**{attr: val}, shape=(100, 10))


def test_create_from_sparse_df():
    s = sp.random(20, 30, density=0.2, format="csr")
    obs_names = [f"obs{i}" for i in range(20)]
    var_names = [f"var{i}" for i in range(30)]
    df = pd.DataFrame.sparse.from_spmatrix(s, index=obs_names, columns=var_names)
    a = AnnData(df)
    b = AnnData(s, obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=var_names))
    assert_equal(a, b)
    assert issparse(a.X)


def test_create_from_df_with_obs_and_var():
    df = pd.DataFrame(np.ones((3, 2)), index=["a", "b", "c"], columns=["A", "B"])
    obs = pd.DataFrame(np.ones((3, 1)), index=df.index, columns=["C"])
    var = pd.DataFrame(np.ones((2, 1)), index=df.columns, columns=["D"])
    ad = AnnData(df, obs=obs, var=var)
    assert df.values.tolist() == ad.X.tolist()
    assert df.columns.tolist() == ad.var_names.tolist()
    assert df.index.tolist() == ad.obs_names.tolist()
    assert obs.equals(ad.obs)
    assert var.equals(ad.var)

    with pytest.raises(ValueError, match=r"Index of obs must match index of X."):
        AnnData(df, obs=obs.reset_index())
    with pytest.raises(ValueError, match=r"Index of var must match columns of X."):
        AnnData(df, var=var.reset_index())


def test_matching_int_index():
    adata = AnnData(
        pd.DataFrame(dict(a=[0.0, 0.5]), index=[0, 1]), obs=pd.DataFrame(index=[0, 1])
    )
    pd.testing.assert_index_equal(adata.obs_names, pd.Index(["0", "1"]))


def test_from_df_and_dict():
    df = pd.DataFrame(dict(a=[0.1, 0.2, 0.3], b=[1.1, 1.2, 1.3]))
    adata = AnnData(df, dict(species=pd.Categorical(["a", "b", "a"])))
    assert adata.obs["species"].values.tolist() == ["a", "b", "a"]


def test_df_warnings():
    df = pd.DataFrame(dict(A=[1, 2, 3], B=[1.0, 2.0, 3.0]), index=["a", "b", "c"])
    with pytest.warns(UserWarning, match=r"X.*dtype float64"):
        adata = AnnData(df)
    with pytest.warns(UserWarning, match=r"X.*dtype float64"):
        adata.X = df


@pytest.mark.parametrize("attr", ["X", "layers", "obsm", "varm", "obsp", "varp"])
@pytest.mark.parametrize("when", ["init", "assign"])
def test_convert_matrix(attr, when):
    """Test that initializing or assigning aligned arrays to a np.matrix converts it."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"the matrix.*not.*recommended", PendingDeprecationWarning
        )
        mat = np.matrix([[1, 2], [3, 0]])

    direct = attr in {"X"}

    with pytest.warns(ImplicitModificationWarning, match=r"np\.ndarray"):
        if when == "init":
            adata = (
                AnnData(**{attr: mat})
                if direct
                else AnnData(shape=(2, 2), **{attr: {"a": mat}})
            )
        elif when == "assign":
            adata = AnnData(shape=(2, 2))
            if direct:
                setattr(adata, attr, mat)
            else:
                getattr(adata, attr)["a"] = mat
        else:
            raise ValueError(when)

    arr = getattr(adata, attr) if direct else getattr(adata, attr)["a"]
    assert isinstance(arr, np.ndarray), f"{arr} is not an array"
    assert not isinstance(arr, np.matrix), f"{arr} is still a matrix"


def test_attr_deletion():
    full = gen_adata((30, 30))
    # Empty has just X, obs_names, var_names
    empty = AnnData(None, obs=full.obs[[]], var=full.var[[]])
    for attr in ["X", "obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns"]:
        delattr(full, attr)
        assert_equal(getattr(full, attr), getattr(empty, attr))
    assert_equal(full, empty, exact=True)


def test_names():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=["A", "B"]),
        dict(var_names=["a", "b", "c"]),
    )

    assert adata.obs_names.tolist() == ["A", "B"]
    assert adata.var_names.tolist() == ["a", "b", "c"]

    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]), var=dict(var_names=["a", "b"]))
    assert adata.var_names.tolist() == ["a", "b"]


@pytest.mark.parametrize(
    ("names", "after"),
    [
        pytest.param(["a", "b"], None, id="list"),
        pytest.param(
            pd.Series(["AAD", "CCA"], name="barcodes"), "barcodes", id="Series-str"
        ),
        pytest.param(pd.Series(["x", "y"], name=0), None, id="Series-int"),
    ],
)
@pytest.mark.parametrize("attr", ["obs_names", "var_names"])
def test_setting_index_names(names, after, attr):
    adata = adata_dense.copy()
    assert getattr(adata, attr).name is None
    setattr(adata, attr, names)
    assert getattr(adata, attr).name == after
    if hasattr(names, "name"):
        assert names.name is not None

    # Testing for views
    new = adata[:, :]
    assert new.is_view
    setattr(new, attr, names)
    assert_equal(new, adata, exact=True)
    assert not new.is_view


@pytest.mark.parametrize("attr", ["obs_names", "var_names"])
def test_setting_index_names_error(attr):
    orig = adata_sparse[:2, :2]
    adata = adata_sparse[:2, :2]
    assert getattr(adata, attr).name is None
    with pytest.raises(ValueError, match=rf"AnnData expects \.{attr[:3]}\.index\.name"):
        setattr(adata, attr, pd.Index(["x", "y"], name=0))
    assert adata.is_view
    assert getattr(adata, attr).tolist() != ["x", "y"]
    assert getattr(adata, attr).tolist() == getattr(orig, attr).tolist()
    assert_equal(orig, adata, exact=True)


@pytest.mark.parametrize("dim", ["obs", "var"])
@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset"), [(False, False), (True, True)]
)
def test_setting_dim_index(dim, obs_xdataset, var_xdataset):
    index_attr = f"{dim}_names"
    mapping_attr = f"{dim}m"

    orig = gen_adata((5, 5), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset)
    orig.raw = orig.copy()
    curr = orig.copy()
    view = orig[:, :]
    new_idx = pd.Index(list("abcde"), name="letters")

    setattr(curr, index_attr, new_idx)
    pd.testing.assert_index_equal(getattr(curr, index_attr), new_idx)
    pd.testing.assert_index_equal(getattr(curr, mapping_attr)["df"].index, new_idx)
    pd.testing.assert_index_equal(getattr(curr, mapping_attr).dim_names, new_idx)
    pd.testing.assert_index_equal(curr.obs_names, curr.raw.obs_names)

    # Testing view behaviour
    setattr(view, index_attr, new_idx)
    assert not view.is_view
    pd.testing.assert_index_equal(getattr(view, index_attr), new_idx)
    pd.testing.assert_index_equal(getattr(view, mapping_attr)["df"].index, new_idx)
    pd.testing.assert_index_equal(getattr(view, mapping_attr).dim_names, new_idx)
    with pytest.raises(AssertionError):
        pd.testing.assert_index_equal(
            getattr(view, index_attr), getattr(orig, index_attr)
        )
    assert_equal(view, curr, exact=True)

    # test case in #459
    fake_m = pd.DataFrame(curr.X.T, index=getattr(curr, index_attr))
    getattr(curr, mapping_attr)["df2"] = fake_m


def test_indices_dtypes():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=["A", "B"]),
        dict(var_names=["a", "b", "c"]),
    )
    adata.obs_names = ["ö", "a"]
    assert adata.obs_names.tolist() == ["ö", "a"]


def test_slicing():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    # assert adata[:, 0].X.tolist() == adata.X[:, 0].tolist()  # No longer the case

    assert adata[0, 0].X.tolist() == np.reshape(1, (1, 1)).tolist()
    assert adata[0, :].X.tolist() == np.reshape([1, 2, 3], (1, 3)).tolist()
    assert adata[:, 0].X.tolist() == np.reshape([1, 4], (2, 1)).tolist()

    assert adata[:, [0, 1]].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array([0, 2])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, np.array([False, True, True])].X.tolist() == [
        [2, 3],
        [5, 6],
    ]
    assert adata[:, 1:3].X.tolist() == [[2, 3], [5, 6]]

    assert adata[0:2, :][:, 0:2].X.tolist() == [[1, 2], [4, 5]]
    assert adata[0:1, :][:, 0:2].X.tolist() == np.reshape([1, 2], (1, 2)).tolist()
    assert adata[0, :][:, 0].X.tolist() == np.reshape(1, (1, 1)).tolist()
    assert adata[:, 0:2][0:2, :].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, 0:2][0:1, :].X.tolist() == np.reshape([1, 2], (1, 2)).tolist()
    assert adata[:, 0][0, :].X.tolist() == np.reshape(1, (1, 1)).tolist()


def test_boolean_slicing():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    obs_selector = np.array([True, False], dtype=bool)
    vars_selector = np.array([True, False, False], dtype=bool)
    assert adata[obs_selector, :][:, vars_selector].X.tolist() == [[1]]
    assert adata[:, vars_selector][obs_selector, :].X.tolist() == [[1]]
    assert adata[obs_selector, :][:, 0].X.tolist() == [[1]]
    assert adata[:, 0][obs_selector, :].X.tolist() == [[1]]
    assert adata[0, :][:, vars_selector].X.tolist() == [[1]]
    assert adata[:, vars_selector][0, :].X.tolist() == [[1]]

    obs_selector = np.array([True, False], dtype=bool)
    vars_selector = np.array([True, True, False], dtype=bool)
    assert adata[obs_selector, :][:, vars_selector].X.tolist() == [[1, 2]]
    assert adata[:, vars_selector][obs_selector, :].X.tolist() == [[1, 2]]
    assert adata[obs_selector, :][:, 0:2].X.tolist() == [[1, 2]]
    assert adata[:, 0:2][obs_selector, :].X.tolist() == [[1, 2]]
    assert adata[0, :][:, vars_selector].X.tolist() == [[1, 2]]
    assert adata[:, vars_selector][0, :].X.tolist() == [[1, 2]]

    obs_selector = np.array([True, True], dtype=bool)
    vars_selector = np.array([True, True, False], dtype=bool)
    assert adata[obs_selector, :][:, vars_selector].X.tolist() == [
        [1, 2],
        [4, 5],
    ]
    assert adata[:, vars_selector][obs_selector, :].X.tolist() == [
        [1, 2],
        [4, 5],
    ]
    assert adata[obs_selector, :][:, 0:2].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, 0:2][obs_selector, :].X.tolist() == [[1, 2], [4, 5]]
    assert adata[0:2, :][:, vars_selector].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, vars_selector][0:2, :].X.tolist() == [[1, 2], [4, 5]]


def test_oob_boolean_slicing():
    len1, len2 = np.random.choice(100, 2, replace=False)
    with pytest.raises(IndexError) as e:
        AnnData(np.empty((len1, 100)))[np.random.randint(0, 2, len2, dtype=bool), :]
    assert str(len1) in str(e.value)
    assert str(len2) in str(e.value)

    len1, len2 = np.random.choice(100, 2, replace=False)
    with pytest.raises(IndexError) as e:
        AnnData(np.empty((100, len1)))[:, np.random.randint(0, 2, len2, dtype=bool)]
    assert str(len1) in str(e.value)
    assert str(len2) in str(e.value)


def test_slicing_strings():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=["A", "B"]),
        dict(var_names=["a", "b", "c"]),
    )

    assert adata["A", "a"].X.tolist() == [[1]]
    assert adata["A", :].X.tolist() == [[1, 2, 3]]
    assert adata[:, "a"].X.tolist() == [[1], [4]]
    assert adata[:, ["a", "b"]].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array(["a", "c"])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, "b":"c"].X.tolist() == [[2, 3], [5, 6]]

    with pytest.raises(KeyError):
        _ = adata[:, "X"]
    with pytest.raises(KeyError):
        _ = adata["X", :]
    with pytest.raises(KeyError):
        _ = adata["A":"X", :]
    with pytest.raises(KeyError):
        _ = adata[:, "a":"X"]

    # Test if errors are helpful
    with pytest.raises(KeyError, match=r"not_in_var"):
        adata[:, ["A", "B", "not_in_var"]]
    with pytest.raises(KeyError, match=r"not_in_obs"):
        adata[["A", "B", "not_in_obs"], :]


def test_slicing_series():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6]]),
        dict(obs_names=["A", "B", "C"]),
        dict(var_names=["a", "b"]),
    )
    df = pd.DataFrame(dict(a=["1", "2", "2"]))
    df1 = pd.DataFrame(dict(b=["1", "2"]))
    assert adata[df["a"].values == "2"].X.tolist() == adata[df["a"] == "2"].X.tolist()
    assert (
        adata[:, df1["b"].values == "2"].X.tolist()
        == adata[:, df1["b"] == "2"].X.tolist()
    )


def test_strings_to_categoricals():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6], [7, 8]]), dict(k=["a", "a", "b", "b"])
    )
    adata.strings_to_categoricals()
    assert adata.obs["k"].cat.categories.tolist() == ["a", "b"]


def test_slicing_remove_unused_categories():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6], [7, 8]]), dict(k=["a", "a", "b", "b"])
    )
    adata._sanitize()
    assert adata[2:4].obs["k"].cat.categories.tolist() == ["b"]


def test_slicing_dont_remove_unused_categories():
    with settings.override(remove_unused_categories=False):
        adata = AnnData(
            np.array([[1, 2], [3, 4], [5, 6], [7, 8]]), dict(k=["a", "a", "b", "b"])
        )
        adata._sanitize()
        assert adata[2:4].obs["k"].cat.categories.tolist() == ["a", "b"]


def test_no_uniqueness_check_gives_repeat_indices():
    with settings.override(check_uniqueness=False):
        obs_names = ["0", "0", "1", "1"]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            adata = AnnData(
                np.array([[1, 2], [3, 4], [5, 6], [7, 8]]),
                obs=pd.DataFrame(index=obs_names),
            )
    assert adata.obs_names.values.tolist() == obs_names


def test_get_subset_annotation():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(S=["A", "B"]),
        dict(F=["a", "b", "c"]),
    )

    assert adata[0, 0].obs["S"].tolist() == ["A"]
    assert adata[0, 0].var["F"].tolist() == ["a"]


def test_append_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.obs["new"] = [1, 2]
    # this worked in the initial AnnData, but not with a dataframe
    # adata.obs[['new2', 'new3']] = [['A', 'B'], ['c', 'd']]

    with pytest.raises(
        ValueError, match="Length of values.*does not match length of index"
    ):
        adata.obs["new4"] = ["far", "too", "long"]


def test_delete_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]), dict(o1=[1, 2], o2=[3, 4]))
    assert adata.obs_keys() == ["o1", "o2"]

    del adata.obs["o1"]
    assert adata.obs_keys() == ["o2"]
    assert adata.obs["o2"].tolist() == [3, 4]


def test_set_obs():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.obs = pd.DataFrame(dict(a=[3, 4]))
    assert adata.obs_names.tolist() == ["0", "1"]

    with pytest.raises(ValueError, match="`shape` is inconsistent with `obs`"):
        adata.obs = pd.DataFrame(dict(a=[3, 4, 5]))


def test_multicol():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    # 'c' keeps the columns as should be
    adata.obsm["c"] = np.array([[0.0, 1.0], [2, 3]])
    assert adata.obsm_keys() == ["c"]
    assert adata.obsm["c"].tolist() == [[0.0, 1.0], [2, 3]]


def test_n_obs():
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    assert adata.n_obs == 3
    adata1 = adata[:2]
    assert adata1.n_obs == 2


def test_equality_comparisons():
    adata1 = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    adata2 = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    with pytest.raises(NotImplementedError):
        adata1 == adata1  # noqa: B015, PLR0124
    with pytest.raises(NotImplementedError):
        adata1 == adata2  # noqa: B015
    with pytest.raises(NotImplementedError):
        adata1 != adata2  # noqa: B015
    with pytest.raises(NotImplementedError):
        adata1 == 1  # noqa: B015
    with pytest.raises(NotImplementedError):
        adata1 != 1  # noqa: B015


def test_rename_categories():
    X = np.ones((6, 3))
    obs = pd.DataFrame(dict(cat_anno=pd.Categorical(["a", "a", "a", "a", "b", "a"])))
    adata = AnnData(X=X, obs=obs)
    adata.uns["tool"] = {}
    adata.uns["tool"]["cat_array"] = np.rec.fromarrays(
        [np.ones(2) for cat in adata.obs["cat_anno"].cat.categories],
        dtype=[(cat, "float32") for cat in adata.obs["cat_anno"].cat.categories],
    )
    adata.uns["tool"]["params"] = dict(groupby="cat_anno")

    new_categories = ["c", "d"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        adata.rename_categories("cat_anno", new_categories)

    assert list(adata.obs["cat_anno"].cat.categories) == new_categories
    assert list(adata.uns["tool"]["cat_array"].dtype.names) == new_categories


def test_pickle():
    import pickle

    adata = AnnData()
    adata2 = pickle.loads(pickle.dumps(adata))
    assert adata2.obsm.parent is adata2


def test_to_df_dense():
    X_df = adata_dense.to_df()
    layer_df = adata_dense.to_df(layer="test")

    np.testing.assert_array_equal(adata_dense.layers["test"], layer_df.values)
    np.testing.assert_array_equal(adata_dense.X, X_df.values)
    pd.testing.assert_index_equal(X_df.columns, layer_df.columns)
    pd.testing.assert_index_equal(X_df.index, layer_df.index)


def test_convenience():
    adata = adata_sparse.copy()
    adata.layers["x2"] = adata.X * 2
    adata.var["anno2"] = ["p1", "p2", "p3"]
    adata.raw = adata.copy()
    adata.X = adata.X / 2
    adata_dense = adata.copy()
    adata_dense.X = adata_dense.X.toarray()

    def assert_same_op_result(a1, a2, op):
        r1 = op(a1)
        r2 = op(a2)
        assert np.all(r1 == r2)
        assert type(r1) is type(r2)

    assert np.allclose(adata.obs_vector("b"), np.array([1.0, 2.5]))
    assert np.allclose(adata.raw.obs_vector("c"), np.array([3, 6]))
    assert np.all(adata.obs_vector("anno1") == np.array(["c1", "c2"]))
    assert np.allclose(adata.var_vector("s1"), np.array([0, 1.0, 1.5]))
    assert np.allclose(adata.raw.var_vector("s2"), np.array([0, 5, 6]))

    for obs_k, layer in product(["a", "b", "c", "anno1"], [None, "x2"]):
        assert_same_op_result(
            adata, adata_dense, partial(AnnData.obs_vector, k=obs_k, layer=layer)
        )

    for obs_k in ["a", "b", "c"]:
        assert_same_op_result(
            adata.raw, adata_dense.raw, partial(Raw.obs_vector, k=obs_k)
        )

    for var_k, layer in product(["s1", "s2", "anno2"], [None, "x2"]):
        assert_same_op_result(
            adata, adata_dense, partial(AnnData.var_vector, k=var_k, layer=layer)
        )

    for var_k in ["s1", "s2", "anno2"]:
        assert_same_op_result(
            adata.raw, adata_dense.raw, partial(Raw.var_vector, k=var_k)
        )


def test_1d_slice_dtypes():
    N, M = 10, 20
    obs_df = pd.DataFrame(
        dict(
            cat=pd.Categorical(np.arange(N, dtype=int)),
            int=np.arange(N, dtype=int),
            float=np.arange(N, dtype=float),
            obj=[str(i) for i in np.arange(N, dtype=int)],
        ),
        index=[f"cell{i}" for i in np.arange(N, dtype=int)],
    )
    var_df = pd.DataFrame(
        dict(
            cat=pd.Categorical(np.arange(M, dtype=int)),
            int=np.arange(M, dtype=int),
            float=np.arange(M, dtype=float),
            obj=[str(i) for i in np.arange(M, dtype=int)],
        ),
        index=[f"gene{i}" for i in np.arange(M, dtype=int)],
    )
    adata = AnnData(X=np.random.random((N, M)), obs=obs_df, var=var_df)

    new_obs_df = pd.DataFrame(index=adata.obs_names)
    for k in obs_df.columns:
        new_obs_df[k] = adata.obs_vector(k)
        assert new_obs_df[k].dtype == obs_df[k].dtype
    assert np.all(new_obs_df == obs_df)
    new_var_df = pd.DataFrame(index=adata.var_names)
    for k in var_df.columns:
        new_var_df[k] = adata.var_vector(k)
        assert new_var_df[k].dtype == var_df[k].dtype
    assert np.all(new_var_df == var_df)


def test_to_df_sparse():
    X = adata_sparse.X.toarray()
    df = adata_sparse.to_df()
    assert df.values.tolist() == X.tolist()


def test_to_df_no_X():
    adata = AnnData(
        obs=pd.DataFrame(index=[f"cell-{i:02}" for i in range(20)]),
        var=pd.DataFrame(index=[f"gene-{i:02}" for i in range(30)]),
        layers={"present": np.ones((20, 30))},
    )
    v = adata[:10]

    with pytest.raises(ValueError, match=r"X is None"):
        _ = adata.to_df()
    with pytest.raises(ValueError, match=r"X is None"):
        _ = v.to_df()

    expected = pd.DataFrame(
        np.ones(adata.shape), index=adata.obs_names, columns=adata.var_names
    )
    actual = adata.to_df(layer="present")

    pd.testing.assert_frame_equal(actual, expected)

    view_expected = pd.DataFrame(
        np.ones(v.shape), index=v.obs_names, columns=v.var_names
    )
    view_actual = v.to_df(layer="present")

    pd.testing.assert_frame_equal(view_actual, view_expected)


def test_copy():
    adata_copy = adata_sparse.copy()

    def assert_eq_not_id(a, b):
        assert a is not b
        assert issparse(a) == issparse(b)
        if issparse(a):
            assert np.all(a.data == b.data)
            assert np.all(a.indices == b.indices)
            assert np.all(a.indptr == b.indptr)
        else:
            assert np.all(a == b)

    assert adata_sparse is not adata_copy
    assert_eq_not_id(adata_sparse.X, adata_copy.X)
    for attr in ["layers", "var", "obs", "obsm", "varm"]:
        map_sprs = getattr(adata_sparse, attr)
        map_copy = getattr(adata_copy, attr)
        assert map_sprs is not map_copy
        if attr not in {"obs", "var"}:
            # check that we don’t create too many references
            assert getattr(adata_copy, f"_{attr}") is map_copy._data
        assert_eq_not_id(map_sprs.keys(), map_copy.keys())
        for key in map_sprs:
            assert_eq_not_id(map_sprs[key], map_copy[key])


def test_to_memory_no_copy():
    adata = gen_adata((3, 5), **GEN_ADATA_NO_XARRAY_ARGS)
    mem = adata.to_memory()

    assert mem.X is adata.X
    # Currently does not hold for `obs`, `var`, but should in future
    for key in adata.layers:
        assert mem.layers[key] is adata.layers[key]
    for key in adata.obsm:
        assert mem.obsm[key] is adata.obsm[key]
    for key in adata.varm:
        assert mem.varm[key] is adata.varm[key]
    for key in adata.obsp:
        assert mem.obsp[key] is adata.obsp[key]
    for key in adata.varp:
        assert mem.varp[key] is adata.varp[key]


@pytest.mark.parametrize("axis", ["obs", "var"])
@pytest.mark.parametrize("elem_type", ["p", "m"])
def test_create_adata_from_single_axis_elem(
    axis: Literal["obs", "var"], elem_type: Literal["m", "p"], tmp_path: Path
):
    d = dict(
        a=np.zeros((10, 10)),
    )
    in_memory = AnnData(**{f"{axis}{elem_type}": d})
    assert in_memory.shape == (10, 0) if axis == "obs" else (0, 10)
    in_memory.write_h5ad(tmp_path / "adata.h5ad")
    from_disk = ad.read_h5ad(tmp_path / "adata.h5ad")
    ad.tests.helpers.assert_equal(from_disk, in_memory)
