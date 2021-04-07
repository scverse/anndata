from itertools import product

import numpy as np
from numpy import ma
import pandas as pd
import pytest
from scipy import sparse as sp
from scipy.sparse import csr_matrix, issparse

from anndata import AnnData
from anndata.tests.helpers import assert_equal, gen_adata


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
    AnnData(sp.eye(2))
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

    with pytest.raises(ValueError):
        AnnData(np.array([[1, 2], [3, 4]]), dict(TooLong=[1, 2, 3, 4]))

    # init with empty data matrix
    shape = (3, 5)
    adata = AnnData(None, uns=dict(test=np.array((3, 3))), shape=shape)
    assert adata.X is None
    assert adata.shape == shape
    assert "test" in adata.uns


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


def test_create_from_sparse_df():
    s = sp.random(20, 30, density=0.2)
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

    assert adata.obs_names.tolist() == "A B".split()
    assert adata.var_names.tolist() == "a b c".split()

    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]), var=dict(var_names=["a", "b"]))
    assert adata.var_names.tolist() == ["a", "b"]


@pytest.mark.parametrize(
    "names,after",
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
    with pytest.raises(ValueError, match=fr"AnnData expects \.{attr[:3]}\.index\.name"):
        setattr(adata, attr, pd.Index(["x", "y"], name=0))
    assert adata.is_view
    assert getattr(adata, attr).tolist() != ["x", "y"]
    assert getattr(adata, attr).tolist() == getattr(orig, attr).tolist()
    assert_equal(orig, adata, exact=True)


@pytest.mark.parametrize("dim", ["obs", "var"])
def test_setting_dim_index(dim):
    index_attr = f"{dim}_names"
    mapping_attr = f"{dim}m"

    orig = gen_adata((5, 5))
    orig.raw = orig
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


def test_slicing_graphs():
    # Testing for deprecated behaviour of connectivity matrices in .uns["neighbors"]
    with pytest.warns(FutureWarning, match=r".obsp\['connectivities'\]"):
        adata = AnnData(
            np.array([[1, 2], [3, 4], [5, 6]]),
            uns=dict(neighbors=dict(connectivities=sp.csr_matrix(np.ones((3, 3))))),
        )

    adata_sub = adata[[0, 1], :]
    with pytest.warns(FutureWarning):
        assert adata_sub.uns["neighbors"]["connectivities"].shape[0] == 2
        assert adata.uns["neighbors"]["connectivities"].shape[0] == 3
        assert adata_sub.copy().uns["neighbors"]["connectivities"].shape[0] == 2


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

    with pytest.raises(ValueError):
        adata.obs["new4"] = "far too long".split()


def test_delete_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]), dict(o1=[1, 2], o2=[3, 4]))
    assert ["o1", "o2"] == adata.obs_keys()

    del adata.obs["o1"]
    assert ["o2"] == adata.obs_keys()
    assert [3, 4] == adata.obs["o2"].tolist()


def test_set_obs():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.obs = pd.DataFrame(dict(a=[3, 4]))
    assert adata.obs_names.tolist() == [0, 1]

    with pytest.raises(ValueError):
        adata.obs = pd.DataFrame(dict(a=[3, 4, 5]))
        adata.obs = dict(a=[1, 2])


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
        adata1 == adata1
    with pytest.raises(NotImplementedError):
        adata1 == adata2
    with pytest.raises(NotImplementedError):
        adata1 != adata2
    with pytest.raises(NotImplementedError):
        adata1 == 1
    with pytest.raises(NotImplementedError):
        adata1 != 1


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
    adata.raw = adata
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
            adata, adata_dense, lambda x: x.obs_vector(obs_k, layer=layer)
        )

    for obs_k in ["a", "b", "c"]:
        assert_same_op_result(adata, adata_dense, lambda x: x.raw.obs_vector(obs_k))

    for var_k, layer in product(["s1", "s2", "anno2"], [None, "x2"]):
        assert_same_op_result(
            adata, adata_dense, lambda x: x.var_vector(var_k, layer=layer)
        )

    for var_k in ["s1", "s2", "anno2"]:
        assert_same_op_result(adata, adata_dense, lambda x: x.raw.var_vector(var_k))


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
    for attr in "layers var obs obsm varm".split():
        map_sprs = getattr(adata_sparse, attr)
        map_copy = getattr(adata_copy, attr)
        assert map_sprs is not map_copy
        assert_eq_not_id(map_sprs.keys(), map_copy.keys())
        for key in map_sprs.keys():
            assert_eq_not_id(map_sprs[key], map_copy[key])
