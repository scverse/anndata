from collections.abc import Hashable
from copy import deepcopy
from itertools import chain, product
from functools import partial

import numpy as np
from numpy import ma
import pandas as pd
from pandas.api.types import is_categorical_dtype
import pytest
from scipy import sparse
from boltons.iterutils import research, remap, default_exit


from anndata import AnnData, Raw
from anndata._core.index import _subset
from anndata.tests import helpers
from anndata.tests.helpers import asarray, assert_equal, gen_adata


@pytest.fixture(
    params=[asarray, sparse.csr_matrix, sparse.csc_matrix],
    ids=["np_array", "scipy_csr", "scipy_csc"],
)
def array_type(request):
    return request.param


@pytest.fixture(params=["inner", "outer"])
def join_type(request):
    return request.param


@pytest.fixture(params=[0, np.nan, np.pi])
def fill_val(request):
    return request.param


def fix_known_differences(orig, result):
    """
    Helper function for reducing anndata's to only the elements we expect to be
    equivalent after concatenation.

    Only for the case where orig is the ground truth result of what concatenation should be.
    """
    orig = orig.copy()
    result = result.copy()

    result.obs.drop(columns=["batch"], inplace=True)
    result.strings_to_categoricals()  # Should this be implicit in concatenation?

    # TODO
    # * merge varm, varp similar to uns
    # * merge obsp, but some information should be lost
    del orig.varm
    del orig.varp
    del orig.obsp  # TODO

    # Possibly need to fix this, ordered categoricals lose orderedness
    for k, dtype in orig.obs.dtypes.items():
        if is_categorical_dtype(dtype) and dtype.ordered:
            result.obs[k] = result.obs[k].astype(dtype)

    return orig, result


def test_concatenate_roundtrip(join_type, array_type):
    adata = gen_adata((100, 10), X_type=array_type)

    remaining = adata.obs_names
    subsets = []
    while len(remaining) > 0:
        n = min(len(remaining), np.random.choice(50))
        subset_idx = np.random.choice(remaining, n, replace=False)
        subsets.append(adata[subset_idx])
        remaining = remaining.difference(subset_idx)

    result = subsets[0].concatenate(
        subsets[1:], join=join_type, uns_merge="same", index_unique=None
    )

    # Correcting for known differences
    orig, result = fix_known_differences(adata, result)

    assert_equal(result[orig.obs_names].copy(), orig)


def test_concatenate_dense():
    # dense data
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2, X_3=X3),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2, X_3=X3),
        layers={"Xs": X2},
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"], annoB=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2),
        layers=dict(Xs=X3),
    )

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    X_combined = [[2, 3], [5, 6], [3, 2], [6, 5], [3, 2], [6, 5]]
    assert adata.X.astype(int).tolist() == X_combined
    assert adata.layers["Xs"].astype(int).tolist() == X_combined
    assert adata.obs_keys() == ["anno1", "anno2", "batch"]
    assert adata.var_keys() == ["annoA-0", "annoA-1", "annoB-2"]
    assert adata.var.values.tolist() == [[1, 2, 2], [2, 1, 1]]
    assert adata.obsm_keys() == ["X_1", "X_2"]
    assert adata.obsm["X_1"].tolist() == np.concatenate([X1, X1, X1]).tolist()

    # with batch_key and batch_categories
    adata = adata1.concatenate(adata2, adata3, batch_key="batch1")
    assert adata.obs_keys() == ["anno1", "anno2", "batch1"]
    adata = adata1.concatenate(adata2, adata3, batch_categories=["a1", "a2", "a3"])
    assert adata.obs["batch"].cat.categories.tolist() == ["a1", "a2", "a3"]
    assert adata.var_names.tolist() == ["b", "c"]

    # outer join
    adata = adata1.concatenate(adata2, adata3, join="outer")

    X_ref = np.array(
        [
            [1.0, 2.0, 3.0, np.nan],
            [4.0, 5.0, 6.0, np.nan],
            [np.nan, 3.0, 2.0, 1.0],
            [np.nan, 6.0, 5.0, 4.0],
            [np.nan, 3.0, 2.0, 1.0],
            [np.nan, 6.0, 5.0, 4.0],
        ]
    )
    np.testing.assert_equal(adata.X, X_ref)
    var_ma = ma.masked_invalid(adata.var.values.tolist())
    var_ma_ref = ma.masked_invalid(
        np.array(
            [
                [0.0, np.nan, np.nan],
                [1.0, 2.0, 2.0],
                [2.0, 1.0, 1.0],
                [np.nan, 0.0, 0.0],
            ]
        )
    )
    assert np.array_equal(var_ma.mask, var_ma_ref.mask)
    assert np.allclose(var_ma.compressed(), var_ma_ref.compressed())


def test_concatenate_layers(array_type, join_type):
    adatas = []
    for _ in range(5):
        a = array_type(sparse.random(100, 200, format="csr"))
        adatas.append(AnnData(X=a, layers={"a": a}))

    merged = adatas[0].concatenate(adatas[1:], join=join_type)
    assert_equal(merged.X, merged.layers["a"])


@pytest.fixture
def obsm_adatas():
    def gen_index(n):
        return [f"cell{i}" for i in range(n)]

    return [
        AnnData(
            X=sparse.csr_matrix((3, 5)),
            obs=pd.DataFrame(index=gen_index(3)),
            obsm={
                "dense": np.arange(6).reshape(3, 2),
                "sparse": sparse.csr_matrix(np.arange(6).reshape(3, 2)),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(3),
                        "b": list("abc"),
                        "c": pd.Categorical(list("aab")),
                    },
                    index=gen_index(3),
                ),
            },
        ),
        AnnData(
            X=sparse.csr_matrix((4, 10)),
            obs=pd.DataFrame(index=gen_index(4)),
            obsm={
                "dense": np.arange(12).reshape(4, 3),
                "df": pd.DataFrame({"a": np.arange(3, 7),}, index=gen_index(4),),
            },
        ),
        AnnData(
            X=sparse.csr_matrix((2, 100)),
            obs=pd.DataFrame(index=gen_index(2)),
            obsm={
                "sparse": np.arange(8).reshape(2, 4),
                "dense": np.arange(4, 8).reshape(2, 2),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(7, 9),
                        "b": list("cd"),
                        "c": pd.Categorical(list("ab")),
                    },
                    index=gen_index(2),
                ),
            },
        ),
    ]


def test_concatenate_obsm_inner(obsm_adatas):
    adata = obsm_adatas[0].concatenate(obsm_adatas[1:], join="inner")

    assert set(adata.obsm.keys()) == {"dense", "df"}
    assert adata.obsm["dense"].shape == (9, 2)
    assert adata.obsm["dense"].tolist() == [
        [0, 1],
        [2, 3],
        [4, 5],
        [0, 1],
        [3, 4],
        [6, 7],
        [9, 10],
        [4, 5],
        [6, 7],
    ]

    assert adata.obsm["df"].columns == ["a"]
    assert adata.obsm["df"]["a"].tolist() == list(range(9))
    # fmt: off
    true_df = (
        pd.concat([a.obsm["df"] for a in obsm_adatas], join="inner")
        .reset_index(drop=True)
    )
    # fmt: on
    cur_df = adata.obsm["df"].reset_index(drop=True)
    pd.testing.assert_frame_equal(true_df, cur_df)


def test_concatenate_obsm_outer(obsm_adatas, fill_val):
    outer = obsm_adatas[0].concatenate(
        obsm_adatas[1:], join="outer", fill_value=fill_val
    )

    inner = obsm_adatas[0].concatenate(obsm_adatas[1:], join="inner")
    for k, inner_v in inner.obsm.items():
        assert np.array_equal(
            _subset(outer.obsm[k], (slice(None), slice(None, inner_v.shape[1]))),
            inner_v,
        )

    assert set(outer.obsm.keys()) == {"dense", "df", "sparse"}
    assert isinstance(outer.obsm["dense"], np.ndarray)
    np.testing.assert_equal(
        outer.obsm["dense"],
        np.array(
            [
                [0, 1, fill_val],
                [2, 3, fill_val],
                [4, 5, fill_val],
                [0, 1, 2],
                [3, 4, 5],
                [6, 7, 8],
                [9, 10, 11],
                [4, 5, fill_val],
                [6, 7, fill_val],
            ]
        ),
    )

    assert isinstance(outer.obsm["sparse"], sparse.spmatrix)
    np.testing.assert_equal(
        outer.obsm["sparse"].toarray(),
        np.array(
            [
                [0, 1, fill_val, fill_val],
                [2, 3, fill_val, fill_val],
                [4, 5, fill_val, fill_val],
                [fill_val, fill_val, fill_val, fill_val],
                [fill_val, fill_val, fill_val, fill_val],
                [fill_val, fill_val, fill_val, fill_val],
                [fill_val, fill_val, fill_val, fill_val],
                [0, 1, 2, 3],
                [4, 5, 6, 7],
            ]
        ),
    )

    # fmt: off
    true_df = (
        pd.concat([a.obsm["df"] for a in obsm_adatas], join="outer")
        .reset_index(drop=True)
    )
    # fmt: on
    cur_df = outer.obsm["df"].reset_index(drop=True)
    pd.testing.assert_frame_equal(true_df, cur_df)


def test_concatenate_layers_misaligned(array_type, join_type):
    adatas = []
    for _ in range(5):
        a = array_type(sparse.random(100, 200, format="csr"))
        adata = AnnData(X=a, layers={"a": a})
        adatas.append(
            adata[:, np.random.choice(adata.var_names, 150, replace=False)].copy()
        )

    merged = adatas[0].concatenate(adatas[1:], join=join_type)
    assert_equal(merged.X, merged.layers["a"])


def test_concatenate_layers_outer(array_type, fill_val):
    # Testing that issue #368 is fixed
    a = AnnData(
        X=np.ones((10, 20)),
        layers={"a": array_type(sparse.random(10, 20, format="csr"))},
    )
    b = AnnData(X=np.ones((10, 20)))

    c = a.concatenate(b, join="outer", fill_value=fill_val, batch_categories=["a", "b"])

    np.testing.assert_array_equal(
        asarray(c[c.obs["batch"] == "b"].layers["a"]), fill_val
    )


def test_concatenate_fill_value(fill_val):
    def get_obs_els(adata):
        return {
            "X": adata.X,
            **{f"layer_{k}": adata.layers[k] for k in adata.layers},
            **{f"obsm_{k}": adata.obsm[k] for k in adata.obsm},
        }

    adata1 = gen_adata((10, 10))
    adata1.obsm = {
        k: v for k, v in adata1.obsm.items() if not isinstance(v, pd.DataFrame)
    }
    adata2 = gen_adata((10, 5))
    adata2.obsm = {
        k: v[:, : v.shape[1] // 2]
        for k, v in adata2.obsm.items()
        if not isinstance(v, pd.DataFrame)
    }
    adata3 = gen_adata((7, 3))
    adata3.obsm = {
        k: v[:, : v.shape[1] // 3]
        for k, v in adata3.obsm.items()
        if not isinstance(v, pd.DataFrame)
    }
    joined = adata1.concatenate([adata2, adata3], join="outer", fill_value=fill_val)

    ptr = 0
    for orig in [adata1, adata2, adata3]:
        cur = joined[ptr : ptr + orig.n_obs]
        cur_els = get_obs_els(cur)
        orig_els = get_obs_els(orig)
        for k, cur_v in cur_els.items():
            orig_v = orig_els.get(k, sparse.csr_matrix((orig.n_obs, 0)))
            assert_equal(cur_v[:, : orig_v.shape[1]], orig_v)
            np.testing.assert_equal(asarray(cur_v[:, orig_v.shape[1] :]), fill_val)
        ptr += orig.n_obs


def test_concatenate_dense_duplicates():
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    # inner join duplicates
    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoC=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.0],
        ),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoC=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.0],
        ),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.1],
        ),
    )

    adata = adata1.concatenate(adata2, adata3)
    assert adata.var_keys() == [
        "annoA",
        "annoB",
        "annoC-0",
        "annoD-0",
        "annoC-1",
        "annoD-1",
        "annoD-2",
    ]


def test_concatenate_sparse():
    # sparse data
    from scipy.sparse import csr_matrix

    X1 = csr_matrix([[0, 2, 3], [0, 5, 6]])
    X2 = csr_matrix([[0, 2, 3], [0, 5, 6]])
    X3 = csr_matrix([[1, 2, 0], [0, 5, 6]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"]),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"]),
        layers=dict(Xs=X2),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s5", "s6"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"]),
        layers=dict(Xs=X3),
    )

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    X_combined = [[2, 3], [5, 6], [3, 2], [6, 5], [0, 2], [6, 5]]
    assert adata.X.toarray().astype(int).tolist() == X_combined
    assert adata.layers["Xs"].toarray().astype(int).tolist() == X_combined

    # outer join
    adata = adata1.concatenate(adata2, adata3, join="outer")
    assert adata.X.toarray().tolist() == [
        [0.0, 2.0, 3.0, 0.0],
        [0.0, 5.0, 6.0, 0.0],
        [0.0, 3.0, 2.0, 0.0],
        [0.0, 6.0, 5.0, 0.0],
        [0.0, 0.0, 2.0, 1.0],
        [0.0, 6.0, 5.0, 0.0],
    ]


def test_concatenate_mixed():
    X1 = sparse.csr_matrix(np.array([[1, 2, 0], [4, 0, 6], [0, 0, 9]]))
    X2 = sparse.csr_matrix(np.array([[0, 2, 3], [4, 0, 0], [7, 0, 9]]))
    X3 = sparse.csr_matrix(np.array([[1, 0, 3], [0, 0, 6], [0, 8, 0]]))
    X4 = np.array([[0, 2, 3], [4, 0, 0], [7, 0, 9]])
    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2", "s3"], anno1=["c1", "c2", "c3"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        layers=dict(counts=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s4", "s5", "s6"], anno1=["c3", "c4", "c5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(counts=X4),  # sic
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s7", "s8", "s9"], anno2=["d3", "d4", "d5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 2, 3], annoB=[0, 1, 2]),
        layers=dict(counts=X3),
    )
    adata4 = AnnData(
        X4,
        dict(obs_names=["s4", "s5", "s6"], anno1=["c3", "c4", "c5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(counts=X2),  # sic
    )

    adata_all = AnnData.concatenate(adata1, adata2, adata3, adata4)
    assert isinstance(adata_all.X, sparse.csr_matrix)
    assert isinstance(adata_all.layers["counts"], sparse.csr_matrix)


def test_concatenate_with_raw():
    # dense data
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    X4 = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(Xs=X2),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"], annoB=[0, 1, 2]),
        layers=dict(Xs=X3),
    )

    adata4 = AnnData(
        X4,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c", "z"], annoA=[0, 1, 2, 3]),
        layers=dict(Xs=X4),
    )

    adata1.raw = adata1
    adata2.raw = adata2
    adata3.raw = adata3

    adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == {"b", "c"}
    assert_equal(adata_all.raw.to_adata().obs, adata_all.obs)
    assert np.array_equal(adata_all.raw.X, adata_all.X)

    adata_all = AnnData.concatenate(adata1, adata2, adata3, join="outer")
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == set("abcd")
    assert_equal(adata_all.raw.to_adata().obs, adata_all.obs)
    assert np.array_equal(np.nan_to_num(adata_all.raw.X), np.nan_to_num(adata_all.X))

    adata3.raw = adata4
    adata_all = AnnData.concatenate(adata1, adata2, adata3, join="outer")
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == set("abcdz")
    assert set(adata_all.var_names) == set("abcd")
    assert not np.array_equal(
        np.nan_to_num(adata_all.raw.X), np.nan_to_num(adata_all.X)
    )

    del adata3.raw
    with pytest.warns(
        UserWarning,
        match=(
            "Only some AnnData objects have `.raw` attribute, "
            "not concatenating `.raw` attributes."
        ),
    ):
        adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert adata_all.raw is None

    del adata1.raw
    del adata2.raw
    assert all(_adata.raw is None for _adata in (adata1, adata2, adata3))
    adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert adata_all.raw is None


def test_merge_unique():
    from anndata._core.merge import merge_unique

    # Simple cases
    assert merge_unique([{"a": "b"}, {"a": "b"}]) == {"a": "b"}
    assert merge_unique([{"a": {"b": "c"}}, {"a": {"b": "c"}}]) == {"a": {"b": "c"}}
    assert merge_unique([{"a": {"b": "c"}}, {"a": {"b": "d"}}]) == {}
    assert merge_unique([{"a": {"b": "c", "d": "e"}}, {"a": {"b": "c", "d": "f"}}]) == {
        "a": {"b": "c"}
    }

    assert merge_unique(
        [{"a": {"b": {"c": {"d": "e"}}}}, {"a": {"b": {"c": {"d": "e"}}}}]
    ) == {"a": {"b": {"c": {"d": "e"}}}}
    assert (
        merge_unique(
            [
                {"a": {"b": {"c": {"d": "e"}}}},
                {"a": {"b": {"c": {"d": "f"}}}},
                {"a": {"b": {"c": {"d": "e"}}}},
            ]
        )
        == {}
    )

    assert merge_unique([{"a": 1}, {"b": 2}]) == {"a": 1, "b": 2}
    assert merge_unique([{"a": 1}, {"b": 2}, {"a": 1, "b": {"c": 2, "d": 3}}]) == {
        "a": 1
    }

    # Test equivalency between arrays and lists
    assert list(
        merge_unique([{"a": np.ones(5)}, {"a": list(np.ones(5))}])["a"]
    ) == list(np.ones(5))
    assert merge_unique([{"a": np.ones(5)}, {"a": list(np.ones(4))}]) == {}


def test_merge_same():
    from anndata._core.merge import merge_same

    # Same as unique for a number of cases:
    assert merge_same([{"a": "b"}, {"a": "b"}]) == {"a": "b"}
    assert merge_same([{"a": {"b": "c"}}, {"a": {"b": "c"}}]) == {"a": {"b": "c"}}
    assert merge_same([{"a": {"b": "c"}}, {"a": {"b": "d"}}]) == {}
    assert merge_same([{"a": {"b": "c", "d": "e"}}, {"a": {"b": "c", "d": "f"}}]) == {
        "a": {"b": "c"}
    }

    assert merge_same([{"a": {"b": "c"}, "d": "e"}, {"a": {"b": "c"}, "d": 2}]) == {
        "a": {"b": "c"}
    }
    assert merge_same(
        [{"a": {"b": {"c": {"d": "e"}}}}, {"a": {"b": {"c": {"d": "e"}}}}]
    ) == {"a": {"b": {"c": {"d": "e"}}}}

    assert merge_same([{"a": 1}, {"b": 2}]) == {}
    assert merge_same([{"a": 1}, {"b": 2}, {"a": 1, "b": {"c": 2, "d": 3}}]) == {}

    # Test equivalency between arrays and lists
    assert list(merge_same([{"a": np.ones(5)}, {"a": list(np.ones(5))}])["a"]) == list(
        np.ones(5)
    )


def test_merge_first():
    from anndata._core.merge import merge_first

    assert merge_first([{"a": "b"}, {"a": "b"}]) == {"a": "b"}
    assert merge_first([{"a": {"b": "c"}}, {"a": {"b": "c"}}]) == {"a": {"b": "c"}}
    assert merge_first([{"a": 1}, {"a": 2}]) == {"a": 1}

    assert merge_first([{"a": 1}, {"a": {"b": {"c": {"d": "e"}}}}]) == {"a": 1}
    assert merge_first([{"a": {"b": {"c": {"d": "e"}}}}, {"a": 1}]) == {
        "a": {"b": {"c": {"d": "e"}}}
    }


# Helpers for test_concatenate_uns


def uns_ad(uns):
    return AnnData(np.zeros((10, 10)), uns=uns)


def map_values(mapping, path, key, old_parent, new_parent, new_items):
    ret = default_exit(path, key, old_parent, new_parent, new_items)
    for k, v in ret.items():
        if isinstance(v, Hashable) and v in mapping:
            ret[k] = mapping[v]
    return ret


def permute_nested_values(dicts: "List[dict]", gen_val: "Callable[[int], Any]"):
    """
    This function permutes the values of a nested mapping, for testing that out merge
    method work regardless of the values types.

    Assumes the intial dictionary had integers for values.
    """
    dicts = deepcopy(dicts)
    initial_values = [
        x[1] for x in research(dicts, query=lambda p, k, v: isinstance(v, int))
    ]
    mapping = {k: gen_val(k) for k in initial_values}
    return [remap(d, exit=partial(map_values, mapping)) for d in dicts]


def gen_df(n):
    return helpers.gen_typed_df(n)


def gen_array(n):
    return np.random.randn(n)


def gen_list(n):
    return list(gen_array(n))


def gen_sparse(n):
    return sparse.random(np.random.randint(1, 100), np.random.randint(1, 100))


def gen_something(n):
    options = [gen_df, gen_array, gen_list, gen_sparse]
    return np.random.choice(options)(n)


def gen_concat_params(unss, compat2result):
    value_generators = [
        lambda x: x,
        gen_df,
        gen_array,
        gen_list,
        gen_sparse,
        gen_something,
    ]
    for gen, (mode, result) in product(value_generators, compat2result.items()):
        yield pytest.param(unss, mode, result, gen)


@pytest.mark.parametrize(
    ["unss", "merge_strategy", "result", "value_gen"],
    chain(
        gen_concat_params(
            [{"a": 1}, {"a": 2}],
            {None: {}, "first": {"a": 1}, "unique": {}, "same": {}, "only": {}},
        ),
        gen_concat_params(
            [{"a": 1}, {"b": 2}],
            {
                None: {},
                "first": {"a": 1, "b": 2},
                "unique": {"a": 1, "b": 2},
                "same": {},
                "only": {"a": 1, "b": 2},
            },
        ),
        gen_concat_params(
            [{"a": {"b": 1, "c": {"d": 3}}}, {"a": {"b": 1, "c": {"e": 4}}},],
            {
                None: {},
                "first": {"a": {"b": 1, "c": {"d": 3, "e": 4}}},
                "unique": {"a": {"b": 1, "c": {"d": 3, "e": 4}}},
                "same": {"a": {"b": 1}},
                "only": {"a": {"c": {"d": 3, "e": 4}}},
            },
        ),
        gen_concat_params(
            [{"a": 1}, {"a": 1, "b": 2}, {"a": 1, "b": {"b.a": 1}, "c": 3}, {"d": 4},],
            {
                None: {},
                "first": {"a": 1, "b": 2, "c": 3, "d": 4},
                "unique": {"a": 1, "c": 3, "d": 4},
                "same": {},
                "only": {"c": 3, "d": 4},
            },
        ),
        gen_concat_params(
            [{"a": i} for i in range(15)],
            {None: {}, "first": {"a": 0}, "unique": {}, "same": {}, "only": {},},
        ),
        gen_concat_params(
            [{"a": 1} for i in range(10)] + [{"a": 2}],
            {None: {}, "first": {"a": 1}, "unique": {}, "same": {}, "only": {},},
        ),
    ),
)
def test_concatenate_uns(unss, merge_strategy, result, value_gen):
    """
    Test that concatenation works out for different strategies and sets of values.

    Params
    ------
    unss
        Set of patterns for values in uns.
    compat
        Strategy to use for merging uns.
    result
        Pattern we expect to see for the given input and strategy.
    value_gen
        Maps values in unss and results to another set of values. This is for checking that
        we're comparing values correctly. For example `[{"a": 1}, {"a": 1}]` may get mapped
        to `[{"a": [1, 2, 3]}, {"a": [1, 2, 3]}]`.
    """
    # So we can see what the initial pattern was meant to be
    print(merge_strategy, "\n", unss, "\n", result)
    result, *unss = permute_nested_values([result] + unss, value_gen)
    adatas = [uns_ad(uns) for uns in unss]
    assert_equal(
        adatas[0].concatenate(adatas[1:], uns_merge=merge_strategy).uns,
        result,
        elem_name="uns",
    )


# Leaving out for now. See definition of these values for explanation
# def test_concatenate_uns_types():
#     from anndata._core.merge import UNS_STRATEGIES, UNS_STRATEGIES_TYPE
#     assert set(UNS_STRATEGIES.keys()) == set(UNS_STRATEGIES_TYPE.__args__)
