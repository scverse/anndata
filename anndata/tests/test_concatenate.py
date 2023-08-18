from collections.abc import Hashable
from copy import deepcopy
from itertools import chain, product
from functools import partial, singledispatch
from typing import Any, List, Callable
import warnings

import numpy as np
from numpy import ma
import pandas as pd
import pytest
from scipy import sparse
from boltons.iterutils import research, remap, default_exit


from anndata import AnnData, Raw, concat
from anndata._core.index import _subset
from anndata._core import merge
from anndata.tests import helpers
from anndata.tests.helpers import (
    assert_equal,
    as_dense_dask_array,
    gen_adata,
    GEN_ADATA_DASK_ARGS,
    BASE_MATRIX_PARAMS,
    DASK_MATRIX_PARAMS,
    CUPY_MATRIX_PARAMS,
)
from anndata.utils import asarray
from anndata.compat import DaskArray, AwkArray


@singledispatch
def filled_like(a, fill_value=None):
    raise NotImplementedError()


@filled_like.register(np.ndarray)
def _filled_array_np(a, fill_value=None):
    if fill_value is None:
        fill_value = np.nan
    return np.broadcast_to(fill_value, a.shape)


@filled_like.register(DaskArray)
def _filled_array(a, fill_value=None):
    return as_dense_dask_array(_filled_array_np(a, fill_value))


@filled_like.register(sparse.spmatrix)
def _filled_sparse(a, fill_value=None):
    if fill_value is None:
        return sparse.csr_matrix(a.shape)
    else:
        return sparse.csr_matrix(np.broadcast_to(fill_value, a.shape))


@filled_like.register(pd.DataFrame)
def _filled_df(a, fill_value=np.nan):
    # dtype from pd.concat can be unintuitive, this returns something close enough
    return a.loc[[], :].reindex(index=a.index, fill_value=fill_value)


def check_filled_like(x, fill_value=None, elem_name=None):
    if fill_value is None:
        assert_equal(x, filled_like(x), elem_name=elem_name)
    else:
        assert_equal(x, filled_like(x, fill_value=fill_value), elem_name=elem_name)


def make_idx_tuple(idx, axis):
    tup = [slice(None), slice(None)]
    tup[axis] = idx
    return tuple(tup)


# Will call func(sparse_matrix) so these types should be sparse compatible
# See array_type if only dense arrays are expected as input.
@pytest.fixture(params=BASE_MATRIX_PARAMS + DASK_MATRIX_PARAMS + CUPY_MATRIX_PARAMS)
def array_type(request):
    return request.param


@pytest.fixture(params=["inner", "outer"])
def join_type(request):
    return request.param


@pytest.fixture(params=[0, np.nan, np.pi])
def fill_val(request):
    return request.param


@pytest.fixture(params=[0, 1])
def axis(request):
    return request.param


@pytest.fixture(params=list(merge.MERGE_STRATEGIES.keys()))
def merge_strategy(request):
    return request.param


def fix_known_differences(orig, result, backwards_compat=True):
    """
    Helper function for reducing anndata's to only the elements we expect to be
    equivalent after concatenation.

    Only for the case where orig is the ground truth result of what concatenation should be.

    If backwards_compat, checks against what `AnnData.concatenate` could do. Otherwise checks for `concat`.
    """
    orig = orig.copy()
    result = result.copy()

    result.strings_to_categoricals()  # Should this be implicit in concatenation?

    # TODO
    # * merge varm, varp similar to uns
    # * merge obsp, but some information should be lost
    del orig.obsp  # TODO

    if backwards_compat:
        del orig.varm
        del orig.varp
        result.obs.drop(columns=["batch"], inplace=True)

    # Possibly need to fix this, ordered categoricals lose orderedness
    for k, dtype in orig.obs.dtypes.items():
        if isinstance(dtype, pd.CategoricalDtype) and dtype.ordered:
            result.obs[k] = result.obs[k].astype(dtype)

    return orig, result


def test_concat_interface_errors():
    adatas = [gen_adata((5, 10)), gen_adata((5, 10))]

    with pytest.raises(ValueError):
        concat(adatas, axis=3)
    with pytest.raises(ValueError):
        concat(adatas, join="not implemented")
    with pytest.raises(ValueError):
        concat([])


@pytest.mark.parametrize(
    ["concat_func", "backwards_compat"],
    [
        (partial(concat, merge="unique"), False),
        (lambda x, **kwargs: x[0].concatenate(x[1:], **kwargs), True),
    ],
)
def test_concatenate_roundtrip(join_type, array_type, concat_func, backwards_compat):
    adata = gen_adata((100, 10), X_type=array_type, **GEN_ADATA_DASK_ARGS)

    remaining = adata.obs_names
    subsets = []
    while len(remaining) > 0:
        n = min(len(remaining), np.random.choice(50))
        subset_idx = np.random.choice(remaining, n, replace=False)
        subsets.append(adata[subset_idx])
        remaining = remaining.difference(subset_idx)

    result = concat_func(subsets, join=join_type, uns_merge="same", index_unique=None)

    # Correcting for known differences
    orig, result = fix_known_differences(
        adata, result, backwards_compat=backwards_compat
    )

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
            obsm=dict(
                dense=np.arange(12).reshape(4, 3),
                df=pd.DataFrame(dict(a=np.arange(3, 7)), index=gen_index(4)),
            ),
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


def test_concat_annot_join(obsm_adatas, join_type):
    adatas = [
        AnnData(sparse.csr_matrix(a.shape), obs=a.obsm["df"], var=a.var)
        for a in obsm_adatas
    ]
    pd.testing.assert_frame_equal(
        concat(adatas, join=join_type).obs,
        pd.concat([a.obs for a in adatas], join=join_type),
    )


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
        k: v
        for k, v in adata1.obsm.items()
        if not isinstance(v, (pd.DataFrame, AwkArray))
    }
    adata2 = gen_adata((10, 5))
    adata2.obsm = {
        k: v[:, : v.shape[1] // 2]
        for k, v in adata2.obsm.items()
        if not isinstance(v, (pd.DataFrame, AwkArray))
    }
    adata3 = gen_adata((7, 3))
    adata3.obsm = {
        k: v[:, : v.shape[1] // 3]
        for k, v in adata3.obsm.items()
        if not isinstance(v, (pd.DataFrame, AwkArray))
    }
    # remove AwkArrays from adata.var, as outer joins are not yet implemented for them
    for tmp_ad in [adata1, adata2, adata3]:
        for k in [k for k, v in tmp_ad.varm.items() if isinstance(v, AwkArray)]:
            del tmp_ad.varm[k]

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


def test_concatenate_awkward(join_type):
    import awkward as ak

    a = ak.Array([[{"a": 1, "b": "foo"}], [{"a": 2, "b": "bar"}, {"a": 3, "b": "baz"}]])
    b = ak.Array(
        [
            [{"a": 4}, {"a": 5}],
            [{"a": 6}],
            [{"a": 7}],
        ]
    )

    adata_a = AnnData(np.zeros((2, 0), dtype=float), obsm={"awk": a})
    adata_b = AnnData(np.zeros((3, 0), dtype=float), obsm={"awk": b})

    if join_type == "inner":
        expected = ak.Array(
            [
                [{"a": 1}],
                [{"a": 2}, {"a": 3}],
                [{"a": 4}, {"a": 5}],
                [{"a": 6}],
                [{"a": 7}],
            ]
        )
    elif join_type == "outer":
        # TODO: This is what we would like to return, but waiting on:
        # * https://github.com/scikit-hep/awkward/issues/2182 and awkward 2.1.0
        # * https://github.com/scikit-hep/awkward/issues/2173
        # expected = ak.Array(
        #     [
        #         [{"a": 1, "b": "foo"}],
        #         [{"a": 2, "b": "bar"}, {"a": 3, "b": "baz"}],
        #         [{"a": 4, "b": None}, {"a": 5, "b": None}],
        #         [{"a": 6, "b": None}],
        #         [{"a": 7, "b": None}],
        #     ]
        # )
        expected = ak.concatenate(
            [  # I don't think I can construct a UnionArray directly
                ak.Array(
                    [
                        [{"a": 1, "b": "foo"}],
                        [{"a": 2, "b": "bar"}, {"a": 3, "b": "baz"}],
                    ]
                ),
                ak.Array(
                    [
                        [{"a": 4}, {"a": 5}],
                        [{"a": 6}],
                        [{"a": 7}],
                    ]
                ),
            ]
        )

    result = concat([adata_a, adata_b], join=join_type).obsm["awk"]

    assert_equal(expected, result)


@pytest.mark.parametrize(
    "other",
    [
        pd.DataFrame({"a": [4, 5, 6], "b": ["foo", "bar", "baz"]}, index=list("cde")),
        np.ones((3, 2)),
        sparse.random(3, 100, format="csr"),
    ],
)
def test_awkward_does_not_mix(join_type, other):
    import awkward as ak

    awk = ak.Array(
        [[{"a": 1, "b": "foo"}], [{"a": 2, "b": "bar"}, {"a": 3, "b": "baz"}]]
    )

    adata_a = AnnData(
        np.zeros((2, 3), dtype=float),
        obs=pd.DataFrame(index=list("ab")),
        obsm={"val": awk},
    )
    adata_b = AnnData(
        np.zeros((3, 3), dtype=float),
        obs=pd.DataFrame(index=list("cde")),
        obsm={"val": other},
    )

    with pytest.raises(
        NotImplementedError,
        match="Cannot concatenate an AwkwardArray with other array types",
    ):
        concat([adata_a, adata_b], join=join_type)


def test_pairwise_concat(axis, array_type):
    dim_sizes = [[100, 200, 50], [50, 50, 50]]
    if axis:
        dim_sizes.reverse()
    Ms, Ns = dim_sizes
    dim = ("obs", "var")[axis]
    alt = ("var", "obs")[axis]
    dim_attr = f"{dim}p"
    alt_attr = f"{alt}p"

    def gen_dim_array(m):
        return array_type(sparse.random(m, m, format="csr", density=0.1))

    adatas = {
        k: AnnData(
            **{
                "X": sparse.csr_matrix((m, n)),
                "obsp": {"arr": gen_dim_array(m)},
                "varp": {"arr": gen_dim_array(n)},
            }
        )
        for k, m, n in zip("abc", Ms, Ns)
    }

    w_pairwise = concat(adatas, axis=axis, label="orig", pairwise=True)
    wo_pairwise = concat(adatas, axis=axis, label="orig", pairwise=False)

    # Check that argument controls whether elements are included
    assert getattr(wo_pairwise, dim_attr) == {}
    assert getattr(w_pairwise, dim_attr) != {}

    # Check values of included elements
    full_inds = np.arange(w_pairwise.shape[axis])
    groups = getattr(w_pairwise, dim).groupby("orig").indices
    for k, inds in groups.items():
        orig_arr = getattr(adatas[k], dim_attr)["arr"]
        full_arr = getattr(w_pairwise, dim_attr)["arr"]

        # Check original values are intact
        assert_equal(orig_arr, _subset(full_arr, (inds, inds)))
        # Check that entries are filled with zeroes
        assert_equal(
            sparse.csr_matrix((len(inds), len(full_inds) - len(inds))),
            _subset(full_arr, (inds, np.setdiff1d(full_inds, inds))),
        )
        assert_equal(
            sparse.csr_matrix((len(full_inds) - len(inds), len(inds))),
            _subset(full_arr, (np.setdiff1d(full_inds, inds), inds)),
        )

    # Check that argument does not affect alternative axis
    assert "arr" in getattr(
        concat(adatas, axis=axis, pairwise=False, merge="first"), alt_attr
    )


def test_nan_merge(axis, join_type, array_type):
    # concat_dim = ("obs", "var")[axis]
    alt_dim = ("var", "obs")[axis]
    mapping_attr = f"{alt_dim}m"
    adata_shape = (20, 10)

    arr = array_type(
        sparse.random(adata_shape[1 - axis], 10, density=0.1, format="csr")
    )
    arr_nan = arr.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=sparse.SparseEfficiencyWarning)
        for _ in range(10):
            arr_nan[
                np.random.choice(arr.shape[0]), np.random.choice(arr.shape[1])
            ] = np.nan

    _data = {"X": sparse.csr_matrix(adata_shape), mapping_attr: {"arr": arr_nan}}
    orig1 = AnnData(**_data)
    orig2 = AnnData(**_data)
    result = concat([orig1, orig2], axis=axis, merge="same")

    assert_equal(getattr(orig1, mapping_attr), getattr(result, mapping_attr))

    orig_nonan = AnnData(
        **{"X": sparse.csr_matrix(adata_shape), mapping_attr: {"arr": arr}}
    )
    result_nonan = concat([orig1, orig_nonan], axis=axis, merge="same")

    assert len(getattr(result_nonan, mapping_attr)) == 0


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
            [
                {"a": {"b": 1, "c": {"d": 3}}},
                {"a": {"b": 1, "c": {"e": 4}}},
            ],
            {
                None: {},
                "first": {"a": {"b": 1, "c": {"d": 3, "e": 4}}},
                "unique": {"a": {"b": 1, "c": {"d": 3, "e": 4}}},
                "same": {"a": {"b": 1}},
                "only": {"a": {"c": {"d": 3, "e": 4}}},
            },
        ),
        gen_concat_params(
            [
                {"a": 1},
                {"a": 1, "b": 2},
                {"a": 1, "b": {"b.a": 1}, "c": 3},
                {"d": 4},
            ],
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
            {None: {}, "first": {"a": 0}, "unique": {}, "same": {}, "only": {}},
        ),
        gen_concat_params(
            [{"a": 1} for i in range(10)] + [{"a": 2}],
            {None: {}, "first": {"a": 1}, "unique": {}, "same": {}, "only": {}},
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
    with pytest.warns(FutureWarning, match=r"concatenate method is deprecated"):
        merged = AnnData.concatenate(*adatas, uns_merge=merge_strategy).uns
    assert_equal(merged, result, elem_name="uns")


def test_transposed_concat(array_type, axis, join_type, merge_strategy, fill_val):
    lhs = gen_adata((10, 10), X_type=array_type, **GEN_ADATA_DASK_ARGS)
    rhs = gen_adata((10, 12), X_type=array_type, **GEN_ADATA_DASK_ARGS)

    a = concat([lhs, rhs], axis=axis, join=join_type, merge=merge_strategy)
    b = concat(
        [lhs.T, rhs.T], axis=abs(axis - 1), join=join_type, merge=merge_strategy
    ).T

    assert_equal(a, b)


def test_batch_key(axis, array_type):
    """Test that concat only adds a label if the key is provided"""

    def get_annot(adata):
        return getattr(adata, ("obs", "var")[axis])

    lhs = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    rhs = gen_adata((10, 12), **GEN_ADATA_DASK_ARGS)

    # There is probably a prettier way to do this
    annot = get_annot(concat([lhs, rhs], axis=axis))
    assert (
        list(
            annot.columns.difference(
                get_annot(lhs).columns.union(get_annot(rhs).columns)
            )
        )
        == []
    )

    batch_annot = get_annot(concat([lhs, rhs], axis=axis, label="batch"))
    assert list(
        batch_annot.columns.difference(
            get_annot(lhs).columns.union(get_annot(rhs).columns)
        )
    ) == ["batch"]


def test_concat_categories_from_mapping():
    mapping = {
        "a": gen_adata((10, 10)),
        "b": gen_adata((10, 10)),
    }
    keys = list(mapping.keys())
    adatas = list(mapping.values())

    mapping_call = partial(concat, mapping)
    iter_call = partial(concat, adatas, keys=keys)

    assert_equal(mapping_call(), iter_call())
    assert_equal(mapping_call(label="batch"), iter_call(label="batch"))
    assert_equal(mapping_call(index_unique="-"), iter_call(index_unique="-"))
    assert_equal(
        mapping_call(label="group", index_unique="+"),
        iter_call(label="group", index_unique="+"),
    )


def test_concat_categories_maintain_dtype():
    a = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat": pd.Categorical(list("aabcc")),
                "cat_ordered": pd.Categorical(list("aabcc"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5)],
        ),
    )
    b = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat": pd.Categorical(list("bccdd")),
                "cat_ordered": pd.Categorical(list("bccdd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5, 10)],
        ),
    )
    c = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("bccdd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5, 10)],
        ),
    )

    result = concat({"a": a, "b": b, "c": c}, join="outer")

    assert isinstance(
        result.obs["cat"].dtype, pd.CategoricalDtype
    ), f"Was {result.obs['cat'].dtype}"
    assert pd.api.types.is_string_dtype(result.obs["cat_ordered"])


def test_concat_ordered_categoricals_retained():
    a = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("aabcd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5)],
        ),
    )
    b = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("abcdd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5, 10)],
        ),
    )

    c = concat([a, b])

    assert isinstance(c.obs["cat_ordered"].dtype, pd.CategoricalDtype)
    assert c.obs["cat_ordered"].cat.ordered


def test_bool_promotion():
    np_bool = AnnData(
        np.ones((5, 1)),
        obs=pd.DataFrame({"bool": [True] * 5}, index=[f"cell{i:02}" for i in range(5)]),
    )
    missing = AnnData(
        np.ones((5, 1)),
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(5, 10)]),
    )
    result = concat({"np_bool": np_bool, "b": missing}, join="outer", label="batch")

    assert pd.api.types.is_bool_dtype(result.obs["bool"])
    assert pd.isnull(result.obs.loc[result.obs["batch"] == "missing", "bool"]).all()

    # Check that promotion doesn't occur if it doesn't need to:
    np_bool_2 = AnnData(
        np.ones((5, 1)),
        obs=pd.DataFrame(
            {"bool": [True] * 5}, index=[f"cell{i:02}" for i in range(5, 10)]
        ),
    )
    result = concat(
        {"np_bool": np_bool, "np_bool_2": np_bool_2}, join="outer", label="batch"
    )

    assert result.obs["bool"].dtype == np.dtype(bool)


def test_concat_names(axis):
    def get_annot(adata):
        return getattr(adata, ("obs", "var")[axis])

    lhs = gen_adata((10, 10))
    rhs = gen_adata((10, 10))

    assert not get_annot(concat([lhs, rhs], axis=axis)).index.is_unique
    assert get_annot(concat([lhs, rhs], axis=axis, index_unique="-")).index.is_unique


def axis_labels(adata, axis):
    return (adata.obs_names, adata.var_names)[axis]


def expected_shape(a, b, axis, join):
    labels = partial(axis_labels, axis=abs(axis - 1))
    shape = [None, None]

    shape[axis] = a.shape[axis] + b.shape[axis]
    if join == "inner":
        shape[abs(axis - 1)] = len(labels(a).intersection(labels(b)))
    elif join == "outer":
        shape[abs(axis - 1)] = len(labels(a).union(labels(b)))
    else:
        raise ValueError()

    return tuple(shape)


@pytest.mark.parametrize(
    "shape", [pytest.param((8, 0), id="no_var"), pytest.param((0, 10), id="no_obs")]
)
def test_concat_size_0_dim(axis, join_type, merge_strategy, shape):
    # https://github.com/scverse/anndata/issues/526
    a = gen_adata((5, 7))
    b = gen_adata(shape)
    alt_axis = 1 - axis
    dim = ("obs", "var")[axis]

    expected_size = expected_shape(a, b, axis=axis, join=join_type)
    result = concat(
        {"a": a, "b": b},
        axis=axis,
        join=join_type,
        merge=merge_strategy,
        pairwise=True,
        index_unique="-",
    )
    assert result.shape == expected_size

    if join_type == "outer":
        # Check new entries along axis of concatenation
        axis_new_inds = axis_labels(result, axis).str.endswith("-b")
        altaxis_new_inds = ~axis_labels(result, alt_axis).isin(axis_labels(a, alt_axis))
        axis_idx = make_idx_tuple(axis_new_inds, axis)
        altaxis_idx = make_idx_tuple(altaxis_new_inds, 1 - axis)

        check_filled_like(result.X[axis_idx], elem_name="X")
        check_filled_like(result.X[altaxis_idx], elem_name="X")
        for k, elem in getattr(result, "layers").items():
            check_filled_like(elem[axis_idx], elem_name=f"layers/{k}")
            check_filled_like(elem[altaxis_idx], elem_name=f"layers/{k}")

        if shape[axis] > 0:
            b_result = result[axis_idx].copy()
            mapping_elem = f"{dim}m"
            setattr(b_result, f"{dim}_names", getattr(b, f"{dim}_names"))
            for k, result_elem in getattr(b_result, mapping_elem).items():
                elem_name = f"{mapping_elem}/{k}"
                # pd.concat can have unintuitive return types. is similar to numpy promotion
                if isinstance(result_elem, pd.DataFrame):
                    assert_equal(
                        getattr(b, mapping_elem)[k].astype(object),
                        result_elem.astype(object),
                        elem_name=elem_name,
                    )
                else:
                    assert_equal(
                        getattr(b, mapping_elem)[k],
                        result_elem,
                        elem_name=elem_name,
                    )


@pytest.mark.parametrize("elem", ["sparse", "array", "df", "da"])
def test_concat_outer_aligned_mapping(elem):
    a = gen_adata((5, 5), **GEN_ADATA_DASK_ARGS)
    b = gen_adata((3, 5), **GEN_ADATA_DASK_ARGS)
    del b.obsm[elem]

    concated = concat({"a": a, "b": b}, join="outer", label="group")
    result = concated.obsm[elem][concated.obs["group"] == "b"]

    check_filled_like(result, elem_name=f"obsm/{elem}")


def test_concatenate_size_0_dim():
    # https://github.com/scverse/anndata/issues/526

    a = gen_adata((5, 10))
    b = gen_adata((5, 0))

    # Mostly testing that this doesn't error
    a.concatenate([b]).shape == (10, 0)
    b.concatenate([a]).shape == (10, 0)


def test_concat_null_X():
    adatas_orig = {k: gen_adata((20, 10)) for k in list("abc")}
    adatas_no_X = {}
    for k, v in adatas_orig.items():
        v = v.copy()
        del v.X
        adatas_no_X[k] = v

    orig = concat(adatas_orig, index_unique="-")
    no_X = concat(adatas_no_X, index_unique="-")
    del orig.X

    assert_equal(no_X, orig)


# https://github.com/scverse/ehrapy/issues/151#issuecomment-1016753744
def test_concat_X_dtype():
    adatas_orig = {k: AnnData(np.ones((20, 10), dtype=np.int8)) for k in list("abc")}
    for adata in adatas_orig.values():
        adata.raw = AnnData(np.ones((20, 30), dtype=np.float64))

    result = concat(adatas_orig, index_unique="-")

    assert result.X.dtype == np.int8
    assert result.raw.X.dtype == np.float64


# Leaving out for now. See definition of these values for explanation
# def test_concatenate_uns_types():
#     from anndata._core.merge import UNS_STRATEGIES, UNS_STRATEGIES_TYPE
#     assert set(UNS_STRATEGIES.keys()) == set(UNS_STRATEGIES_TYPE.__args__)


# Tests how dask plays with other types on concatenation.
def test_concat_different_types_dask(merge_strategy, array_type):
    from scipy import sparse
    import anndata as ad
    import dask.array as da

    varm_array = sparse.random(5, 20, density=0.5, format="csr")

    ad1 = ad.AnnData(X=np.ones((5, 5)), varm={"a": varm_array})
    ad1_other = ad.AnnData(X=np.ones((5, 5)), varm={"a": array_type(varm_array)})
    ad2 = ad.AnnData(X=np.zeros((5, 5)), varm={"a": da.ones(5, 20)})

    result1 = ad.concat([ad1, ad2], merge=merge_strategy)
    target1 = ad.concat([ad1_other, ad2], merge=merge_strategy)
    result2 = ad.concat([ad2, ad1], merge=merge_strategy)
    target2 = ad.concat([ad2, ad1_other], merge=merge_strategy)

    assert_equal(result1, target1)
    assert_equal(result2, target2)


def test_outer_concat_with_missing_value_for_df():
    # https://github.com/scverse/anndata/issues/901
    # TODO: Extend this test to cover all cases of missing values
    # TODO: Check values
    a_idx = ["a", "b", "c", "d", "e"]
    b_idx = ["f", "g", "h", "i", "j", "k", "l", "m"]
    a = AnnData(
        np.ones((5, 5)),
        obs=pd.DataFrame(index=a_idx),
    )
    b = AnnData(
        np.zeros((8, 9)),
        obs=pd.DataFrame(index=b_idx),
        obsm={"df": pd.DataFrame({"col": np.arange(8)}, index=b_idx)},
    )

    concat([a, b], join="outer")


def test_outer_concat_outputs_nullable_bool_writable(tmp_path):
    a = gen_adata((5, 5), obsm_types=(pd.DataFrame,))
    b = gen_adata((3, 5), obsm_types=(pd.DataFrame,))

    del b.obsm["df"]

    adatas = concat({"a": a, "b": b}, join="outer", label="group")
    adatas.write(tmp_path / "test.h5ad")


def test_concat_duplicated_columns(join_type):
    # https://github.com/scverse/anndata/issues/483
    a = AnnData(
        obs=pd.DataFrame(
            np.ones((5, 2)), columns=["a", "a"], index=[str(x) for x in range(5)]
        )
    )
    b = AnnData(
        obs=pd.DataFrame(
            np.ones((5, 1)), columns=["a"], index=[str(x) for x in range(5, 10)]
        )
    )

    with pytest.raises(pd.errors.InvalidIndexError, match=r"'a'"):
        concat([a, b], join=join_type)
