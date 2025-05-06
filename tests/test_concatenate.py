from __future__ import annotations

import warnings
from collections.abc import Hashable
from copy import deepcopy
from functools import partial, singledispatch
from itertools import chain, permutations, product
from operator import attrgetter
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import scipy
import xarray as xr
from boltons.iterutils import default_exit, remap, research
from numpy import ma
from packaging.version import Version
from scipy import sparse

from anndata import AnnData, Raw, concat
from anndata._core import merge
from anndata._core.index import _subset
from anndata.compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    CupySparseMatrix,
    DaskArray,
    XDataset,
)
from anndata.tests import helpers
from anndata.tests.helpers import (
    BASE_MATRIX_PARAMS,
    CUPY_MATRIX_PARAMS,
    DASK_MATRIX_PARAMS,
    DEFAULT_COL_TYPES,
    GEN_ADATA_DASK_ARGS,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
    gen_vstr_recarray,
)
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

mark_legacy_concatenate = pytest.mark.filterwarnings(
    r"ignore:.*AnnData\.concatenate is deprecated:FutureWarning"
)


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


@filled_like.register(CSMatrix)
def _filled_sparse(a, fill_value=None):
    if fill_value is None:
        return sparse.csr_matrix(a.shape)
    else:
        return sparse.csr_matrix(np.broadcast_to(fill_value, a.shape))


@filled_like.register(CSArray)
def _filled_sparse_array(a, fill_value=None):
    return sparse.csr_array(filled_like(sparse.csr_matrix(a)))


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


@pytest.fixture(params=BASE_MATRIX_PARAMS + DASK_MATRIX_PARAMS)
def cpu_array_type(request):
    return request.param


@pytest.fixture(params=["inner", "outer"])
def join_type(request):
    return request.param


@pytest.fixture(params=[0, np.nan, np.pi])
def fill_val(request):
    return request.param


@pytest.fixture(params=["obs", "var"])
def axis_name(request) -> Literal["obs", "var"]:
    return request.param


@pytest.fixture(params=list(merge.MERGE_STRATEGIES.keys()))
def merge_strategy(request):
    return request.param


def fix_known_differences(
    orig: AnnData, result: AnnData, *, backwards_compat: bool = True
):
    """
    Helper function for reducing anndata's to only the elements we expect to be
    equivalent after concatenation.

    Only for the case where orig is the ground truth result of what concatenation should be.

    If backwards_compat, checks against what `AnnData.concatenate` could do. Otherwise checks for `concat`.
    """
    orig = orig.copy()
    result = result.copy()

    if backwards_compat:
        del orig.varm
        del orig.varp
        if isinstance(result.obs, XDataset):
            result.obs = result.obs.drop_vars(["batch"])
        else:
            result.obs.drop(columns=["batch"], inplace=True)

    for attrname in ("obs", "var"):
        if isinstance(getattr(result, attrname), XDataset):
            for adata in (orig, result):
                df = getattr(adata, attrname).to_dataframe()
                df.index.name = "index"
                setattr(adata, attrname, df)
            resattr = getattr(result, attrname)
            origattr = getattr(orig, attrname)
            for colname, col in resattr.items():
                # concatenation of XDatasets happens via Dask arrays and those don't know about Pandas Extension arrays
                # so categoricals and nullable arrays are all converted to other dtypes
                if col.dtype != origattr[
                    colname
                ].dtype and pd.api.types.is_extension_array_dtype(
                    origattr[colname].dtype
                ):
                    resattr[colname] = col.astype(origattr[colname].dtype)

    result.strings_to_categoricals()  # Should this be implicit in concatenation?

    # TODO
    # * merge varm, varp similar to uns
    # * merge obsp, but some information should be lost
    del orig.obsp  # TODO

    # Possibly need to fix this, ordered categoricals lose orderedness
    for get_df in [lambda k: k.obs, lambda k: k.obsm["df"]]:
        str_to_df_converted = get_df(result)
        for k, dtype in get_df(orig).dtypes.items():
            if isinstance(dtype, pd.CategoricalDtype) and dtype.ordered:
                str_to_df_converted[k] = str_to_df_converted[k].astype(dtype)

    return orig, result


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset"), [(False, False), (True, True)]
)
def test_concat_interface_errors(obs_xdataset, var_xdataset):
    adatas = [
        gen_adata((5, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset),
        gen_adata((5, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset),
    ]

    with pytest.raises(ValueError, match="`axis` must be.*0, 1, 'obs', or 'var'"):
        concat(adatas, axis=3)
    with pytest.raises(ValueError, match="'inner' or 'outer'"):
        concat(adatas, join="not implemented")
    with pytest.raises(ValueError, match="No objects to concatenate"):
        concat([])


@mark_legacy_concatenate
@pytest.mark.parametrize(
    ("concat_func", "backwards_compat"),
    [
        (partial(concat, merge="unique"), False),
        (lambda x, **kwargs: x[0].concatenate(x[1:], **kwargs), True),
    ],
)
@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concatenate_roundtrip(
    join_type,
    array_type,
    concat_func,
    backwards_compat,
    obs_xdataset,
    var_xdataset,
    force_lazy,
):
    if backwards_compat and force_lazy:
        pytest.skip("unsupported")
    adata = gen_adata(
        (100, 10),
        X_type=array_type,
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )

    remaining = adata.obs_names
    subsets = []
    while len(remaining) > 0:
        n = min(len(remaining), np.random.choice(50))
        subset_idx = np.random.choice(remaining, n, replace=False)
        subsets.append(adata[subset_idx])
        remaining = remaining.difference(subset_idx)

    if (
        backwards_compat
        and (obs_xdataset or var_xdataset)
        and Version(xr.__version__) < Version("2025.4.0")
    ):
        pytest.xfail("https://github.com/pydata/xarray/issues/10218")
    result = concat_func(subsets, join=join_type, uns_merge="same", index_unique=None)
    if backwards_compat and var_xdataset:
        result.var = xr.Dataset.from_dataframe(
            result.var
        )  # backwards compat always returns a dataframe

    # Correcting for known differences
    orig, result = fix_known_differences(
        adata, result, backwards_compat=backwards_compat
    )

    assert_equal(result[orig.obs_names].copy(), orig)
    base_type = type(orig.X)
    if sparse.issparse(orig.X):
        base_type = CSArray if isinstance(orig.X, CSArray) else CSMatrix
    if isinstance(orig.X, CupySparseMatrix):
        base_type = CupySparseMatrix
    assert isinstance(result.X, base_type)


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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

    assert isinstance(outer.obsm["sparse"], CSMatrix)
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


@pytest.mark.parametrize(
    ("axis", "axis_name"),
    [("obs", 0), ("var", 1)],
)
def test_concat_axis_param(axis, axis_name):
    a, b = gen_adata((10, 10)), gen_adata((10, 10))
    assert_equal(concat([a, b], axis=axis), concat([a, b], axis=axis_name))


def test_concat_annot_join(obsm_adatas, join_type):
    adatas = [
        AnnData(sparse.csr_matrix(a.shape), obs=a.obsm["df"], var=a.var)
        for a in obsm_adatas
    ]
    pd.testing.assert_frame_equal(
        concat(adatas, join=join_type).obs,
        pd.concat([a.obs for a in adatas], join=join_type),
    )


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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
        if not isinstance(v, pd.DataFrame | AwkArray | XDataset)
    }
    adata2 = gen_adata((10, 5))
    adata2.obsm = {
        k: v[:, : v.shape[1] // 2]
        for k, v in adata2.obsm.items()
        if not isinstance(v, pd.DataFrame | AwkArray | XDataset)
    }
    adata3 = gen_adata((7, 3))
    adata3.obsm = {
        k: v[:, : v.shape[1] // 3]
        for k, v in adata3.obsm.items()
        if not isinstance(v, pd.DataFrame | AwkArray | XDataset)
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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


@mark_legacy_concatenate
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

    adata1.raw = adata1.copy()
    adata2.raw = adata2.copy()
    adata3.raw = adata3.copy()

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

    adata3.raw = adata4.copy()
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
        match=r"Cannot concatenate an AwkwardArray with other array types",
    ):
        concat([adata_a, adata_b], join=join_type)


def test_pairwise_concat(axis_name, array_type):
    axis, axis_name = merge._resolve_axis(axis_name)
    _, alt_axis_name = merge._resolve_axis(1 - axis)
    axis_sizes = [[100, 200, 50], [50, 50, 50]]
    if axis_name == "var":
        axis_sizes.reverse()
    Ms, Ns = axis_sizes
    axis_attr = f"{axis_name}p"
    alt_attr = f"{alt_axis_name}p"

    def gen_axis_array(m):
        return array_type(sparse.random(m, m, format="csr", density=0.1))

    adatas = {
        k: AnnData(
            X=sparse.csr_matrix((m, n)),
            obsp={"arr": gen_axis_array(m)},
            varp={"arr": gen_axis_array(n)},
        )
        for k, m, n in zip("abc", Ms, Ns, strict=True)
    }

    w_pairwise = concat(adatas, axis=axis, label="orig", pairwise=True)
    wo_pairwise = concat(adatas, axis=axis, label="orig", pairwise=False)

    # Check that argument controls whether elements are included
    assert getattr(wo_pairwise, axis_attr) == {}
    assert getattr(w_pairwise, axis_attr) != {}

    # Check values of included elements
    full_inds = np.arange(w_pairwise.shape[axis])
    obs_var: pd.DataFrame = getattr(w_pairwise, axis_name)
    groups = obs_var.groupby("orig", observed=True).indices
    for k, inds in groups.items():
        orig_arr = getattr(adatas[k], axis_attr)["arr"]
        full_arr = getattr(w_pairwise, axis_attr)["arr"]

        if isinstance(full_arr, DaskArray):
            full_arr = full_arr.compute()

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


def test_nan_merge(axis_name, join_type, array_type):
    axis, _ = merge._resolve_axis(axis_name)
    alt_axis, alt_axis_name = merge._resolve_axis(1 - axis)
    mapping_attr = f"{alt_axis_name}m"
    adata_shape = (20, 10)

    arr = array_type(
        sparse.random(adata_shape[alt_axis], 10, density=0.1, format="csr")
    )
    arr_nan = arr.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=sparse.SparseEfficiencyWarning)
        for _ in range(10):
            arr_nan[np.random.choice(arr.shape[0]), np.random.choice(arr.shape[1])] = (
                np.nan
            )

    _data = {"X": sparse.csr_matrix(adata_shape), mapping_attr: {"arr": arr_nan}}
    orig1 = AnnData(**_data)
    orig2 = AnnData(**_data)
    result = concat([orig1, orig2], axis=axis, join=join_type, merge="same")

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


def permute_nested_values(dicts: list[dict], gen_val: Callable[[int], Any]):
    """
    This function permutes the values of a nested mapping, for testing that out merge
    method work regardless of the values types.

    Assumes the initial dictionary had integers for values.
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
    return sparse.random(
        np.random.randint(1, 100), np.random.randint(1, 100), format="csr"
    )


def gen_something(n):
    options = [gen_df, gen_array, gen_list, gen_sparse]
    return np.random.choice(options)(n)


def gen_3d_numeric_array(n):
    return np.random.randn(n, n, n)


def gen_3d_recarray(_):
    # Ignoring n as it can get quite slow
    return gen_vstr_recarray(8, 3).reshape(2, 2, 2)


def gen_concat_params(unss, compat2result):
    value_generators = [
        lambda x: x,
        gen_df,
        gen_array,
        gen_list,
        gen_sparse,
        gen_something,
        gen_3d_numeric_array,
        gen_3d_recarray,
    ]
    for gen, (mode, result) in product(value_generators, compat2result.items()):
        yield pytest.param(unss, mode, result, gen)


@pytest.mark.parametrize(
    ("unss", "merge_strategy", "result", "value_gen"),
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
    result, *unss = permute_nested_values([result, *unss], value_gen)
    adatas = [uns_ad(uns) for uns in unss]
    with pytest.warns(FutureWarning, match=r"concatenate is deprecated"):
        merged = AnnData.concatenate(*adatas, uns_merge=merge_strategy).uns
    assert_equal(merged, result, elem_name="uns")


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_transposed_concat(
    array_type,
    axis_name,
    join_type,
    merge_strategy,
    obs_xdataset,
    var_xdataset,
    force_lazy,
):
    axis, axis_name = merge._resolve_axis(axis_name)
    alt_axis = 1 - axis
    lhs = gen_adata(
        (10, 10),
        X_type=array_type,
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )
    rhs = gen_adata(
        (10, 12),
        X_type=array_type,
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )

    a = concat(
        [lhs, rhs],
        axis=axis,
        join=join_type,
        merge=merge_strategy,
        force_lazy=force_lazy,
    )
    b = concat(
        [lhs.T, rhs.T],
        axis=alt_axis,
        join=join_type,
        merge=merge_strategy,
        force_lazy=force_lazy,
    ).T

    assert_equal(a, b)


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_batch_key(axis_name, obs_xdataset, var_xdataset, force_lazy):
    """Test that concat only adds a label if the key is provided"""

    get_annot = attrgetter(axis_name)

    lhs = gen_adata(
        (10, 10),
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )
    rhs = gen_adata(
        (10, 12),
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )

    # There is probably a prettier way to do this
    annot = get_annot(concat([lhs, rhs], axis=axis_name, force_lazy=force_lazy))
    assert (
        list(
            annot.columns.difference(
                get_annot(lhs).columns.union(get_annot(rhs).columns)
            )
        )
        == []
    )

    batch_annot = get_annot(
        concat([lhs, rhs], axis=axis_name, label="batch", force_lazy=force_lazy)
    )
    assert list(
        batch_annot.columns.difference(
            get_annot(lhs).columns.union(get_annot(rhs).columns)
        )
    ) == ["batch"]


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concat_categories_from_mapping(obs_xdataset, var_xdataset, force_lazy):
    mapping = {
        "a": gen_adata((10, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset),
        "b": gen_adata((10, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset),
    }
    keys = list(mapping.keys())
    adatas = list(mapping.values())

    mapping_call = partial(concat, mapping, force_lazy=force_lazy)
    iter_call = partial(concat, adatas, keys=keys, force_lazy=force_lazy)

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

    assert isinstance(result.obs["cat"].dtype, pd.CategoricalDtype), (
        f"Was {result.obs['cat'].dtype}"
    )
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


def test_concat_categorical_dtype_promotion():
    """https://github.com/scverse/anndata/issues/1170

    When concatenating categorical with other dtype, defer to pandas.
    """
    a = AnnData(
        np.ones((3, 3)),
        obs=pd.DataFrame(
            {"col": pd.Categorical(["a", "a", "b"])},
            index=[f"cell_{i:02d}" for i in range(3)],
        ),
    )
    b = AnnData(
        np.ones((3, 3)),
        obs=pd.DataFrame(
            {"col": ["c", "c", "c"]},
            index=[f"cell_{i:02d}" for i in range(3, 6)],
        ),
    )

    result = concat([a, b])
    expected = pd.concat([a.obs, b.obs])

    assert_equal(result.obs, expected)


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


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concat_names(axis_name, obs_xdataset, var_xdataset, force_lazy):
    get_annot = attrgetter(axis_name)

    lhs = gen_adata((10, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset)
    rhs = gen_adata((10, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset)

    assert not get_annot(
        concat([lhs, rhs], axis=axis_name, force_lazy=force_lazy)
    ).index.is_unique
    assert get_annot(
        concat([lhs, rhs], axis=axis_name, index_unique="-", force_lazy=force_lazy)
    ).index.is_unique


def axis_labels(adata: AnnData, axis: Literal[0, 1]) -> pd.Index:
    return (adata.obs_names, adata.var_names)[axis]


def expected_shape(
    a: AnnData, b: AnnData, axis: Literal[0, 1], join: Literal["inner", "outer"]
) -> tuple[int, int]:
    alt_axis = 1 - axis
    labels = partial(axis_labels, axis=alt_axis)
    shape = [None, None]

    shape[axis] = a.shape[axis] + b.shape[axis]
    if join == "inner":
        shape[alt_axis] = len(labels(a).intersection(labels(b)))
    elif join == "outer":
        shape[alt_axis] = len(labels(a).union(labels(b)))
    else:
        raise ValueError()

    return tuple(shape)


@pytest.mark.parametrize(
    "shape", [pytest.param((8, 0), id="no_var"), pytest.param((0, 10), id="no_obs")]
)
@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concat_size_0_axis(
    axis_name, join_type, merge_strategy, shape, obs_xdataset, var_xdataset, force_lazy
):
    """Regression test for https://github.com/scverse/anndata/issues/526"""
    axis, axis_name = merge._resolve_axis(axis_name)
    alt_axis = 1 - axis
    col_dtypes = (*DEFAULT_COL_TYPES, pd.StringDtype)
    a = gen_adata(
        (5, 7),
        obs_dtypes=col_dtypes,
        var_dtypes=col_dtypes,
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
    )
    b = gen_adata(
        shape,
        obs_dtypes=col_dtypes,
        var_dtypes=col_dtypes,
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
    )

    expected_size = expected_shape(a, b, axis=axis, join=join_type)

    result = concat(
        {"a": a, "b": b},
        axis=axis,
        join=join_type,
        merge=merge_strategy,
        pairwise=True,
        index_unique="-",
        force_lazy=force_lazy,
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
        for k, elem in result.layers.items():
            check_filled_like(elem[axis_idx], elem_name=f"layers/{k}")
            check_filled_like(elem[altaxis_idx], elem_name=f"layers/{k}")

        if shape[axis] > 0:
            b_result = result[axis_idx].copy()
            mapping_elem = f"{axis_name}m"
            setattr(b_result, f"{axis_name}_names", getattr(b, f"{axis_name}_names"))
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
@pytest.mark.parametrize("axis", ["obs", "var"])
@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concat_outer_aligned_mapping(
    elem, axis, obs_xdataset, var_xdataset, force_lazy
):
    a = gen_adata(
        (5, 5),
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )
    b = gen_adata(
        (3, 5),
        obs_xdataset=obs_xdataset,
        var_xdataset=var_xdataset,
        **GEN_ADATA_DASK_ARGS,
    )
    del getattr(b, f"{axis}m")[elem]

    concated = concat(
        {"a": a, "b": b}, join="outer", label="group", axis=axis, force_lazy=force_lazy
    )

    mask = getattr(concated, axis)["group"] == "b"
    result = getattr(
        concated[(mask, slice(None)) if axis == "obs" else (slice(None), mask)],
        f"{axis}m",
    )[elem]

    check_filled_like(result, elem_name=f"{axis}m/{elem}")


@mark_legacy_concatenate
def test_concatenate_size_0_axis():
    # https://github.com/scverse/anndata/issues/526

    a = gen_adata((5, 10))
    b = gen_adata((5, 0))

    # Mostly testing that this doesn't error
    assert a.concatenate([b]).shape == (10, 0)
    assert b.concatenate([a]).shape == (10, 0)


@pytest.mark.parametrize(
    ("obs_xdataset", "var_xdataset", "force_lazy"),
    [(False, False, False), (True, True, False), (True, True, True)],
)
def test_concat_null_X(obs_xdataset, var_xdataset, force_lazy):
    adatas_orig = {
        k: gen_adata((20, 10), obs_xdataset=obs_xdataset, var_xdataset=var_xdataset)
        for k in list("abc")
    }
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
@pytest.mark.parametrize("sparse_indexer_type", [np.int64, np.int32])
def test_concat_X_dtype(cpu_array_type, sparse_indexer_type):
    adatas_orig = {
        k: AnnData(cpu_array_type(np.ones((20, 10), dtype=np.int8)))
        for k in list("abc")
    }
    for adata in adatas_orig.values():
        adata.raw = AnnData(cpu_array_type(np.ones((20, 30), dtype=np.float64)))
        if sparse.issparse(adata.X):
            adata.X.indptr = adata.X.indptr.astype(sparse_indexer_type)
            adata.X.indices = adata.X.indices.astype(sparse_indexer_type)

    result = concat(adatas_orig, index_unique="-")

    assert result.X.dtype == np.int8
    assert result.raw.X.dtype == np.float64
    if sparse.issparse(result.X):
        # https://github.com/scipy/scipy/issues/20389 was merged in 1.15 but is still an issue with matrix
        if sparse_indexer_type == np.int64 and (
            (
                (issubclass(cpu_array_type, CSArray) or adata.X.format == "csc")
                and Version(scipy.__version__) < Version("1.15.0")
            )
            or issubclass(cpu_array_type, CSMatrix)
        ):
            pytest.xfail(
                "Data type int64 is not maintained for sparse matrices or csc array"
            )
        assert result.X.indptr.dtype == sparse_indexer_type, result.X
        assert result.X.indices.dtype == sparse_indexer_type


# Leaving out for now. See definition of these values for explanation
# def test_concatenate_uns_types():
#     from anndata._core.merge import UNS_STRATEGIES, UNS_STRATEGIES_TYPE
#     assert set(UNS_STRATEGIES.keys()) == set(UNS_STRATEGIES_TYPE.__args__)


# Tests how dask plays with other types on concatenation.
def test_concat_different_types_dask(merge_strategy, array_type):
    import dask.array as da
    from scipy import sparse

    import anndata as ad

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


def test_concat_missing_elem_dask_join(join_type):
    import dask.array as da

    import anndata as ad

    ad1 = ad.AnnData(X=np.ones((5, 10)))
    ad2 = ad.AnnData(X=np.zeros((5, 5)), layers={"a": da.ones((5, 5))})
    ad_in_memory_with_layers = ad2.to_memory()

    result1 = ad.concat([ad1, ad2], join=join_type)
    result2 = ad.concat([ad1, ad_in_memory_with_layers], join=join_type)
    assert_equal(result1, result2)


def test_impute_dask(axis_name):
    import dask.array as da

    from anndata._core.merge import _resolve_axis, missing_element

    axis, _ = _resolve_axis(axis_name)
    els = [da.ones((5, 5))]
    missing = missing_element(6, els, axis=axis, off_axis_size=17)
    assert isinstance(missing, DaskArray)
    in_memory = missing.compute()
    assert np.all(np.isnan(in_memory))
    assert in_memory.shape[axis] == 6
    assert in_memory.shape[axis - 1] == 17


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


@pytest.mark.gpu
def test_error_on_mixed_device():
    """https://github.com/scverse/anndata/issues/1083"""
    import cupy
    import cupyx.scipy.sparse as cupy_sparse

    cp_adata = AnnData(
        cupy.random.randn(10, 10),
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(10)]),
    )
    cp_sparse_adata = AnnData(
        cupy_sparse.random(10, 10, format="csr", density=0.2),
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(10, 20)]),
    )
    np_adata = AnnData(
        np.random.randn(10, 10),
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(20, 30)]),
    )
    sparse_adata = AnnData(
        sparse.random(10, 10, format="csr", density=0.2),
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(30, 40)]),
    )

    adatas = {
        "cupy": cp_adata,
        "cupy_sparse": cp_sparse_adata,
        "numpy": np_adata,
        "sparse": sparse_adata,
    }

    for p in map(dict, permutations(adatas.items())):
        print(list(p.keys()))
        with pytest.raises(
            NotImplementedError, match=r"Cannot concatenate a cupy array with other"
        ):
            concat(p)

    for p in permutations([cp_adata, cp_sparse_adata]):
        concat(p)


def test_concat_on_var_outer_join(array_type):
    # https://github.com/scverse/anndata/issues/1286
    a = AnnData(
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(10)]),
        var=pd.DataFrame(index=[f"gene_{i:02d}" for i in range(10)]),
        layers={
            "X": array_type(np.ones((10, 10))),
        },
    )
    b = AnnData(
        obs=pd.DataFrame(index=[f"cell_{i:02d}" for i in range(10)]),
        var=pd.DataFrame(index=[f"gene_{i:02d}" for i in range(10, 20)]),
    )

    # This shouldn't error
    # TODO: specify expected result while accounting for null value
    _ = concat([a, b], join="outer", axis=1)


def test_concat_dask_sparse_matches_memory(join_type, merge_strategy):
    import dask.array as da

    X = sparse.random(50, 20, density=0.5, format="csr")
    X_dask = da.from_array(X, chunks=(5, 20))
    var_names_1 = [f"gene_{i}" for i in range(20)]
    var_names_2 = [f"gene_{i}{'_foo' if (i % 2) else ''}" for i in range(20, 40)]

    ad1 = AnnData(X=X, var=pd.DataFrame(index=var_names_1))
    ad2 = AnnData(X=X, var=pd.DataFrame(index=var_names_2))

    ad1_dask = AnnData(X=X_dask, var=pd.DataFrame(index=var_names_1))
    ad2_dask = AnnData(X=X_dask, var=pd.DataFrame(index=var_names_2))

    res_in_memory = concat([ad1, ad2], join=join_type, merge=merge_strategy)
    res_dask = concat([ad1_dask, ad2_dask], join=join_type, merge=merge_strategy)

    assert_equal(res_in_memory, res_dask)


def test_1d_concat():
    adata = AnnData(np.ones((5, 20)), obsm={"1d-array": np.ones(5)})
    concated = concat([adata, adata])
    assert concated.obsm["1d-array"].shape == (10, 1)
