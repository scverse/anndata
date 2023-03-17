from collections.abc import Hashable
from copy import deepcopy
from itertools import chain, product
from functools import partial, singledispatch
from typing import Any, List, Callable
import typing
import warnings

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
import pytest
from scipy import sparse
from boltons.iterutils import research, remap, default_exit


from anndata import AnnData, Raw, concat_on_disk, concat, read_zarr
from anndata._core.index import _subset
from anndata._core import merge
from anndata.tests import helpers
from anndata.tests.helpers import (
    assert_equal,
    as_dense_dask_array,
    gen_adata,
    GEN_ADATA_DASK_ARGS,
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
        assert_equal(x, filled_like(x, fill_value=fill_value),
                     elem_name=elem_name)


def make_idx_tuple(idx, axis):
    tup = [slice(None), slice(None)]
    tup[axis] = idx
    return tuple(tup)


# Will call func(sparse_matrix) so these types should be sparse compatible
# See array_type if only dense arrays are expected as input.
@pytest.fixture(
    params=[asarray, sparse.csr_matrix,
            sparse.csc_matrix, as_dense_dask_array],
    ids=["np_array", "scipy_csr", "scipy_csc", "dask_array"],
)
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
        if is_categorical_dtype(dtype) and dtype.ordered:
            result.obs[k] = result.obs[k].astype(dtype)

    return orig, result


def adatas_to_zarr_paths(adatas, tmp_path):
    """
    Gets list of adatas, writes them and returns their paths as zarr
    """
    paths = None
    if isinstance(adatas, typing.Mapping):
        paths = {}
        for k, v in adatas.items():
            p = tmp_path/f"{k}.zarr"
            v.write_zarr(p)
            paths[k] = p
    else:
        paths = []
        for i, a in enumerate(adatas):
            p = tmp_path/f"{i}.zarr"
            a.write_zarr(p)
            paths += [p]
    return paths


def test_concat_interface_errors_disk(tmp_path):
    adatas = [gen_adata((5, 10)), gen_adata((5, 10))]
    paths = adatas_to_zarr_paths(adatas, tmp_path)

    with pytest.raises(ValueError):
        concat_on_disk(paths, out_path=tmp_path/'test', axis=3)
    with pytest.raises(ValueError):
        concat_on_disk(paths, out_path=tmp_path/'test', join="not implemented")
    with pytest.raises(ValueError):
        concat_on_disk([], tmp_path/"test")


def assert_eq_concat_on_disk(adatas, tmp_path, *args, **kwargs):

    # create one from the concat function
    res1 = concat(adatas, *args, **kwargs)
    # create one from the on disk concat function
    paths = adatas_to_zarr_paths(adatas, tmp_path)
    out_name = tmp_path/"out.zarr"
    res_path = concat_on_disk(paths, out_name, *args, **kwargs)
    res2 = read_zarr(res_path)
    assert_equal(res1, res2)


def test_concatenate_dense(join_type, tmp_path):
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

    assert_eq_concat_on_disk([adata1, adata2, adata3],
                             tmp_path, join=join_type)


def test_concatenate_layers(array_type, join_type, tmp_path):
    adatas = []
    for _ in range(5):
        a = array_type(sparse.random(100, 200, format="csr"))
        adatas.append(AnnData(X=a, layers={"a": a}))

    assert_eq_concat_on_disk(adatas, tmp_path, join=join_type)


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


def test_concatenate_obsm(obsm_adatas, tmp_path, join_type):
    assert_eq_concat_on_disk(obsm_adatas, tmp_path, join=join_type)


def test_concatenate_layers_misaligned(array_type, join_type, tmp_path):
    adatas = []
    for _ in range(5):
        a = array_type(sparse.random(100, 200, format="csr"))
        adata = AnnData(X=a, layers={"a": a})
        adatas.append(
            adata[:, np.random.choice(
                adata.var_names, 150, replace=False)].copy()
        )

    assert_eq_concat_on_disk(adatas, tmp_path, join=join_type)


def test_concatenate_layers_outer(array_type, fill_val, tmp_path):
    # Testing that issue #368 is fixed
    a = AnnData(
        X=np.ones((10, 20)),
        layers={"a": array_type(sparse.random(10, 20, format="csr"))},
    )
    b = AnnData(X=np.ones((10, 20)))

    assert_eq_concat_on_disk(
        [a, b], tmp_path, join="outer", fill_value=fill_val)


def test_concatenate_dense_duplicates(tmp_path):
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

    assert_eq_concat_on_disk([adata1, adata2, adata3], tmp_path)


def test_concatenate_sparse(tmp_path):
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

    assert_eq_concat_on_disk([adata1, adata2, adata3], tmp_path)


def test_concatenate_mixed(tmp_path):
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

    assert_eq_concat_on_disk([adata1, adata2, adata3, adata4], tmp_path)


def test_concatenate_awkward(join_type, tmp_path):
    import awkward as ak

    a = ak.Array([[{"a": 1, "b": "foo"}], [
                 {"a": 2, "b": "bar"}, {"a": 3, "b": "baz"}]])
    b = ak.Array(
        [
            [{"a": 4}, {"a": 5}],
            [{"a": 6}],
            [{"a": 7}],
        ]
    )

    adata1 = AnnData(np.zeros((2, 0), dtype=float), obsm={"awk": a})
    adata2 = AnnData(np.zeros((3, 0), dtype=float), obsm={"awk": b})

    assert_eq_concat_on_disk([adata1, adata2], tmp_path, join=join_type)


@pytest.mark.parametrize(
    "other",
    [
        pd.DataFrame(
            {"a": [4, 5, 6], "b": ["foo", "bar", "baz"]}, index=list("cde")),
        np.ones((3, 2)),
        sparse.random(3, 100, format="csr"),
    ],
)
def test_awkward_does_not_mix(join_type, other, tmp_path):
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
        adatas = [adata_a, adata_b]
        paths = adatas_to_zarr_paths(adatas, tmp_path/"awk")
        concat_on_disk(paths, tmp_path/"out", join=join_type)


def test_pairwise_concat(axis, array_type, tmp_path):
    dim_sizes = [[100, 200, 50], [50, 50, 50]]
    if axis:
        dim_sizes.reverse()
    Ms, Ns = dim_sizes

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
    assert_eq_concat_on_disk(
        adatas, tmp_path, axis=axis, label="orig", pairwise=True)
    assert_eq_concat_on_disk(adatas, tmp_path, axis=axis,
                             label="orig", pairwise=False)


def test_nan_merge(axis, join_type, array_type, tmp_path):
    # concat_dim = ("obs", "var")[axis]
    alt_dim = ("var", "obs")[axis]
    mapping_attr = f"{alt_dim}m"
    adata_shape = (20, 10)

    arr = array_type(
        sparse.random(adata_shape[1 - axis], 10, density=0.1, format="csr")
    )
    arr_nan = arr.copy()
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore", category=sparse.SparseEfficiencyWarning)
        for _ in range(10):
            arr_nan[
                np.random.choice(arr.shape[0]), np.random.choice(arr.shape[1])
            ] = np.nan

    _data = {"X": sparse.csr_matrix(adata_shape), mapping_attr: {
        "arr": arr_nan}}
    orig1 = AnnData(**_data)
    orig2 = AnnData(**_data)
    orig_nonan = AnnData(
        **{"X": sparse.csr_matrix(adata_shape), mapping_attr: {"arr": arr}}
    )
    assert_eq_concat_on_disk(
        [orig1, orig2], tmp_path, axis=axis, merge="same")
    assert_eq_concat_on_disk(
        [orig1, orig_nonan], tmp_path, axis=axis, merge="same")


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
def test_concatenate_uns(unss, merge_strategy, result, value_gen, tmp_path):
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
    assert_eq_concat_on_disk(adatas, tmp_path, uns_merge=merge_strategy)


def test_transposed_concat(array_type, axis, join_type, merge_strategy, tmp_path):
    lhs = gen_adata((10, 10), X_type=array_type, **GEN_ADATA_DASK_ARGS)
    rhs = gen_adata((10, 12), X_type=array_type, **GEN_ADATA_DASK_ARGS)

    assert_eq_concat_on_disk([lhs, rhs],
                             tmp_path,
                             axis=axis,
                             join=join_type, merge=merge_strategy)
    assert_eq_concat_on_disk(
        [lhs.T, rhs.T],
        tmp_path,
        axis=abs(axis - 1), join=join_type, merge=merge_strategy
    )


def test_batch_key(axis, tmp_path):
    """Test that concat only adds a label if the key is provided"""

    lhs = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    rhs = gen_adata((10, 12), **GEN_ADATA_DASK_ARGS)

    assert_eq_concat_on_disk([lhs, rhs], tmp_path, axis=axis)
    assert_eq_concat_on_disk([lhs, rhs], tmp_path, axis=axis, label="batch")


def test_concat_categories_from_mapping(tmp_path):
    mapping = {
        "a": gen_adata((10, 10)),
        "b": gen_adata((10, 10)),
    }
    keys = list(mapping.keys())
    adatas = list(mapping.values())

    mapping_call = partial(assert_eq_concat_on_disk, mapping, tmp_path)
    iter_call = partial(assert_eq_concat_on_disk, adatas, tmp_path, keys=keys)

    mapping_call()
    iter_call()
    mapping_call(label="batch")
    iter_call(label="batch"),
    mapping_call(index_unique="-")
    iter_call(index_unique="-")
    mapping_call(label="group", index_unique="+")
    iter_call(label="group", index_unique="+")


def test_concat_categories_maintain_dtype(tmp_path):
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

    assert_eq_concat_on_disk({"a": a, "b": b, "c": c}, tmp_path, join="outer")


def test_concat_ordered_categoricals_retained(tmp_path):
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

    assert_eq_concat_on_disk([a, b], tmp_path)


def test_concat_names(axis, tmp_path):

    lhs = gen_adata((10, 10))
    rhs = gen_adata((10, 10))

    assert_eq_concat_on_disk([lhs, rhs], tmp_path, axis=axis)
    assert_eq_concat_on_disk([lhs, rhs], tmp_path,
                             axis=axis, index_unique="-")


@pytest.mark.parametrize(
    "shape", [pytest.param((8, 0), id="no_var"),
              pytest.param((0, 10), id="no_obs")]
)
def test_concat_size_0_dim(axis, join_type, merge_strategy, shape, tmp_path):
    # https://github.com/scverse/anndata/issues/526
    a = gen_adata((5, 7))
    b = gen_adata(shape)

    assert_eq_concat_on_disk(
        {"a": a, "b": b},
        tmp_path,
        axis=axis,
        join=join_type,
        merge=merge_strategy,
        pairwise=True,
        index_unique="-",
    )


@pytest.mark.parametrize("elem", ["sparse", "array", "df", "da"])
def test_concat_outer_aligned_mapping(elem, tmp_path):
    a = gen_adata((5, 5), **GEN_ADATA_DASK_ARGS)
    b = gen_adata((3, 5), **GEN_ADATA_DASK_ARGS)
    del b.obsm[elem]

    assert_eq_concat_on_disk(
        {"a": a, "b": b}, tmp_path, join="outer", label="group")


def test_concatenate_size_0_dim(tmp_path):
    # https://github.com/scverse/anndata/issues/526

    a = gen_adata((5, 10))
    b = gen_adata((5, 0))

    assert_eq_concat_on_disk([a, b], tmp_path)
    assert_eq_concat_on_disk([b, a], tmp_path)


def test_concat_null_X(tmp_path):
    adatas_orig = {k: gen_adata((20, 10)) for k in list("abc")}
    adatas_no_X = {}
    for k, v in adatas_orig.items():
        v = v.copy()
        del v.X
        adatas_no_X[k] = v

    assert_eq_concat_on_disk(adatas_orig, tmp_path, index_unique="-")
    assert_eq_concat_on_disk(adatas_no_X, tmp_path, index_unique="-")


# https://github.com/scverse/ehrapy/issues/151#issuecomment-1016753744
def test_concat_X_dtype(tmp_path):
    adatas_orig = {k: AnnData(np.ones((20, 10), dtype=np.int8))
                   for k in list("abc")}
    for adata in adatas_orig.values():
        adata.raw = AnnData(np.ones((20, 30), dtype=np.float64))

    assert_eq_concat_on_disk(adatas_orig, tmp_path, index_unique="-")


# Tests how dask plays with other types on concatenation.
def test_concat_different_types_dask(merge_strategy, array_type, tmp_path):
    from scipy import sparse
    import anndata as ad
    import dask.array as da

    varm_array = sparse.random(5, 20, density=0.5, format="csr")

    ad1 = ad.AnnData(X=np.ones((5, 5)), varm={"a": varm_array})
    ad1_other = ad.AnnData(X=np.ones((5, 5)), varm={
                           "a": array_type(varm_array)})
    ad2 = ad.AnnData(X=np.zeros((5, 5)), varm={"a": da.ones(5, 20)})

    assert_eq_concat_on_disk([ad1, ad2], tmp_path, merge=merge_strategy)
    assert_eq_concat_on_disk([ad1_other, tmp_path, ad2], merge=merge_strategy)
    assert_eq_concat_on_disk([ad2, ad1], tmp_path, merge=merge_strategy)
    assert_eq_concat_on_disk([ad2, ad1_other], tmp_path, merge=merge_strategy)


def test_outer_concat_with_missing_value_for_df(tmp_path):
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

    assert_eq_concat_on_disk([a, b], tmp_path, join="outer")
