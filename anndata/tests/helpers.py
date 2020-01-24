from functools import singledispatch, wraps
from string import ascii_letters
from typing import Tuple
from collections.abc import Mapping

import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import pytest
from scipy import sparse

from anndata import AnnData
from anndata._core.views import ArrayView
from anndata._core.sparse_dataset import SparseDataset
from anndata._core.aligned_mapping import AlignedMapping


@singledispatch
def asarray(x):
    """Convert x to a numpy array"""
    return np.asarray(x)


@asarray.register(sparse.spmatrix)
def asarray_sparse(x):
    return x.toarray()


@asarray.register(SparseDataset)
def asarray_sparse_dataset(x):
    return asarray(x.value)


@asarray.register(h5py.Dataset)
def asarray_h5py_dataset(x):
    return x[...]


def gen_vstr_recarray(m, n, dtype=None):
    size = m * n
    lengths = np.random.randint(3, 5, size)
    letters = np.array(list(ascii_letters))
    gen_word = lambda l: "".join(np.random.choice(letters, l))
    arr = np.array([gen_word(l) for l in lengths]).reshape(m, n)
    return pd.DataFrame(arr, columns=[gen_word(5) for i in range(n)]).to_records(
        index=False, column_dtypes=dtype
    )


def gen_typed_df(n, index=None):
    # TODO: Think about allowing index to be passed for n
    letters = np.fromiter(iter(ascii_letters), "U1")
    if n > len(letters):
        letters = letters[: n // 2]  # Make sure categories are repeated
    return pd.DataFrame(
        dict(
            cat=pd.Categorical(np.random.choice(letters, n)),
            cat_ordered=pd.Categorical(np.random.choice(letters, n), ordered=True),
            int64=np.random.randint(-50, 50, n),
            float64=np.random.random(n),
            uint8=np.random.randint(255, size=n, dtype="uint8"),
        ),
        index=index,
    )


def gen_typed_df_t2_size(m, n, index=None, columns=None) -> pd.DataFrame:
    s = 0
    df = pd.DataFrame()
    new_vals = gen_typed_df(m)
    while s < (n / new_vals.shape[1]):
        new_vals = gen_typed_df(m, index=index)
        new_vals.columns = new_vals.columns + "_" + str(s)
        df[new_vals.columns] = new_vals
        s += 1
    df = df.iloc[:m, :n].copy()
    if columns is not None:
        df.columns = columns
    return df


# TODO: Use hypothesis for this?
def gen_adata(
    shape: Tuple[int, int],
    X_type=sparse.csr_matrix,
    X_dtype=np.float32,
    # obs_dtypes,
    # var_dtypes,
    obsm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame,),
    varm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame,),
    layers_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame,),
) -> AnnData:
    """\
    Helper function to generate a random AnnData for testing purposes.

    Note: For `obsm_types`, `varm_types`, and `layers_types` these currently
    just filter already created objects.
    In future, these should choose which objects are created.

    Params
    ------
    shape
        What shape you want the anndata to be.
    X_type
        What kind of container should `X` be? This will be called on a randomly
        generated 2d array.
    X_dtype
        What should the dtype of the `.X` container be?
    obsm_types
        What kinds of containers should be in `.obsm`?
    varm_types
        What kinds of containers should be in `.varm`?
    layers_types
        What kinds of containers should be in `.layers`?
    """
    M, N = shape
    obs_names = pd.Index(f"cell{i}" for i in range(shape[0]))
    var_names = pd.Index(f"gene{i}" for i in range(shape[1]))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    # For #147
    obs.rename(columns=dict(cat="obs_cat"), inplace=True)
    var.rename(columns=dict(cat="var_cat"), inplace=True)

    obsm = dict(
        array=np.random.random((M, 50)),
        sparse=sparse.random(M, 100, format="csr"),
        df=gen_typed_df(M, obs_names),
    )
    obsm = {k: v for k, v in obsm.items() if type(v) in obsm_types}
    varm = dict(
        array=np.random.random((N, 50)),
        sparse=sparse.random(N, 100, format="csr"),
        df=gen_typed_df(N, var_names),
    )
    varm = {k: v for k, v in varm.items() if type(v) in varm_types}
    layers = dict(
        array=np.random.random((M, N)), sparse=sparse.random(M, N, format="csr")
    )
    layers = {k: v for k, v in layers.items() if type(v) in layers_types}
    obsp = dict(
        array=np.random.random((M, M)), sparse=sparse.random(M, M, format="csr")
    )
    varp = dict(
        array=np.random.random((N, N)), sparse=sparse.random(N, N, format="csr")
    )
    uns = dict(
        O_recarray=gen_vstr_recarray(N, 5),
        # U_recarray=gen_vstr_recarray(N, 5, "U4")
    )
    adata = AnnData(
        X=X_type(np.random.binomial(100, 0.005, (M, N)).astype(X_dtype)),
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        layers=layers,
        obsp=obsp,
        varp=varp,
        dtype=X_dtype,
        uns=uns,
    )
    return adata


def array_bool_subset(index, min_size=2):
    b = np.zeros(len(index), dtype=bool)
    selected = np.random.choice(
        range(len(index)),
        size=np.random.randint(min_size, len(index), ()),
        replace=False,
    )
    b[selected] = True
    return b


def array_subset(index, min_size=2):
    if len(index) < min_size:
        raise ValueError(
            f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        )
    return np.random.choice(
        index, size=np.random.randint(min_size, len(index), ()), replace=False
    )


def array_int_subset(index, min_size=2):
    if len(index) < min_size:
        raise ValueError(
            f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        )
    return np.random.choice(
        np.arange(len(index)),
        size=np.random.randint(min_size, len(index), ()),
        replace=False,
    )


def slice_subset(index, min_size=2):
    while True:
        points = np.random.choice(np.arange(len(index) + 1), size=2, replace=False)
        s = slice(*sorted(points))
        if len(range(*s.indices(len(index)))) >= min_size:
            break
    return s


def single_subset(index):
    return index[np.random.randint(0, len(index), size=())]


@pytest.fixture(
    params=[
        array_subset,
        slice_subset,
        single_subset,
        array_int_subset,
        array_bool_subset,
    ]
)
def subset_func(request):
    return request.param


###################
# Checking equality
###################


def format_msg(elem_name):
    if elem_name is not None:
        return f"Error raised from element {elem_name!r}."
    else:
        return ""


# TODO: it would be better to modify the other exception
def report_name(func):
    """Report name of element being tested if test fails."""

    @wraps(func)
    def func_wrapper(*args, _elem_name=None, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if _elem_name is not None and not hasattr(e, "_name_attached"):
                msg = format_msg(_elem_name)
                args = list(e.args)
                if len(args) == 0:
                    args = [msg]
                else:
                    args[0] = f"{args[0]}\n\n{msg}"
                e.args = tuple(args)
                e._name_attached = True
            raise e

    return func_wrapper


@report_name
def _assert_equal(a, b):
    """Allows reporting elem name for simple assertion."""
    assert a == b


@singledispatch
def assert_equal(a, b, exact=False, elem_name=None):
    _assert_equal(a, b, _elem_name=elem_name)


@assert_equal.register(np.ndarray)
def assert_equal_ndarray(a, b, exact=False, elem_name=None):
    b = asarray(b)
    if not exact and is_numeric_dtype(a) and is_numeric_dtype(b):
        assert a.shape == b.shape, format_msg(elem_name)
        assert np.allclose(a, b, equal_nan=True), format_msg(elem_name)
    elif (  # Structured dtype
        not exact
        and hasattr(a, "dtype")
        and hasattr(b, "dtype")
        and len(a.dtype) > 1
        and len(b.dtype) > 0
    ):
        assert_equal(pd.DataFrame(a), pd.DataFrame(b), exact, elem_name)
    else:
        assert np.all(a == b), format_msg(elem_name)


@assert_equal.register(ArrayView)
def assert_equal_arrayview(a, b, exact=False, elem_name=None):
    assert_equal(asarray(a), asarray(b), exact=exact, elem_name=elem_name)


@assert_equal.register(SparseDataset)
@assert_equal.register(sparse.spmatrix)
def assert_equal_sparse(a, b, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact, elem_name=elem_name)


@assert_equal.register(h5py.Dataset)
def assert_equal_h5py_dataset(a, b, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact, elem_name=elem_name)


@assert_equal.register(pd.DataFrame)
def are_equal_dataframe(a, b, exact=False, elem_name=None):
    if not isinstance(b, pd.DataFrame):
        assert_equal(b, a, exact, elem_name)  # , a.values maybe?

    report_name(pd.testing.assert_frame_equal)(
        a,
        b,
        check_index_type=exact,
        check_exact=exact,
        _elem_name=elem_name,
        check_frame_type=False,
    )


@assert_equal.register(Mapping)
def assert_equal_mapping(a, b, exact=False, elem_name=None):
    assert set(a.keys()) == set(b.keys()), format_msg(elem_name)
    for k in a.keys():
        if elem_name is None:
            elem_name = ""
        assert_equal(a[k], b[k], exact, f"{elem_name}/{k}")


@assert_equal.register(AlignedMapping)
def assert_equal_aligned_mapping(a, b, exact=False, elem_name=None):
    a_indices = (a.parent.obs_names, a.parent.var_names)
    b_indices = (b.parent.obs_names, b.parent.var_names)
    for axis_idx in a.axes:
        assert_equal(
            a_indices[axis_idx], b_indices[axis_idx], exact=exact, elem_name=axis_idx,
        )
    assert a.attrname == b.attrname, format_msg(elem_name)
    assert_equal_mapping(a, b, exact=exact, elem_name=elem_name)


@assert_equal.register(pd.Index)
def assert_equal_index(a, b, exact=False, elem_name=None):
    if not exact:
        report_name(pd.testing.assert_index_equal)(
            a, b, check_names=False, check_categorical=False, _elem_name=elem_name,
        )
    else:
        report_name(pd.testing.assert_index_equal)(a, b, _elem_name=elem_name)


@assert_equal.register(AnnData)
def assert_adata_equal(a: AnnData, b: AnnData, exact: bool = False):
    """\
    Check whether two AnnData objects are equivalent,
    raising an AssertionError if they aren’t.

    Params
    ------
    a
    b
    exact
        Whether comparisons should be exact or not. This has a somewhat flexible
        meaning and should probably get refined in the future.
    """
    # There may be issues comparing views, since np.allclose
    # can modify ArrayViews if they contain `nan`s
    assert_equal(a.obs_names, b.obs_names, exact, elem_name="obs_names")
    assert_equal(a.var_names, b.var_names, exact, elem_name="var_names")
    if not exact:
        # Reorder all elements if neccesary
        idx = [slice(None), slice(None)]
        # Since it’s a pain to compare a list of pandas objects
        change_flag = False
        if not np.all(a.obs_names == b.obs_names):
            idx[0] = a.obs_names
            change_flag = True
        if not np.all(a.var_names == b.var_names):
            idx[1] = a.var_names
            change_flag = True
        if change_flag:
            b = b[tuple(idx)].copy()
    assert_equal(a.obs, b.obs, exact, elem_name="obs")
    assert_equal(a.var, b.var, exact, elem_name="var")
    assert_equal(a.X, b.X, exact, elem_name="X")
    for mapping_attr in ["obsm", "varm", "layers", "uns", "obsp", "varp"]:
        assert_equal(
            getattr(a, mapping_attr),
            getattr(b, mapping_attr),
            exact,
            elem_name=mapping_attr,
        )
    if a.raw is not None:
        assert_equal(a.raw.X, b.raw.X, exact, elem_name="raw/X")
        assert_equal(a.raw.var, b.raw.var, exact, elem_name="raw/var")
        assert_equal(a.raw.varm, b.raw.varm, exact, elem_name="raw/varm")
