from __future__ import annotations

from contextlib import contextmanager
from functools import singledispatch, wraps, partial
import re
from string import ascii_letters
from typing import Tuple, Optional, Type
from collections.abc import Mapping, Collection
import warnings

import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import pytest
from scipy import sparse
import random

from anndata import AnnData, Raw
from anndata._core.views import ArrayView
from anndata._core.sparse_dataset import BaseCompressedSparseDataset
from anndata._core.aligned_mapping import AlignedMapping
from anndata.utils import asarray
from anndata.compat import (
    AwkArray,
    DaskArray,
    CupySparseMatrix,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
)

# Give this to gen_adata when dask array support is expected.
GEN_ADATA_DASK_ARGS = dict(
    obsm_types=(
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
        DaskArray,
    ),
    varm_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame, DaskArray),
    layers_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame, DaskArray),
)


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
        {
            "cat": pd.Categorical(np.random.choice(letters, n)),
            "cat_ordered": pd.Categorical(np.random.choice(letters, n), ordered=True),
            "int64": np.random.randint(-50, 50, n),
            "float64": np.random.random(n),
            "uint8": np.random.randint(255, size=n, dtype="uint8"),
            "bool": np.random.randint(0, 2, size=n, dtype=bool),
            "nullable-bool": pd.arrays.BooleanArray(
                np.random.randint(0, 2, size=n, dtype=bool),
                mask=np.random.randint(0, 2, size=n, dtype=bool),
            ),
            "nullable-int": pd.arrays.IntegerArray(
                np.random.randint(0, 1000, size=n, dtype=np.int32),
                mask=np.random.randint(0, 2, size=n, dtype=bool),
            ),
        },
        index=index,
    )


def _gen_awkward_inner(shape, rng, dtype):
    # the maximum length a ragged dimension can take
    MAX_RAGGED_DIM_LEN = 20
    if not len(shape):
        # abort condition -> no dimension left, return an actual value instead
        return dtype(rng.randrange(1000))
    else:
        curr_dim_len = shape[0]
        lil = []
        if curr_dim_len is None:
            # ragged dimension, set random length
            curr_dim_len = rng.randrange(MAX_RAGGED_DIM_LEN)

        for _ in range(curr_dim_len):
            lil.append(_gen_awkward_inner(shape[1:], rng, dtype))

        return lil


def gen_awkward(shape, dtype=np.int32):
    """Function to generate an awkward array with random values.

    Awkward array dimensions can either be fixed-length ("regular") or variable length ("ragged")
    (the first dimension is always fixed-length).


    Parameters
    ----------
    shape
        shape of the array to be generated. Any dimension specified as `None` will be simulated as ragged.
    """
    import awkward as ak

    if shape[0] is None:
        raise ValueError("The first dimension must be fixed-length.")

    rng = random.Random(123)
    shape = np.array(shape)

    if np.any(shape == 0):
        # use empty numpy array for fixed dimensions, then add empty singletons for ragged dimensions
        var_dims = [i for i, s in enumerate(shape) if s is None]
        shape = [s for s in shape if s is not None]
        arr = ak.Array(np.empty(shape, dtype=dtype))
        for d in var_dims:
            arr = ak.singletons(arr, axis=d - 1)
        return arr
    else:
        lil = _gen_awkward_inner(shape, rng, dtype)
        arr = ak.values_astype(AwkArray(lil), dtype)

    # make fixed-length dimensions regular
    for i, d in enumerate(shape):
        if d is not None:
            arr = ak.to_regular(arr, i)

    return arr


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
    obsm_types: "Collection[Type]" = (
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
        AwkArray,
    ),
    varm_types: "Collection[Type]" = (
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
        AwkArray,
    ),
    layers_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
    sparse_fmt: str = "csr",
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
    sparse_fmt
        What sparse format should be used for sparse matrices?
        (csr, csc)
    """
    import dask.array as da

    M, N = shape
    obs_names = pd.Index(f"cell{i}" for i in range(shape[0]))
    var_names = pd.Index(f"gene{i}" for i in range(shape[1]))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    # For #147
    obs.rename(columns=dict(cat="obs_cat"), inplace=True)
    var.rename(columns=dict(cat="var_cat"), inplace=True)

    if X_type is None:
        X = None
    else:
        X = X_type(np.random.binomial(100, 0.005, (M, N)).astype(X_dtype))
    obsm = dict(
        array=np.random.random((M, 50)),
        sparse=sparse.random(M, 100, format=sparse_fmt),
        df=gen_typed_df(M, obs_names),
        awk_2d_ragged=gen_awkward((M, None)),
        da=da.random.random((M, 50)),
    )
    obsm = {k: v for k, v in obsm.items() if type(v) in obsm_types}
    varm = dict(
        array=np.random.random((N, 50)),
        sparse=sparse.random(N, 100, format=sparse_fmt),
        df=gen_typed_df(N, var_names),
        awk_2d_ragged=gen_awkward((N, None)),
        da=da.random.random((N, 50)),
    )
    varm = {k: v for k, v in varm.items() if type(v) in varm_types}
    layers = dict(
        array=np.random.random((M, N)),
        sparse=sparse.random(M, N, format=sparse_fmt),
        da=da.random.random((M, N)),
    )
    layers = {k: v for k, v in layers.items() if type(v) in layers_types}
    obsp = dict(
        array=np.random.random((M, M)), sparse=sparse.random(M, M, format=sparse_fmt)
    )
    varp = dict(
        array=np.random.random((N, N)), sparse=sparse.random(N, N, format=sparse_fmt)
    )
    uns = dict(
        O_recarray=gen_vstr_recarray(N, 5),
        nested=dict(
            scalar_str="str",
            scalar_int=42,
            scalar_float=3.0,
            nested_further=dict(array=np.arange(5)),
        ),
        awkward_regular=gen_awkward((10, 5)),
        awkward_ragged=gen_awkward((12, None, None)),
        # U_recarray=gen_vstr_recarray(N, 5, "U4")
    )
    adata = AnnData(
        X=X,
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        layers=layers,
        obsp=obsp,
        varp=varp,
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


def matrix_bool_subset(index, min_size=2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PendingDeprecationWarning)
        indexer = np.matrix(
            array_bool_subset(index, min_size=min_size).reshape(len(index), 1)
        )
    return indexer


def spmatrix_bool_subset(index, min_size=2):
    return sparse.csr_matrix(
        array_bool_subset(index, min_size=min_size).reshape(len(index), 1)
    )


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
    return index[np.random.randint(0, len(index))]


@pytest.fixture(
    params=[
        array_subset,
        slice_subset,
        single_subset,
        array_int_subset,
        array_bool_subset,
        matrix_bool_subset,
        spmatrix_bool_subset,
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


@assert_equal.register(CupyArray)
def assert_equal_cupy(a, b, exact=False, elem_name=None):
    assert_equal(b, a.get(), exact, elem_name)


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


@assert_equal.register(BaseCompressedSparseDataset)
@assert_equal.register(sparse.spmatrix)
def assert_equal_sparse(a, b, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact, elem_name=elem_name)


@assert_equal.register(CupySparseMatrix)
def assert_equal_cupy_sparse(a, b, exact=False, elem_name=None):
    a = a.toarray()
    assert_equal(b, a, exact, elem_name=elem_name)


@assert_equal.register(h5py.Dataset)
def assert_equal_h5py_dataset(a, b, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact, elem_name=elem_name)


@assert_equal.register(DaskArray)
def assert_equal_dask_array(a, b, exact=False, elem_name=None):
    from dask.array.utils import assert_eq

    if exact:
        assert_eq(a, b, check_dtype=True, check_type=True, check_graph=False)
    else:
        # TODO: Why does it fail when check_graph=True
        assert_eq(a, b, check_dtype=False, check_type=False, check_graph=False)


@assert_equal.register(pd.DataFrame)
def are_equal_dataframe(a, b, exact=False, elem_name=None):
    if not isinstance(b, pd.DataFrame):
        assert_equal(b, a, exact, elem_name)  # , a.values maybe?

    report_name(pd.testing.assert_frame_equal)(
        a,
        b,
        check_exact=exact,
        check_column_type=exact,
        check_index_type=exact,
        _elem_name=elem_name,
        check_frame_type=False,
    )


@assert_equal.register(AwkArray)
def assert_equal_awkarray(a, b, exact=False, elem_name=None):
    import awkward as ak

    if exact:
        assert a.type == b.type, f"{a.type} != {b.type}, {format_msg(elem_name)}"
    assert ak.to_list(a) == ak.to_list(b), format_msg(elem_name)


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
            a_indices[axis_idx], b_indices[axis_idx], exact=exact, elem_name=axis_idx
        )
    assert a.attrname == b.attrname, format_msg(elem_name)
    assert_equal_mapping(a, b, exact=exact, elem_name=elem_name)


@assert_equal.register(pd.Index)
def assert_equal_index(a, b, exact=False, elem_name=None):
    if not exact:
        report_name(pd.testing.assert_index_equal)(
            a, b, check_names=False, check_categorical=False, _elem_name=elem_name
        )
    else:
        report_name(pd.testing.assert_index_equal)(a, b, _elem_name=elem_name)


@assert_equal.register(pd.api.extensions.ExtensionArray)
def assert_equal_extension_array(a, b, exact=False, elem_name=None):
    report_name(pd.testing.assert_extension_array_equal)(
        a,
        b,
        check_dtype=exact,
        check_exact=exact,
        _elem_name=elem_name,
    )


@assert_equal.register(Raw)
def assert_equal_raw(a, b, exact=False, elem_name=None):
    def assert_is_not_none(x):  # can't put an assert in a lambda
        assert x is not None

    report_name(assert_is_not_none)(b, _elem_name=elem_name)
    for attr in ["X", "var", "varm", "obs_names"]:
        assert_equal(
            getattr(a, attr),
            getattr(b, attr),
            exact=exact,
            elem_name=f"{elem_name}/{attr}",
        )


@assert_equal.register(AnnData)
def assert_adata_equal(
    a: AnnData, b: AnnData, exact: bool = False, elem_name: Optional[str] = None
):
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

    def fmt_name(x):
        if elem_name is None:
            return x
        else:
            return f"{elem_name}/{x}"

    # There may be issues comparing views, since np.allclose
    # can modify ArrayViews if they contain `nan`s
    assert_equal(a.obs_names, b.obs_names, exact, elem_name=fmt_name("obs_names"))
    assert_equal(a.var_names, b.var_names, exact, elem_name=fmt_name("var_names"))
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
    for attr in [
        "X",
        "obs",
        "var",
        "obsm",
        "varm",
        "layers",
        "uns",
        "obsp",
        "varp",
        "raw",
    ]:
        assert_equal(
            getattr(a, attr),
            getattr(b, attr),
            exact,
            elem_name=fmt_name(attr),
        )


@singledispatch
def as_dense_dask_array(a):
    import dask.array as da

    return da.asarray(a)


@as_dense_dask_array.register(sparse.spmatrix)
def _(a):
    return as_dense_dask_array(a.toarray())


@contextmanager
def pytest_8_raises(exc_cls, *, match: str | re.Pattern = None):
    """Error handling using pytest 8's support for __notes__.

    See: https://github.com/pytest-dev/pytest/pull/11227

    Remove once pytest 8 is out!
    """

    with pytest.raises(exc_cls) as exc_info:
        yield exc_info

    check_error_or_notes_match(exc_info, match)


def check_error_or_notes_match(e: pytest.ExceptionInfo, pattern: str | re.Pattern):
    """
    Checks whether the printed error message or the notes contains the given pattern.

    DOES NOT WORK IN IPYTHON - because of the way IPython handles exceptions
    """
    import traceback

    message = "".join(traceback.format_exception_only(e.type, e.value))
    assert re.search(
        pattern, message
    ), f"Could not find pattern: '{pattern}' in error:\n\n{message}\n"


def as_cupy_type(val, typ=None):
    """
    Rough conversion function

    Will try to infer target type from input type if not specified.
    """
    if typ is None:
        input_typ = type(val)
        if issubclass(input_typ, np.ndarray):
            typ = CupyArray
        elif issubclass(input_typ, sparse.csr_matrix):
            typ = CupyCSRMatrix
        elif issubclass(input_typ, sparse.csc_matrix):
            typ = CupyCSCMatrix
        else:
            raise NotImplementedError(
                f"No default target type for input type {input_typ}"
            )

    if issubclass(typ, CupyArray):
        import cupy as cp

        if isinstance(val, sparse.spmatrix):
            val = val.toarray()
        return cp.array(val)
    elif issubclass(typ, CupyCSRMatrix):
        import cupyx.scipy.sparse as cpsparse
        import cupy as cp

        if isinstance(val, np.ndarray):
            return cpsparse.csr_matrix(cp.array(val))
        else:
            return cpsparse.csr_matrix(val)
    elif issubclass(typ, CupyCSCMatrix):
        import cupyx.scipy.sparse as cpsparse
        import cupy as cp

        if isinstance(val, np.ndarray):
            return cpsparse.csc_matrix(cp.array(val))
        else:
            return cpsparse.csc_matrix(val)
    else:
        raise NotImplementedError(
            f"Conversion from {type(val)} to {typ} not implemented"
        )


BASE_MATRIX_PARAMS = [
    pytest.param(asarray, id="np_array"),
    pytest.param(sparse.csr_matrix, id="scipy_csr"),
    pytest.param(sparse.csc_matrix, id="scipy_csc"),
]

DASK_MATRIX_PARAMS = [
    pytest.param(as_dense_dask_array, id="dask_array"),
]

CUPY_MATRIX_PARAMS = [
    pytest.param(
        partial(as_cupy_type, typ=CupyArray), id="cupy_array", marks=pytest.mark.gpu
    ),
    pytest.param(
        partial(as_cupy_type, typ=CupyCSRMatrix),
        id="cupy_csr",
        marks=pytest.mark.gpu,
    ),
    pytest.param(
        partial(as_cupy_type, typ=CupyCSCMatrix),
        id="cupy_csc",
        marks=pytest.mark.gpu,
    ),
]
