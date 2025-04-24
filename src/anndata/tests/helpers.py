from __future__ import annotations

import itertools
import random
import warnings
from collections import Counter, defaultdict
from collections.abc import Mapping
from functools import partial, singledispatch, wraps
from string import ascii_letters
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
from pandas.api.types import is_numeric_dtype
from scipy import sparse

import anndata
from anndata import AnnData, ExperimentalFeatureWarning, Raw
from anndata._core.aligned_mapping import AlignedMappingBase
from anndata._core.sparse_dataset import BaseCompressedSparseDataset
from anndata._core.views import ArrayView
from anndata.compat import (
    AwkArray,
    CSArray,
    CSMatrix,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    CupySparseMatrix,
    DaskArray,
    XDataArray,
    XDataset,
    ZarrArray,
    is_zarr_v2,
)
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Iterable
    from typing import Literal, TypeGuard, TypeVar

    from zarr.abc.store import ByteRequest
    from zarr.core.buffer import BufferPrototype

    from .._types import ArrayStorageType

    DT = TypeVar("DT")


try:
    from pandas.core.arrays.integer import IntegerDtype
except ImportError:
    IntegerDtype = (
        *(pd.Int8Dtype, pd.Int16Dtype, pd.Int32Dtype, pd.Int64Dtype),
        *(pd.UInt8Dtype, pd.UInt16Dtype, pd.UInt32Dtype, pd.UInt64Dtype),
    )


DEFAULT_KEY_TYPES = (
    sparse.csr_matrix,
    np.ndarray,
    pd.DataFrame,
    sparse.csr_array,
)


DEFAULT_COL_TYPES = (
    pd.CategoricalDtype(ordered=False),
    pd.CategoricalDtype(ordered=True),
    np.int64,
    np.float64,
    np.uint8,
    np.bool_,
    pd.BooleanDtype,
    pd.Int32Dtype,
)


# Give this to gen_adata when dask array support is expected.
GEN_ADATA_DASK_ARGS = dict(
    obsm_types=(*DEFAULT_KEY_TYPES, DaskArray),
    varm_types=(*DEFAULT_KEY_TYPES, DaskArray),
    layers_types=(*DEFAULT_KEY_TYPES, DaskArray),
)

GEN_ADATA_NO_XARRAY_ARGS = dict(
    obsm_types=(*DEFAULT_KEY_TYPES, AwkArray), varm_types=(*DEFAULT_KEY_TYPES, AwkArray)
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


def issubdtype(
    a: np.dtype | pd.api.extensions.ExtensionDtype | type,
    b: type[DT] | tuple[type[DT], ...],
) -> TypeGuard[DT]:
    if isinstance(b, tuple):
        return any(issubdtype(a, t) for t in b)
    if isinstance(a, type) and issubclass(a, pd.api.extensions.ExtensionDtype):
        return issubclass(a, b)
    if isinstance(a, pd.api.extensions.ExtensionDtype):
        return isinstance(a, b)
    try:
        return np.issubdtype(a, b)
    except TypeError:  # pragma: no cover
        pytest.fail(f"issubdtype can’t handle everything yet: {a} {b}")


def gen_random_column(  # noqa: PLR0911
    n: int, dtype: np.dtype | pd.api.extensions.ExtensionDtype
) -> tuple[str, np.ndarray | pd.api.extensions.ExtensionArray]:
    if issubdtype(dtype, pd.CategoricalDtype):
        # TODO: Think about allowing index to be passed for n
        letters = np.fromiter(iter(ascii_letters), "U1")
        if n > len(letters):
            letters = letters[: n // 2]  # Make sure categories are repeated
        key = "cat" if dtype.ordered else "cat_unordered"
        return key, pd.Categorical(np.random.choice(letters, n), dtype=dtype)
    if issubdtype(dtype, pd.BooleanDtype):
        return (
            "nullable-bool",
            pd.arrays.BooleanArray(
                np.random.randint(0, 2, size=n, dtype=bool),
                mask=np.random.randint(0, 2, size=n, dtype=bool),
            ),
        )
    if issubdtype(dtype, IntegerDtype):
        return (
            "nullable-int",
            pd.arrays.IntegerArray(
                np.random.randint(0, 1000, size=n, dtype=np.int32),
                mask=np.random.randint(0, 2, size=n, dtype=bool),
            ),
        )
    if issubdtype(dtype, pd.StringDtype):
        letters = np.fromiter(iter(ascii_letters), "U1")
        array = pd.array(np.random.choice(letters, n), dtype=pd.StringDtype())
        array[np.random.randint(0, 2, size=n, dtype=bool)] = pd.NA
        return "string", array
    # if issubdtype(dtype, pd.DatetimeTZDtype):
    #    return "datetime", pd.to_datetime(np.random.randint(0, 1000, size=n))
    if issubdtype(dtype, np.bool_):
        return "bool", np.random.randint(0, 2, size=n, dtype=dtype)

    if not issubdtype(dtype, np.number):  # pragma: no cover
        pytest.fail(f"Unexpected dtype: {dtype}")

    n_bits = 8 * (dtype().itemsize if isinstance(dtype, type) else dtype.itemsize)

    if issubdtype(dtype, np.unsignedinteger):
        return f"uint{n_bits}", np.random.randint(0, 255, n, dtype=dtype)
    if issubdtype(dtype, np.signedinteger):
        return f"int{n_bits}", np.random.randint(-50, 50, n, dtype=dtype)
    if issubdtype(dtype, np.floating):
        return f"float{n_bits}", np.random.random(n).astype(dtype)

    pytest.fail(f"Unexpected numeric dtype: {dtype}")  # pragma: no cover


def gen_typed_df(
    n: int,
    index: pd.Index[str] | None = None,
    dtypes: Collection[np.dtype | pd.api.extensions.ExtensionDtype] = DEFAULT_COL_TYPES,
):
    columns = [gen_random_column(n, dtype) for dtype in dtypes]
    col_names = [n for n, _ in columns]
    assert len(col_names) == len(set(col_names)), "Duplicate column names generated!"
    return pd.DataFrame(dict(columns), index=index)


def _gen_awkward_inner(shape, rng, dtype):
    # the maximum length a ragged dimension can take
    MAX_RAGGED_DIM_LEN = 20
    if not len(shape):
        # abort condition -> no dimension left, return an actual value instead
        return dtype(rng.randrange(1000))
    else:
        curr_dim_len = shape[0]
        if curr_dim_len is None:
            # ragged dimension, set random length
            curr_dim_len = rng.randrange(MAX_RAGGED_DIM_LEN)

        return [_gen_awkward_inner(shape[1:], rng, dtype) for _ in range(curr_dim_len)]


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
        msg = "The first dimension must be fixed-length."
        raise ValueError(msg)

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


def maybe_add_sparse_array(
    mapping: Mapping,
    types: Collection[type],
    format: Literal["csr", "csc"],
    random_state: np.random.Generator,
    shape: tuple[int, int],
):
    if sparse.csr_array in types or sparse.csr_matrix in types:
        mapping["sparse_array"] = sparse.csr_array(
            sparse.random(*shape, format=format, random_state=random_state)
        )
    return mapping


# TODO: Use hypothesis for this?
def gen_adata(  # noqa: PLR0913
    shape: tuple[int, int],
    X_type: Callable[[np.ndarray], object] = sparse.csr_matrix,
    *,
    X_dtype: np.dtype = np.float32,
    obs_dtypes: Collection[
        np.dtype | pd.api.extensions.ExtensionDtype
    ] = DEFAULT_COL_TYPES,
    var_dtypes: Collection[
        np.dtype | pd.api.extensions.ExtensionDtype
    ] = DEFAULT_COL_TYPES,
    obs_xdataset: bool = False,
    var_xdataset: bool = False,
    obsm_types: Collection[type] = (*DEFAULT_KEY_TYPES, AwkArray, XDataset),
    varm_types: Collection[type] = (*DEFAULT_KEY_TYPES, AwkArray, XDataset),
    layers_types: Collection[type] = DEFAULT_KEY_TYPES,
    random_state: np.random.Generator | None = None,
    sparse_fmt: Literal["csr", "csc"] = "csr",
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
    import xarray as xr

    if random_state is None:
        random_state = np.random.default_rng()

    M, N = shape
    obs_names = pd.Index(f"cell{i}" for i in range(shape[0]))
    var_names = pd.Index(f"gene{i}" for i in range(shape[1]))
    obs = gen_typed_df(M, obs_names, dtypes=obs_dtypes)
    var = gen_typed_df(N, var_names, dtypes=var_dtypes)
    # For #147
    obs.rename(columns=dict(cat="obs_cat"), inplace=True)
    var.rename(columns=dict(cat="var_cat"), inplace=True)

    if obs_xdataset:
        obs = XDataset.from_dataframe(obs)
    if var_xdataset:
        var = XDataset.from_dataframe(var)

    if X_type is None:
        X = None
    else:
        X = X_type(random_state.binomial(100, 0.005, (M, N)).astype(X_dtype))

    obsm = dict(
        array=np.random.random((M, 50)),
        sparse=sparse.random(M, 100, format=sparse_fmt, random_state=random_state),
        df=gen_typed_df(M, obs_names, dtypes=obs_dtypes),
        awk_2d_ragged=gen_awkward((M, None)),
        da=da.random.random((M, 50)),
        xdataset=xr.Dataset.from_dataframe(
            gen_typed_df(M, obs_names, dtypes=obs_dtypes)
        ),
    )
    obsm = {k: v for k, v in obsm.items() if type(v) in obsm_types}
    obsm = maybe_add_sparse_array(
        mapping=obsm,
        types=obsm_types,
        format=sparse_fmt,
        random_state=random_state,
        shape=(M, 100),
    )
    varm = dict(
        array=np.random.random((N, 50)),
        sparse=sparse.random(N, 100, format=sparse_fmt, random_state=random_state),
        df=gen_typed_df(N, var_names, dtypes=var_dtypes),
        awk_2d_ragged=gen_awkward((N, None)),
        da=da.random.random((N, 50)),
        xdataset=xr.Dataset.from_dataframe(
            gen_typed_df(N, var_names, dtypes=var_dtypes)
        ),
    )
    varm = {k: v for k, v in varm.items() if type(v) in varm_types}
    varm = maybe_add_sparse_array(
        mapping=varm,
        types=varm_types,
        format=sparse_fmt,
        random_state=random_state,
        shape=(N, 100),
    )
    layers = dict(
        array=np.random.random((M, N)),
        sparse=sparse.random(M, N, format=sparse_fmt, random_state=random_state),
        da=da.random.random((M, N)),
    )
    layers = maybe_add_sparse_array(
        mapping=layers,
        types=layers_types,
        format=sparse_fmt,
        random_state=random_state,
        shape=(M, N),
    )
    layers = {k: v for k, v in layers.items() if type(v) in layers_types}
    obsp = dict(
        array=np.random.random((M, M)),
        sparse=sparse.random(M, M, format=sparse_fmt, random_state=random_state),
    )
    obsp["sparse_array"] = sparse.csr_array(
        sparse.random(M, M, format=sparse_fmt, random_state=random_state)
    )
    varp = dict(
        array=np.random.random((N, N)),
        sparse=sparse.random(N, N, format=sparse_fmt, random_state=random_state),
    )
    varp["sparse_array"] = sparse.csr_array(
        sparse.random(N, N, format=sparse_fmt, random_state=random_state)
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
    # https://github.com/zarr-developers/zarr-python/issues/2134
    # zarr v3 on-disk does not write structured dtypes
    if anndata.settings.zarr_write_format == 3:
        del uns["O_recarray"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ExperimentalFeatureWarning)
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


def list_bool_subset(index, min_size=2):
    return array_bool_subset(index, min_size=min_size).tolist()


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


def sparray_bool_subset(index, min_size=2):
    return sparse.csr_array(
        array_bool_subset(index, min_size=min_size).reshape(len(index), 1)
    )


def array_subset(index, min_size=2):
    if len(index) < min_size:
        msg = f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        raise ValueError(msg)
    return np.random.choice(
        index, size=np.random.randint(min_size, len(index), ()), replace=False
    )


def array_int_subset(index, min_size=2):
    if len(index) < min_size:
        msg = f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        raise ValueError(msg)
    return np.random.choice(
        np.arange(len(index)),
        size=np.random.randint(min_size, len(index), ()),
        replace=False,
    )


def list_int_subset(index, min_size=2):
    return array_int_subset(index, min_size=min_size).tolist()


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
        list_int_subset,
        array_bool_subset,
        list_bool_subset,
        matrix_bool_subset,
        spmatrix_bool_subset,
        sparray_bool_subset,
    ]
)
def subset_func(request):
    return request.param


###################
# Checking equality
###################


def format_msg(elem_name: str | None) -> str:
    if elem_name is not None:
        return f"Error raised from element {elem_name!r}."
    else:
        return ""


# TODO: it would be better to modify the other exception
def report_name(func):
    """Report name of element being tested if test fails."""

    @wraps(func)
    def func_wrapper(*args, _elem_name: str | None = None, **kwargs):
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
def assert_equal(
    a: object, b: object, *, exact: bool = False, elem_name: str | None = None
):
    _assert_equal(a, b, _elem_name=elem_name)


@assert_equal.register(CupyArray)
def assert_equal_cupy(
    a: CupyArray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    assert_equal(b, a.get(), exact=exact, elem_name=elem_name)


@assert_equal.register(np.ndarray)
def assert_equal_ndarray(
    a: np.ndarray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    b = asarray(b)
    if not exact and is_numeric_dtype(a) and is_numeric_dtype(b):
        assert a.shape == b.shape, format_msg(elem_name)
        np.testing.assert_allclose(a, b, equal_nan=True, err_msg=format_msg(elem_name))
    elif (  # Structured dtype
        not exact
        and hasattr(a, "dtype")
        and hasattr(b, "dtype")
        and len(a.dtype) > 1
        and len(b.dtype) > 0
    ):
        # Reshaping to allow >2d arrays
        assert a.shape == b.shape, format_msg(elem_name)
        assert_equal(
            pd.DataFrame(a.reshape(-1)),
            pd.DataFrame(b.reshape(-1)),
            exact=exact,
            elem_name=elem_name,
        )
    else:
        assert np.all(a == b), format_msg(elem_name)


@assert_equal.register(ArrayView)
def assert_equal_arrayview(
    a: ArrayView, b: object, *, exact: bool = False, elem_name: str | None = None
):
    assert_equal(asarray(a), asarray(b), exact=exact, elem_name=elem_name)


@assert_equal.register(BaseCompressedSparseDataset)
@assert_equal.register(sparse.spmatrix)
def assert_equal_sparse(
    a: BaseCompressedSparseDataset | sparse.spmatrix,
    b: object,
    *,
    exact: bool = False,
    elem_name: str | None = None,
):
    a = asarray(a)
    assert_equal(b, a, exact=exact, elem_name=elem_name)


@assert_equal.register(CSArray)
def assert_equal_sparse_array(
    a: CSArray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    return assert_equal_sparse(a, b, exact=exact, elem_name=elem_name)


@assert_equal.register(CupySparseMatrix)
def assert_equal_cupy_sparse(
    a: CupySparseMatrix, b: object, *, exact: bool = False, elem_name: str | None = None
):
    a = a.toarray()
    assert_equal(b, a, exact=exact, elem_name=elem_name)


@assert_equal.register(h5py.Dataset)
@assert_equal.register(ZarrArray)
def assert_equal_h5py_dataset(
    a: ArrayStorageType, b: object, *, exact: bool = False, elem_name: str | None = None
):
    a = asarray(a)
    assert_equal(b, a, exact=exact, elem_name=elem_name)


@assert_equal.register(DaskArray)
def assert_equal_dask_array(
    a: DaskArray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    assert_equal(b, a.compute(), exact=exact, elem_name=elem_name)


@assert_equal.register(pd.DataFrame)
def are_equal_dataframe(
    a: pd.DataFrame, b: object, *, exact: bool = False, elem_name: str | None = None
):
    if not isinstance(b, pd.DataFrame):
        assert_equal(b, a, exact=exact, elem_name=elem_name)  # , a.values maybe?

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
def assert_equal_awkarray(
    a: AwkArray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    import awkward as ak

    if exact:
        assert isinstance(b, AwkArray)
        assert a.type == b.type, f"{a.type} != {b.type}, {format_msg(elem_name)}"
    assert ak.to_list(a) == ak.to_list(b), format_msg(elem_name)


@assert_equal.register(Mapping)
def assert_equal_mapping(
    a: Mapping, b: object, *, exact: bool = False, elem_name: str | None = None
):
    assert isinstance(b, Mapping)
    assert set(a) == set(b), format_msg(elem_name)
    for k in a:
        if elem_name is None:
            elem_name = ""
        assert_equal(a[k], b[k], exact=exact, elem_name=f"{elem_name}/{k}")


@assert_equal.register(AlignedMappingBase)
def assert_equal_aligned_mapping(
    a: AlignedMappingBase,
    b: object,
    *,
    exact: bool = False,
    elem_name: str | None = None,
):
    assert isinstance(b, AlignedMappingBase)
    a_indices = (a.parent.obs_names, a.parent.var_names)
    b_indices = (b.parent.obs_names, b.parent.var_names)
    for axis_idx in a.axes:
        assert_equal(
            a_indices[axis_idx], b_indices[axis_idx], exact=exact, elem_name=axis_idx
        )
    assert a.attrname == b.attrname, format_msg(elem_name)
    assert_equal_mapping(a, b, exact=exact, elem_name=elem_name)


@assert_equal.register(pd.Index)
def assert_equal_index(
    a: pd.Index, b: object, *, exact: bool = False, elem_name: str | None = None
):
    params = dict(check_categorical=False) if not exact else {}
    report_name(pd.testing.assert_index_equal)(
        a, b, check_names=False, **params, _elem_name=elem_name
    )


@assert_equal.register(pd.api.extensions.ExtensionArray)
def assert_equal_extension_array(
    a: pd.api.extensions.ExtensionArray,
    b: object,
    *,
    exact: bool = False,
    elem_name: str | None = None,
):
    report_name(pd.testing.assert_extension_array_equal)(
        a,
        b,
        check_dtype=exact,
        check_exact=exact,
        _elem_name=elem_name,
    )


@assert_equal.register(XDataArray)
def assert_equal_xarray(
    a: XDataArray, b: object, *, exact: bool = False, elem_name: str | None = None
):
    report_name(a.equals)(b, _elem_name=elem_name)


@assert_equal.register(Raw)
def assert_equal_raw(
    a: Raw, b: object, *, exact: bool = False, elem_name: str | None = None
):
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
    a: AnnData, b: object, *, exact: bool = False, elem_name: str | None = None
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

    assert isinstance(b, AnnData)

    # There may be issues comparing views, since np.allclose
    # can modify ArrayViews if they contain `nan`s
    assert_equal(a.obs_names, b.obs_names, exact=exact, elem_name=fmt_name("obs_names"))
    assert_equal(a.var_names, b.var_names, exact=exact, elem_name=fmt_name("var_names"))
    if not exact:
        # Reorder all elements if necessary
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
            exact=exact,
            elem_name=fmt_name(attr),
        )


def _half_chunk_size(a: tuple[int, ...]) -> tuple[int, ...]:
    def half_rounded_up(x):
        div, mod = divmod(x, 2)
        return div + (mod > 0)

    return tuple(half_rounded_up(x) for x in a)


@singledispatch
def as_dense_dask_array(a):
    import dask.array as da

    a = asarray(a)
    return da.asarray(a, chunks=_half_chunk_size(a.shape))


@as_dense_dask_array.register(CSMatrix)
def _(a):
    return as_dense_dask_array(a.toarray())


@as_dense_dask_array.register(DaskArray)
def _(a):
    return a.map_blocks(asarray, dtype=a.dtype, meta=np.ndarray)


@singledispatch
def as_sparse_dask_array(a) -> DaskArray:
    import dask.array as da

    return da.from_array(sparse.csr_matrix(a), chunks=_half_chunk_size(a.shape))


@as_sparse_dask_array.register(CSMatrix)
def _(a):
    import dask.array as da

    return da.from_array(a, _half_chunk_size(a.shape))


@as_sparse_dask_array.register(CSArray)
def _(a):
    import dask.array as da

    return da.from_array(sparse.csr_matrix(a), _half_chunk_size(a.shape))


@as_sparse_dask_array.register(DaskArray)
def _(a):
    return a.map_blocks(sparse.csr_matrix)


@singledispatch
def as_dense_cupy_dask_array(a):
    import cupy as cp

    return as_dense_dask_array(a).map_blocks(
        cp.array, meta=cp.array((1.0), dtype=a.dtype), dtype=a.dtype
    )


@as_dense_cupy_dask_array.register(CupyArray)
def _(a):
    import cupy as cp
    import dask.array as da

    return da.from_array(
        a,
        chunks=_half_chunk_size(a.shape),
        meta=cp.array((1.0), dtype=a.dtype),
    )


@as_dense_cupy_dask_array.register(DaskArray)
def _(a):
    import cupy as cp

    if isinstance(a._meta, cp.ndarray):
        return a.copy()
    return a.map_blocks(
        partial(as_cupy, typ=CupyArray),
        dtype=a.dtype,
        meta=cp.array((1.0), dtype=a.dtype),
    )


try:
    import cupyx.scipy.sparse as cpsparse

    format_to_memory_class = {"csr": cpsparse.csr_matrix, "csc": cpsparse.csc_matrix}
except ImportError:
    format_to_memory_class = {}


# TODO: If there are chunks which divide along columns, then a coo_matrix is returned by compute
# We should try and fix this upstream in dask/ cupy
@singledispatch
def as_cupy_sparse_dask_array(a, format="csr"):
    memory_class = format_to_memory_class[format]
    cpu_da = as_sparse_dask_array(a)
    return cpu_da.rechunk((cpu_da.chunks[0], -1)).map_blocks(
        memory_class, dtype=a.dtype, meta=memory_class(cpu_da._meta)
    )


@as_cupy_sparse_dask_array.register(CupyArray)
@as_cupy_sparse_dask_array.register(CupySparseMatrix)
def _(a, format="csr"):
    import dask.array as da

    memory_class = format_to_memory_class[format]
    return da.from_array(memory_class(a), chunks=(_half_chunk_size(a.shape)[0], -1))


@as_cupy_sparse_dask_array.register(DaskArray)
def _(a, format="csr"):
    memory_class = format_to_memory_class[format]
    if isinstance(a._meta, memory_class):
        return a.copy()
    return a.rechunk((a.chunks[0], -1)).map_blocks(
        partial(as_cupy, typ=memory_class), dtype=a.dtype
    )


def resolve_cupy_type(val):
    input_typ = type(val) if not isinstance(val, type) else val

    if issubclass(input_typ, np.ndarray):
        typ = CupyArray
    elif issubclass(input_typ, sparse.csr_matrix):
        typ = CupyCSRMatrix
    elif issubclass(input_typ, sparse.csc_matrix):
        typ = CupyCSCMatrix
    else:
        msg = f"No default target type for input type {input_typ}"
        raise NotImplementedError(msg)
    return typ


@singledispatch
def as_cupy(val, typ=None):
    """
    Rough conversion function

    Will try to infer target type from input type if not specified.
    """
    if typ is None:
        typ = resolve_cupy_type(val)

    if issubclass(typ, CupyArray):
        import cupy as cp

        if isinstance(val, CSMatrix):
            val = val.toarray()
        return cp.array(val)
    elif issubclass(typ, CupyCSRMatrix):
        import cupy as cp
        import cupyx.scipy.sparse as cpsparse

        if isinstance(val, np.ndarray):
            return cpsparse.csr_matrix(cp.array(val))
        else:
            return cpsparse.csr_matrix(val)
    elif issubclass(typ, CupyCSCMatrix):
        import cupy as cp
        import cupyx.scipy.sparse as cpsparse

        if isinstance(val, np.ndarray):
            return cpsparse.csc_matrix(cp.array(val))
        else:
            return cpsparse.csc_matrix(val)
    else:
        msg = f"Conversion from {type(val)} to {typ} not implemented"
        raise NotImplementedError(msg)


# TODO: test
@as_cupy.register(DaskArray)
def as_cupy_dask(a, typ=None):
    if typ is None:
        typ = resolve_cupy_type(a._meta)
    return a.map_blocks(partial(as_cupy, typ=typ), dtype=a.dtype)


@singledispatch
def shares_memory(x, y) -> bool:
    return np.shares_memory(x, y)


@shares_memory.register(CSMatrix)
def shares_memory_sparse(x, y):
    return (
        np.shares_memory(x.data, y.data)
        and np.shares_memory(x.indices, y.indices)
        and np.shares_memory(x.indptr, y.indptr)
    )


BASE_MATRIX_PARAMS = [
    pytest.param(asarray, id="np_array"),
    pytest.param(sparse.csr_matrix, id="scipy_csr_matrix"),
    pytest.param(sparse.csc_matrix, id="scipy_csc_matrix"),
    pytest.param(sparse.csr_array, id="scipy_csr_array"),
    pytest.param(sparse.csc_array, id="scipy_csc_array"),
]

DASK_MATRIX_PARAMS = [
    pytest.param(as_dense_dask_array, id="dense_dask_array"),
    pytest.param(as_sparse_dask_array, id="sparse_dask_array"),
]

CUPY_MATRIX_PARAMS = [
    pytest.param(
        partial(as_cupy, typ=CupyArray), id="cupy_array", marks=pytest.mark.gpu
    ),
    pytest.param(
        partial(as_cupy, typ=CupyCSRMatrix),
        id="cupy_csr",
        marks=pytest.mark.gpu,
    ),
    pytest.param(
        partial(as_cupy, typ=CupyCSCMatrix),
        id="cupy_csc",
        marks=pytest.mark.gpu,
    ),
]

DASK_CUPY_MATRIX_PARAMS = [
    pytest.param(
        as_dense_cupy_dask_array,
        id="cupy_dense_dask_array",
        marks=pytest.mark.gpu,
    ),
    pytest.param(
        as_cupy_sparse_dask_array, id="cupy_csr_dask_array", marks=pytest.mark.gpu
    ),
]

if is_zarr_v2():
    from zarr.storage import DirectoryStore as LocalStore
else:
    from zarr.storage import LocalStore


class AccessTrackingStoreBase(LocalStore):
    _access_count: Counter[str]
    _accessed: defaultdict[str, set]
    _accessed_keys: defaultdict[str, list[str]]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._access_count = Counter()
        self._accessed = defaultdict(set)
        self._accessed_keys = defaultdict(list)

    def _check_and_track_key(self, key: str):
        for tracked in self._access_count:
            if tracked in key:
                self._access_count[tracked] += 1
                self._accessed[tracked].add(key)
                self._accessed_keys[tracked] += [key]

    def get_access_count(self, key: str) -> int:
        # access defaultdict when value is not there causes key to be there,
        # which causes it to be tracked
        if key not in self._access_count:
            msg = f"{key} not found among access count"
            raise KeyError(msg)
        return self._access_count[key]

    def get_subkeys_accessed(self, key: str) -> set[str]:
        if key not in self._accessed:
            msg = f"{key} not found among accessed"
            raise KeyError(msg)
        return self._accessed[key]

    def get_accessed_keys(self, key: str) -> list[str]:
        if key not in self._accessed_keys:
            msg = f"{key} not found among accessed keys"
            raise KeyError(msg)
        return self._accessed_keys[key]

    def initialize_key_trackers(self, keys_to_track: Iterable[str]) -> None:
        for k in keys_to_track:
            self._access_count[k] = 0
            self._accessed_keys[k] = []
            self._accessed[k] = set()

    def reset_key_trackers(self) -> None:
        self.initialize_key_trackers(self._access_count.keys())

    def assert_access_count(self, key: str, count: int):
        keys_accessed = self.get_subkeys_accessed(key)
        access_count = self.get_access_count(key)
        assert self.get_access_count(key) == count, (
            f"Found {access_count} accesses at {keys_accessed}"
        )


if is_zarr_v2():

    class AccessTrackingStore(AccessTrackingStoreBase):
        def __getitem__(self, key: str) -> bytes:
            self._check_and_track_key(key)
            return super().__getitem__(key)

else:

    class AccessTrackingStore(AccessTrackingStoreBase):
        async def get(
            self,
            key: str,
            prototype: BufferPrototype | None = None,
            byte_range: ByteRequest | None = None,
        ) -> object:
            self._check_and_track_key(key)
            return await super().get(key, prototype=prototype, byte_range=byte_range)


if is_zarr_v2():

    class AccessTrackingStore(AccessTrackingStoreBase):
        def __getitem__(self, key: str) -> bytes:
            self._check_and_track_key(key)
            return super().__getitem__(key)

else:

    class AccessTrackingStore(AccessTrackingStoreBase):
        async def get(
            self,
            key: str,
            prototype: BufferPrototype | None = None,
            byte_range: ByteRequest | None = None,
        ) -> object:
            self._check_and_track_key(key)
            return await super().get(key, prototype=prototype, byte_range=byte_range)


def get_multiindex_columns_df(shape: tuple[int, int]) -> pd.DataFrame:
    return pd.DataFrame(
        np.random.rand(shape[0], shape[1]),
        columns=pd.MultiIndex.from_tuples(
            list(itertools.product(["a"], range(shape[1] - (shape[1] // 2))))
            + list(itertools.product(["b"], range(shape[1] // 2)))
        ),
    )
