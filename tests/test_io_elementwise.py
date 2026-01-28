"""
Tests that each element in an anndata is written correctly
"""

from __future__ import annotations

import re
from contextlib import ExitStack, nullcontext
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, cast

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from packaging.version import Version
from scipy import sparse

import anndata as ad
from anndata._io.specs import _REGISTRY, IOSpec, get_spec
from anndata._io.specs.registry import IORegistryError
from anndata._io.zarr import open_write_group
from anndata.compat import CSArray, CSMatrix, H5Group, ZarrGroup, _read_attr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem, write_elem
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    as_cupy,
    as_cupy_sparse_dask_array,
    as_dense_cupy_dask_array,
    assert_equal,
    check_all_sharded,
    gen_adata,
    visititems_zarr,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path
    from typing import Literal

    from anndata._types import GroupStorageType
    from anndata.compat import H5Group


PANDAS_3 = Version(version("pandas")) >= Version("3rc0")


@pytest.fixture
def exit_stack() -> Generator[ExitStack, None, None]:
    with ExitStack() as stack:
        yield stack


@pytest.fixture
def store(diskfmt, tmp_path) -> Generator[H5Group | ZarrGroup, None, None]:
    if diskfmt == "h5ad":
        file = h5py.File(tmp_path / "test.h5ad", "w")
        store = cast("H5Group", file["/"])
    elif diskfmt == "zarr":
        store = open_write_group(tmp_path / "test.zarr")
    else:
        pytest.fail(f"Unknown store type: {diskfmt}")

    try:
        yield store
    finally:
        if diskfmt == "h5ad":
            file.close()


sparse_formats = ["csr", "csc"]
SIZE = 50
DEFAULT_SHAPE = (SIZE, SIZE * 2)


@pytest.fixture(params=sparse_formats)
def sparse_format(request: pytest.FixtureRequest) -> Literal["csr", "csc"]:
    return request.param


def create_dense_store(
    store: H5Group | ZarrGroup, *, shape: tuple[int, ...] = DEFAULT_SHAPE
) -> H5Group | ZarrGroup:
    X = np.random.randn(*shape)

    write_elem(store, "X", X)
    return store


def create_sparse_store[G: (H5Group, ZarrGroup)](
    sparse_format: Literal["csc", "csr"], store: G, shape=DEFAULT_SHAPE
) -> G:
    """Returns a store

    Parameters
    ----------
    sparse_format
    store

    Returns
    -------
        A store with a key, `X` that is simply a sparse matrix, and `X_dask` where that same array is wrapped by dask
    """
    import dask.array as da

    X = sparse.random(
        shape[0],
        shape[1],
        format=sparse_format,
        density=0.01,
        random_state=np.random.default_rng(),
    )
    X_dask = da.from_array(
        X,
        chunks=(100 if format == "csr" else SIZE, SIZE * 2 if format == "csr" else 100),
    )

    write_elem(store, "X", X)
    write_elem(store, "X_dask", X_dask)
    return store


@pytest.mark.parametrize(
    ("value", "encoding_type"),
    [
        pytest.param(None, "null", id="none"),
        pytest.param("hello world", "string", id="py_str"),
        pytest.param(np.str_("hello world"), "string", id="np_str"),
        pytest.param(np.array([1, 2, 3]), "array", id="np_arr_int"),
        pytest.param(
            np.array(["hello", "world"], dtype=object), "string-array", id="np_arr_str"
        ),
        pytest.param(1, "numeric-scalar", id="py_int"),
        pytest.param(True, "numeric-scalar", id="py_bool"),
        pytest.param(1.0, "numeric-scalar", id="py_float"),
        pytest.param({"a": 1}, "dict", id="py_dict"),
        pytest.param(
            lambda: gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS),
            "anndata",
            id="anndata",
        ),
        pytest.param(
            sparse.random(5, 3, format="csr", density=0.5),
            "csr_matrix",
            id="sp_mat_csr",
        ),
        pytest.param(
            sparse.random(5, 3, format="csc", density=0.5),
            "csc_matrix",
            id="sp_mat_csc",
        ),
        pytest.param(pd.DataFrame({"a": [1, 2, 3]}), "dataframe", id="pd_df"),
        pytest.param(
            pd.Categorical([*"aabccedd", pd.NA]),
            "categorical",
            id="pd_cat_np_str",
        ),
        pytest.param(
            pd.Categorical(list("aabccedd"), ordered=True),
            "categorical",
            id="pd_cat_np_str_ord",
        ),
        pytest.param(
            pd.array([*"aabccedd", pd.NA], dtype="string").astype("category"),
            "categorical",
            id="pd_cat_pd_str",
        ),
        pytest.param(
            pd.Categorical([1, 2, 1, 3], ordered=True), "categorical", id="pd_cat_num"
        ),
        pytest.param(
            pd.array(["hello", "world"], dtype="string"),
            "nullable-string-array",
            id="pd_arr_str",
        ),
        pytest.param(
            pd.array(["hello", "world", pd.NA], dtype="string"),
            "nullable-string-array",
            id="pd_arr_str_mask",
        ),
        pytest.param(
            pd.arrays.IntegerArray(
                np.ones(5, dtype=int), mask=np.array([True, False, True, False, True])
            ),
            "nullable-integer",
            id="pd_arr_int_mask",
        ),
        pytest.param(pd.array([1, 2, 3]), "nullable-integer", id="pd_arr_int"),
        pytest.param(
            pd.arrays.BooleanArray(
                np.random.randint(0, 2, size=5, dtype=bool),
                mask=np.random.randint(0, 2, size=5, dtype=bool),
            ),
            "nullable-boolean",
            id="pd_arr_bool_mask",
        ),
        pytest.param(
            pd.array([True, False, True, True]), "nullable-boolean", id="pd_arr_bool"
        ),
        pytest.param(
            zarr.ones((100, 100), chunks=(10, 10)),
            "array",
            id="zarr_dense_array",
        ),
        pytest.param(
            create_dense_store(
                h5py.File("test1.h5", mode="w", driver="core", backing_store=False)
            )["X"],
            "array",
            id="h5_dense_array",
        ),
        # pytest.param(bytes, b"some bytes", "bytes", id="py_bytes"), # Does not work for zarr
        # TODO consider how specific encodings should be. Should we be fully describing the written type?
        # Currently the info we add is: "what you wouldn't be able to figure out yourself"
        # but that's not really a solid rule.
        # pytest.param(bool, True, "bool", id="py_bool"),
        # pytest.param(bool, np.bool_(False), "bool", id="np_bool"),
    ],
)
def test_io_spec(store: GroupStorageType, value, encoding_type) -> None:
    if callable(value):
        value = value()

    key = f"key_for_{encoding_type}"
    with ad.settings.override(allow_write_nullable_strings=True):
        write_elem(store, key, value, dataset_kwargs={})

    assert encoding_type == _read_attr(store[key].attrs, "encoding-type")

    from_disk = read_elem(store[key])
    assert_equal(value, from_disk)
    assert get_spec(store[key]) == _REGISTRY.get_spec(value)


@pytest.mark.parametrize(
    ("value", "encoding_type"),
    [
        pytest.param(np.asarray(1), "numeric-scalar", id="scalar_int"),
        pytest.param(np.asarray(1.0), "numeric-scalar", id="scalar_float"),
        pytest.param(np.asarray(True), "numeric-scalar", id="scalar_bool"),  # noqa: FBT003
        pytest.param(np.asarray("test"), "string", id="scalar_string"),
    ],
)
def test_io_spec_compressed_scalars(
    store: H5Group | ZarrGroup, value: np.ndarray, encoding_type: str
):
    key = f"key_for_{encoding_type}"
    write_elem(
        store, key, value, dataset_kwargs={"compression": "gzip", "compression_opts": 5}
    )

    assert encoding_type == _read_attr(store[key].attrs, "encoding-type")

    from_disk = read_elem(store[key])
    assert_equal(value, from_disk)


# Can't instantiate cupy types at the top level, so converting them within the test
@pytest.mark.gpu
@pytest.mark.parametrize(
    ("value", "encoding_type"),
    [
        pytest.param(np.array([1, 2, 3]), "array", id="np-1d"),
        pytest.param(np.arange(12).reshape(4, 3), "array", id="np-2d"),
        pytest.param(
            sparse.random(5, 3, format="csr", density=0.5), "csr_matrix", id="csr"
        ),
        pytest.param(
            sparse.random(5, 3, format="csc", density=0.5), "csc_matrix", id="csc"
        ),
    ],
)
@pytest.mark.parametrize("as_dask", [False, True], ids=["local", "dask"])
def test_io_spec_cupy(store, value, encoding_type, as_dask):
    if not as_dask:
        value = as_cupy(value)
    elif isinstance(value, CSMatrix | CSArray):
        value = as_cupy_sparse_dask_array(value, format=encoding_type[:3])
    else:
        value = as_dense_cupy_dask_array(value)

    key = f"key_for_{encoding_type}"
    write_elem(store, key, value, dataset_kwargs={})

    assert encoding_type == _read_attr(store[key].attrs, "encoding-type")

    from_disk = as_cupy(read_elem(store[key]))
    assert_equal(value, from_disk)
    assert get_spec(store[key]) == _REGISTRY.get_spec(value)


def test_dask_write_sparse(sparse_format, store):
    x_sparse_store = create_sparse_store(sparse_format, store)
    X_from_disk = read_elem(x_sparse_store["X"])
    X_dask_from_disk = read_elem(x_sparse_store["X_dask"])

    assert_equal(X_from_disk, X_dask_from_disk)
    assert_equal(dict(x_sparse_store["X"].attrs), dict(x_sparse_store["X_dask"].attrs))

    assert x_sparse_store["X_dask/indptr"].dtype == np.int64
    assert x_sparse_store["X_dask/indices"].dtype == np.int64


def test_read_lazy_2d_dask(sparse_format, store):
    arr_store = create_sparse_store(sparse_format, store)
    X_dask_from_disk = read_elem_lazy(arr_store["X"])
    X_from_disk = read_elem(arr_store["X"])

    assert_equal(X_from_disk, X_dask_from_disk)
    random_int_indices = np.random.randint(0, SIZE, (SIZE // 10,))
    random_int_indices.sort()
    index_slice = slice(0, SIZE // 10)
    for index in [random_int_indices, index_slice]:
        assert_equal(X_from_disk[index, :], X_dask_from_disk[index, :])
        assert_equal(X_from_disk[:, index], X_dask_from_disk[:, index])
    random_bool_mask = np.random.randn(SIZE) > 0
    assert_equal(
        X_from_disk[random_bool_mask, :], X_dask_from_disk[random_bool_mask, :]
    )
    random_bool_mask = np.random.randn(SIZE * 2) > 0
    assert_equal(
        X_from_disk[:, random_bool_mask], X_dask_from_disk[:, random_bool_mask]
    )

    assert arr_store["X_dask/indptr"].dtype == np.int64
    assert arr_store["X_dask/indices"].dtype == np.int64


@pytest.mark.parametrize(
    ("n_dims", "chunks"),
    [
        (1, (10,)),
        (1, (40,)),
        (2, (10, 10)),
        (2, (40, 40)),
        (2, (20, 40)),
        (1, None),
        (2, None),
        (2, (40, -1)),
        (2, (40, None)),
    ],
)
def test_read_lazy_subsets_nd_dask(store: H5Group | ZarrGroup, n_dims, chunks) -> None:
    arr_store = create_dense_store(store, shape=DEFAULT_SHAPE[:n_dims])
    X_dask_from_disk = read_elem_lazy(arr_store["X"], chunks=chunks)
    X_from_disk = read_elem(arr_store["X"])
    assert_equal(X_from_disk, X_dask_from_disk)

    random_int_indices = np.random.randint(0, SIZE, (SIZE // 10,))
    random_int_indices.sort()
    random_bool_mask = np.random.randn(SIZE) > 0
    index_slice = slice(0, SIZE // 10)
    for index in [random_int_indices, index_slice, random_bool_mask]:
        assert_equal(X_from_disk[index], X_dask_from_disk[index])


@pytest.mark.xdist_group("dask")
@pytest.mark.dask_distributed
def test_read_lazy_h5_cluster(
    sparse_format: Literal["csr", "csc"],
    tmp_path: Path,
    local_cluster_addr: str,
) -> None:
    import dask.distributed as dd

    with h5py.File(tmp_path / "test.h5", "w") as file:
        store = file["/"]
        arr_store = create_sparse_store(sparse_format, store)
        X_dask_from_disk = read_elem_lazy(arr_store["X"])
        X_from_disk = read_elem(arr_store["X"])
    with dd.Client(local_cluster_addr):
        assert_equal(X_from_disk, X_dask_from_disk)


def test_undersized_shape_to_default(store: H5Group | ZarrGroup) -> None:
    shape = (1000, 50)
    arr_store = create_dense_store(store, shape=shape)
    X_dask_from_disk = read_elem_lazy(arr_store["X"])
    assert all(c <= s for c, s in zip(X_dask_from_disk.chunksize, shape, strict=True))
    assert X_dask_from_disk.shape == shape


@pytest.mark.parametrize(
    ("arr_type", "chunks", "expected_chunksize"),
    [
        ("dense", (10, 10), (10, 10)),
        ("csc", (SIZE, 1), (SIZE, 1)),
        ("csr", (1, SIZE * 2), (1, SIZE * 2)),
        ("csc", None, (SIZE, 100)),
        ("csr", None, DEFAULT_SHAPE),
        ("csr", (10, -1), (10, SIZE * 2)),
        ("csc", (-1, 10), (SIZE, 10)),
        ("csr", (10, None), (10, SIZE * 2)),
        ("csc", (None, 10), (SIZE, 10)),
        ("csc", (None, None), DEFAULT_SHAPE),
        ("csr", (None, None), DEFAULT_SHAPE),
        ("csr", (-1, -1), DEFAULT_SHAPE),
        ("csc", (-1, -1), DEFAULT_SHAPE),
    ],
)
def test_read_lazy_2d_chunk_kwargs(
    store: H5Group | ZarrGroup,
    arr_type: Literal["csr", "csc", "dense"],
    chunks: None | tuple[int | None, int | None],
    expected_chunksize: tuple[int, int],
) -> None:
    if arr_type == "dense":
        arr_store = create_dense_store(store)
        X_dask_from_disk = read_elem_lazy(arr_store["X"], chunks=chunks)
    else:
        arr_store = create_sparse_store(arr_type, store)
        X_dask_from_disk = read_elem_lazy(arr_store["X"], chunks=chunks)
    assert X_dask_from_disk.chunksize == expected_chunksize
    X_from_disk = read_elem(arr_store["X"])
    assert_equal(X_from_disk, X_dask_from_disk)


def test_read_lazy_bad_chunk_kwargs(tmp_path):
    arr_type = "csr"
    with h5py.File(tmp_path / "test.h5", "w") as file:
        store = file["/"]
        arr_store = create_sparse_store(arr_type, store)
        with pytest.raises(
            ValueError, match=r"`chunks` must be a tuple of two integers"
        ):
            read_elem_lazy(arr_store["X"], chunks=(SIZE,))
        with pytest.raises(ValueError, match=r"Only the major axis can be chunked"):
            read_elem_lazy(arr_store["X"], chunks=(SIZE, 10))


@pytest.mark.parametrize("sparse_format", ["csr", "csc"])
def test_write_indptr_dtype_override(store, sparse_format):
    X = sparse.random(
        100,
        100,
        format=sparse_format,
        density=0.1,
        random_state=np.random.default_rng(),
    )

    write_elem(store, "X", X, dataset_kwargs=dict(indptr_dtype="int64"))

    assert store["X/indptr"].dtype == np.int64
    assert X.indptr.dtype == np.int32
    np.testing.assert_array_equal(store["X/indptr"][...], X.indptr)


def test_io_spec_raw(store):
    adata = gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS)
    adata.raw = adata.copy()

    write_elem(store, "adata", adata)

    assert _read_attr(store["adata/raw"].attrs, "encoding-type") == "raw"

    from_disk = read_elem(store["adata"])
    assert_equal(from_disk.raw, adata.raw)


def test_write_anndata_to_root(store):
    adata = gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS)

    write_elem(store, "/", adata)
    # TODO: see https://github.com/zarr-developers/zarr-python/issues/2716
    if isinstance(store, ZarrGroup):
        store = zarr.open(store.store)
    from_disk = read_elem(store)

    assert _read_attr(store.attrs, "encoding-type") == "anndata"
    assert_equal(from_disk, adata)


@pytest.mark.parametrize(
    ("attribute", "value"),
    [
        ("encoding-type", "floob"),
        ("encoding-version", "10000.0"),
    ],
)
def test_read_iospec_not_found(store, attribute, value):
    adata = gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS)

    write_elem(store, "/", adata)
    store["obs"].attrs.update({attribute: value})

    with pytest.raises(IORegistryError) as exc_info:
        read_elem(store)
    msg = str(exc_info.value)

    assert "No read method registered for IOSpec" in msg
    assert f"{attribute.replace('-', '_')}='{value}'" in msg


@pytest.mark.parametrize(
    "obj",
    [(b"x",)],
)
def test_write_io_error(store, obj):
    full_pattern = re.compile(
        rf"No method registered for writing {type(obj)} into .*Group"
    )

    with pytest.raises(IORegistryError, match=r"while writing key '/el'") as exc_info:
        write_elem(store, "/el", obj)

    msg = str(exc_info.value)
    assert re.search(full_pattern, msg)


PAT_IMPLICIT = r"allow_write_nullable_strings.*None.*future\.infer_string.*False"


@pytest.mark.parametrize(
    ("ad_setting", "pd_setting", "expected_missing", "expected_no_missing"),
    [
        # when explicitly disabled, we expect an error when trying to write an array with missing values
        # and expect an array without missing values to be written in the old, non-nullable format
        *(
            pytest.param(
                False,
                pd_ignored,
                (ValueError, r"missing values.*allow_write_nullable_strings.*False"),
                "string-array",
                id=f"off-explicit-{int(pd_ignored)}",
            )
            for pd_ignored in [False, True]
        ),
        # when implicitly disabled, we expect a different message in both cases
        pytest.param(
            None,
            False,
            (RuntimeError, PAT_IMPLICIT),
            (RuntimeError, PAT_IMPLICIT),
            id="off-implicit",
        ),
        # when enabled, we expect arrays to be written in the nullable format
        pytest.param(None, True, *(["nullable-string-array"] * 2), id="on-implicit"),
        pytest.param(True, False, *(["nullable-string-array"] * 2), id="on-explicit-0"),
        pytest.param(True, True, *(["nullable-string-array"] * 2), id="on-explicit-1"),
    ],
)
@pytest.mark.parametrize("missing", [True, False], ids=["missing", "no_missing"])
def test_write_nullable_string(
    *,
    store: GroupStorageType,
    ad_setting: bool | None,
    pd_setting: bool,
    expected_missing: tuple[type[Exception], str] | str,
    expected_no_missing: tuple[type[Exception], str] | str,
    missing: bool,
) -> None:
    expected = expected_missing if missing else expected_no_missing
    with (
        ad.settings.override(allow_write_nullable_strings=ad_setting),
        pd.option_context("future.infer_string", pd_setting),
        (
            nullcontext()
            if isinstance(expected, str)
            else pytest.raises(expected[0], match=expected[1])
        ),
    ):
        write_elem(store, "/el", pd.array([pd.NA if missing else "a"], dtype="string"))

    if isinstance(expected, str):
        assert store["el"].attrs["encoding-type"] == expected


def test_categorical_order_type(store):
    # https://github.com/scverse/anndata/issues/853
    cat = pd.Categorical([0, 1], ordered=True)
    write_elem(store, "ordered", cat)
    write_elem(store, "unordered", cat.set_ordered(False))

    assert isinstance(read_elem(store["ordered"]).ordered, bool)
    assert read_elem(store["ordered"]).ordered is True
    assert isinstance(read_elem(store["unordered"]).ordered, bool)
    assert read_elem(store["unordered"]).ordered is False


def test_override_specification():
    """
    Test that trying to overwrite an existing encoding raises an error.
    """
    from copy import deepcopy

    registry = deepcopy(_REGISTRY)

    with pytest.raises(TypeError):

        @registry.register_write(
            ZarrGroup, ad.AnnData, IOSpec("some new type", "0.1.0")
        )
        def _(store, key, adata):
            pass


@pytest.mark.parametrize(
    "value",
    [
        pytest.param({"a": 1}, id="dict"),
        pytest.param(
            lambda: gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS), id="anndata"
        ),
        pytest.param(sparse.random(5, 3, format="csr", density=0.5), id="csr_matrix"),
        pytest.param(sparse.random(5, 3, format="csc", density=0.5), id="csc_matrix"),
        pytest.param(pd.DataFrame({"a": [1, 2, 3]}), id="dataframe"),
        pytest.param(pd.Categorical(list("aabccedd")), id="categorical"),
        pytest.param(
            pd.Categorical(list("aabccedd"), ordered=True), id="categorical-ordered"
        ),
        pytest.param(
            pd.Categorical([1, 2, 1, 3], ordered=True), id="categorical-numeric"
        ),
        pytest.param(
            pd.arrays.IntegerArray(
                np.ones(5, dtype=int), mask=np.array([True, False, True, False, True])
            ),
            id="nullable-integer",
        ),
        pytest.param(pd.array([1, 2, 3]), id="nullable-integer-no-nulls"),
        pytest.param(
            pd.arrays.BooleanArray(
                np.random.randint(0, 2, size=5, dtype=bool),
                mask=np.random.randint(0, 2, size=5, dtype=bool),
            ),
            id="nullable-boolean",
        ),
        pytest.param(
            pd.array([True, False, True, True]), id="nullable-boolean-no-nulls"
        ),
    ],
)
def test_write_to_root(store: GroupStorageType, value):
    """
    Test that elements which are written as groups can we written to the root group.
    """
    if callable(value):
        value = value()
    write_elem(store, "/", value)
    # See: https://github.com/zarr-developers/zarr-python/issues/2716
    if isinstance(store, ZarrGroup):
        store = zarr.open(store.store)
    result = read_elem(store)

    assert_equal(result, value)


@pytest.mark.parametrize("consolidated", [True, False])
@pytest.mark.zarr_io
def test_read_zarr_from_group(tmp_path, consolidated):
    # https://github.com/scverse/anndata/issues/1056
    pth = tmp_path / "test.zarr"
    adata = gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS)

    z = open_write_group(pth)
    write_elem(z, "table/table", adata)

    if consolidated:
        zarr.consolidate_metadata(z.store)

    read_func = zarr.open_consolidated if consolidated else zarr.open

    z = read_func(pth)
    expected = ad.read_zarr(z["table/table"])
    assert_equal(adata, expected)


def test_dataframe_column_uniqueness(store):
    repeated_cols = pd.DataFrame(np.ones((3, 2)), columns=["a", "a"])

    with pytest.raises(
        ValueError,
        match=r"Found repeated column names: \['a'\]\. Column names must be unique\.",
    ):
        write_elem(store, "repeated_cols", repeated_cols)

    index_shares_col_name = pd.DataFrame(
        {"col_name": [1, 2, 3]}, index=pd.Index([1, 3, 2], name="col_name")
    )

    with pytest.raises(
        ValueError,
        match=r"DataFrame\.index\.name \('col_name'\) is also used by a column whose values are different\.",
    ):
        write_elem(store, "index_shares_col_name", index_shares_col_name)

    index_shared_okay = pd.DataFrame(
        {"col_name": [1, 2, 3]}, index=pd.Index([1, 2, 3], name="col_name")
    )

    write_elem(store, "index_shared_okay", index_shared_okay)
    result = read_elem(store["index_shared_okay"])

    assert_equal(result, index_shared_okay)


@pytest.mark.parametrize(
    "copy_on_write",
    [
        pytest.param(True, id="cow"),
        pytest.param(
            False,
            marks=pytest.mark.skipif(
                PANDAS_3, reason="Can’t disable copy-on-write in pandas 3+."
            ),
            id="nocow",
        ),
    ],
)
def test_io_pd_cow(
    *, exit_stack: ExitStack, store: GroupStorageType, copy_on_write: bool
) -> None:
    """See <https://github.com/zarr-developers/numcodecs/issues/514>."""
    if not PANDAS_3:  # Setting copy_on_write always warns in pandas 3, and does nothing
        exit_stack.enter_context(pd.option_context("mode.copy_on_write", copy_on_write))
    orig = gen_adata((3, 2), **GEN_ADATA_NO_XARRAY_ARGS)
    write_elem(store, "adata", orig)
    from_store = read_elem(store["adata"])
    assert_equal(orig, from_store)


def test_read_sparse_array(
    tmp_path: Path,
    sparse_format: Literal["csr", "csc"],
    diskfmt: Literal["h5ad", "zarr"],
):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a = sparse.random(100, 100, format=sparse_format)
    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")
    ad.io.write_elem(f, "mtx", a)
    ad.settings.use_sparse_array_on_read = True
    mtx = ad.io.read_elem(f["mtx"])
    assert issubclass(type(mtx), CSArray)


@pytest.mark.parametrize(
    ("chunks", "expected_chunks"),
    [((1,), (1,)), ((-1,), (120,)), (None, (25,))],
    ids=["small", "minus_one_uses_full", "none_uses_ondisk_chunking"],
)
@pytest.mark.parametrize(
    "arr", [np.arange(120), np.array(["a"] * 120)], ids=["numeric", "string"]
)
def test_chunking_1d_array(
    store: H5Group | ZarrGroup,
    arr: np.ndarray,
    chunks: tuple[int] | None,
    expected_chunks: tuple[int],
):
    write_elem(store, "foo", arr, dataset_kwargs={"chunks": 25})
    arr = read_elem_lazy(store["foo"], chunks=chunks)
    assert arr.chunksize == expected_chunks


@pytest.mark.parametrize(
    ("chunks", "expected_chunks"),
    [
        ((1, 50), (1, 50)),
        ((10, -1), (10, 50)),
        ((10, None), (10, 50)),
        (None, (25, 25)),
    ],
    ids=[
        "small",
        "minus_one_uses_full",
        "none_on_axis_uses_full",
        "none_uses_ondisk_chunking",
    ],
)
def test_chunking_2d_array(
    store: H5Group | ZarrGroup,
    chunks: tuple[int] | None,
    expected_chunks: tuple[int],
):
    write_elem(
        store,
        "foo",
        np.arange(100 * 50).reshape((100, 50)),
        dataset_kwargs={"chunks": (25, 25)},
    )
    arr = read_elem_lazy(store["foo"], chunks=chunks)
    assert arr.chunksize == expected_chunks


@pytest.mark.parametrize(
    ("shape", "expected_chunks"),
    [((100, 50), (100, 50)), ((1100, 1100), (1000, 1000))],
    ids=["under_default", "over_default"],
)
def test_h5_unchunked(
    shape: tuple[int, int],
    expected_chunks: tuple[int, int],
    tmp_path: Path,
):
    with h5py.File(tmp_path / "foo.h5ad", mode="w") as f:
        write_elem(
            f,
            "foo",
            np.arange(shape[0] * shape[1]).reshape(shape),
        )
        arr = read_elem_lazy(f["foo"])
    assert arr.chunksize == expected_chunks


@pytest.mark.zarr_io
def test_write_auto_sharded(tmp_path: Path):
    path = tmp_path / "check.zarr"
    adata = gen_adata((1000, 100), **GEN_ADATA_NO_XARRAY_ARGS)
    with ad.settings.override(auto_shard_zarr_v3=True, zarr_write_format=3):
        adata.write_zarr(path)

    check_all_sharded(zarr.open(path))


@pytest.mark.zarr_io
def test_write_auto_sharded_against_v2_format():
    with pytest.raises(ValueError, match=r"Cannot shard v2 format data."):  # noqa: PT012, SIM117
        with ad.settings.override(zarr_write_format=2):
            with ad.settings.override(auto_shard_zarr_v3=True):
                pass


@pytest.mark.zarr_io
def test_write_auto_cannot_set_v2_format_after_sharding():
    with pytest.raises(ValueError, match=r"Cannot set `zarr_write_format` to 2"):  # noqa: PT012, SIM117
        with ad.settings.override(zarr_write_format=3):
            with ad.settings.override(auto_shard_zarr_v3=True):
                with ad.settings.override(zarr_write_format=2):
                    pass


@pytest.mark.zarr_io
def test_write_auto_sharded_does_not_override(tmp_path: Path):
    z = open_write_group(tmp_path / "arr.zarr", zarr_format=3)
    X = sparse.random(
        100, 100, density=0.1, format="csr", random_state=np.random.default_rng(42)
    )
    with ad.settings.override(auto_shard_zarr_v3=True, zarr_write_format=3):
        ad.io.write_elem(z, "X_default", X)
        shards_default = z["X_default"]["indices"].shards
        new_shards = shards_default[0] // 2
        new_shards = int(new_shards - new_shards % 2)
        ad.io.write_elem(
            z,
            "X_manually_set",
            X,
            dataset_kwargs={
                "shards": (new_shards,),
                "chunks": (int(new_shards / 2),),
            },
        )

    def visitor(key: str, array: zarr.Array):
        assert array.shards == (new_shards,)

    visititems_zarr(z["X_manually_set"], visitor)
