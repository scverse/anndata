"""Which dataframe columns can be written, how, and how failures are reported.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import copy
import pickle
import re
import warnings
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr

import anndata as ad
from anndata._io.specs.dataframe_columns import (
    UnsupportedColumnTypeError,
    _object_column_to_writable,
    normalize_column,
)
from anndata._io.specs.registry import _REGISTRY, IORegistryError, IOSpec
from anndata._io.specs.validate import iter_unwriteable
from anndata.tests.helpers import (
    adata_with_obs,
    as_dataframe,
    exception_notes,
    obs_column,
    obs_index,
    open_store,
    store_type_of,
    write_and_read,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Mapping
    from pathlib import Path
    from typing import Any, Literal

    from pandas.api.extensions import ExtensionArray
    from zarr.storage import MemoryStore

    from anndata._io.specs.registry import Writer


# families with a writer and their on-disk dtype; factories give fresh input
SUPPORTED: dict[str, tuple[Callable[[], pd.Series], str]] = {
    "int64": (lambda: obs_column([1, 2, 3]), "int64"),
    "float64-nan": (lambda: obs_column([1.5, np.nan, 2.0]), "float64"),
    "bool": (lambda: obs_column([True, False, True]), "bool"),
    "Int64": (lambda: obs_column([1, None, 3], "Int64"), "Int64"),
    "boolean": (lambda: obs_column([True, None, False], "boolean"), "boolean"),
    "string": (lambda: obs_column(["a", None, "a"], "string"), "string"),
    "category": (lambda: obs_column(["a", "b", None]).astype("category"), "category"),
    "category-ordered": (
        lambda: pd.Series(
            pd.Categorical(["lo", "hi", "lo"], categories=["lo", "hi"], ordered=True),
            index=obs_index(),
        ),
        "category",
    ),
}

# families without a writer: error type and the hint in its notes
UNSUPPORTED: dict[str, tuple[Callable[[], pd.Series], type[Exception], str]] = {
    "Float64": (
        lambda: obs_column([1.5, None, 2.0], "Float64"),
        IORegistryError,
        r"no on-disk encoding for nullable floats.*na_value=np\.nan",
    ),
    "Float32": (
        lambda: obs_column([1.5, None, 2.0], "Float32"),
        IORegistryError,
        r"no on-disk encoding for nullable floats",
    ),
    "datetime64": (
        lambda: obs_column(pd.to_datetime(["2020-01-01", "NaT", "2021-01-01"])),
        IORegistryError,
        r"datetime/timedelta columns are not supported yet, see scverse/anndata#2238",
    ),
    "datetime64-tz": (
        lambda: obs_column(
            pd.to_datetime(["2020-01-01", "NaT", "2021-01-01"]).tz_localize("UTC")
        ),
        IORegistryError,
        r"scverse/anndata#2238",
    ),
    "timedelta64": (
        lambda: obs_column(pd.to_timedelta(["1h", None, "2h"])),
        IORegistryError,
        r"scverse/anndata#2238",
    ),
    "int64-pyarrow": (
        lambda: obs_column([1, None, 3], "int64[pyarrow]"),
        IORegistryError,
        r"pyarrow-backed columns are not supported.*numpy_nullable",
    ),
    "period": (
        lambda: obs_column(pd.period_range("2020-01", periods=3, freq="M")),
        IORegistryError,
        r"Convert the column to a numpy, nullable integer/boolean, string, or categorical",
    ),
    "interval": (
        lambda: obs_column(pd.interval_range(0, 3)),
        IORegistryError,
        r"Convert the column to a numpy",
    ),
    "sparse": (
        lambda: obs_column(pd.arrays.SparseArray([0.0, 1.0, 0.0])),
        IORegistryError,
        r"Convert the column to a numpy",
    ),
}


@pytest.mark.parametrize("family", list(SUPPORTED))
def test_supported_column_roundtrip(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    family: str,
) -> None:
    """Preflight and write agree, and the column comes back with the expected dtype."""
    make, expected_dtype = SUPPORTED[family]
    col = make()
    expected_na = col.isna().tolist()
    expected_values = col.dropna().tolist()
    adata = adata_with_obs(col=col)
    assert list(iter_unwriteable(adata)) == []
    assert not adata.unwriteable()

    with ad.settings.override(
        validate_before_write=True, allow_write_nullable_strings=True
    ):
        read = as_dataframe(write_and_read(adata, diskfmt, diskfmt_store).obs)["col"]
    assert read.dtype == expected_dtype
    assert read.isna().tolist() == expected_na
    assert read.dropna().tolist() == expected_values


@pytest.mark.parametrize("family", list(UNSUPPORTED))
def test_unsupported_column_errors(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    family: str,
) -> None:
    """Preflight and write report the same error, naming the column and a way out."""
    make, exc_type, hint = UNSUPPORTED[family]
    store_type = store_type_of(diskfmt)

    adata = adata_with_obs(col=make())
    errors = list(iter_unwriteable(adata, store_type=store_type))
    assert [path for path, _ in errors] == ["/obs/col"]
    assert isinstance(errors[0][1], exc_type)
    assert re.search(hint, exception_notes(errors[0][1]) + str(errors[0][1]), re.DOTALL)
    assert adata.unwriteable(store_type=store_type)

    # the writer raises the same type
    adata = adata_with_obs(col=make())
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(exc_type) as exc_info,
    ):
        write_and_read(adata, diskfmt, diskfmt_store)
    assert re.search(
        hint, exception_notes(exc_info.value) + str(exc_info.value), re.DOTALL
    )
    assert re.search(r"'col'", str(exc_info.value) + exception_notes(exc_info.value))

    # one problem is raised as is, not as a group
    adata = adata_with_obs(col=make())
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(exc_type) as exc_info,
    ):
        write_and_read(adata, diskfmt, diskfmt_store)
    assert "Element: /obs/col" in exception_notes(exc_info.value)


@pytest.mark.parametrize(
    ("values", "expected_dtype"),
    [
        pytest.param(np.array([1, 2, None, 3]), "Int64", id="int-none"),
        pytest.param([1, 2, pd.NA, 3], "Int64", id="int-NA"),
        pytest.param([np.int32(1), 2, None, np.int64(3)], "Int64", id="npint-none"),
        pytest.param([True, False, None, False], "boolean", id="bool-none"),
        pytest.param([True, False, pd.NA, False], "boolean", id="bool-NA"),
        pytest.param([True, np.nan, False, np.True_], "boolean", id="bool-nan"),
        pytest.param([1.5, 2.0, None, 3.0], "float64", id="float-none"),
        pytest.param([1, 2.5, None, 3], "float64", id="int-float-none"),
        pytest.param([1, 2, 3, 4], "int64", id="int-no-missing"),
        pytest.param([True, False, True, False], "bool", id="bool-no-missing"),
        pytest.param([1.5, 2.0, 2.5, 3.0], "float64", id="float-no-missing"),
    ],
)
def test_write_object_column_converted(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: list[object] | np.ndarray,
    expected_dtype: str,
) -> None:
    """Object columns of ints/bools/floats (with or without missing values) are writable.

    See https://github.com/scverse/anndata/issues/1068
    """
    # `dtype=object` keeps pandas from coercing e.g. `[1.5, None]` to float64 up front;
    # `test_write_object_column_realistic` covers the idioms that produce such columns
    col = pd.Series(np.array(values, dtype=object), index=[f"c{i}" for i in range(4)])
    assert col.dtype == object
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}), var=pd.DataFrame(index=col.index))
    # `varm` dataframes go through the same code path
    orig.varm["df"] = pd.DataFrame({"col": col.to_numpy()}, index=orig.var_names)
    assert list(iter_unwriteable(orig)) == []

    with ad.settings.override(validate_before_write=True):
        getattr(orig, f"write_{diskfmt}")(diskfmt_store)
    # The in-memory object is not modified by writing
    assert orig.obs["col"].dtype == object
    curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)

    assert isinstance(curr.obs, pd.DataFrame)
    varm_df = curr.varm["df"]
    assert isinstance(varm_df, pd.DataFrame)
    for read_col in (curr.obs["col"], varm_df["col"]):
        assert read_col.dtype == expected_dtype
        assert read_col.isna().tolist() == col.isna().tolist()
        assert read_col.dropna().tolist() == col.dropna().tolist()


@pytest.mark.parametrize(
    ("values", "match"),
    [
        pytest.param(
            [2**63, 1, None], r"inferred as 'integer'.*'Int64' failed", id="int-too-big"
        ),
        pytest.param([-(2**63) - 1, None], r"'Int64' failed", id="int-too-small"),
        # without missing values, the direct numpy cast must fail the same way
        pytest.param(
            [2**63, 1], r"inferred as 'integer'.*'Int64' failed", id="int-too-big-full"
        ),
        pytest.param(
            [2**53 + 1, 0.5, None],
            r"inferred as 'mixed-integer-float'.*exactly as float64",
            id="int-float-lossy",
        ),
        pytest.param(
            [2**53 + 1, 0.5], r"exactly as float64", id="int-float-lossy-full"
        ),
        pytest.param(
            [-(2**60) - 1, 0.5, None], r"exactly as float64", id="negint-float-lossy"
        ),
        pytest.param(
            [np.int64(2**53 + 1), 0.5, None],
            r"exactly as float64",
            id="npint-float-lossy",
        ),
        pytest.param(
            [np.uint64(2**64 - 1), 0.5, None],
            r"exactly as float64",
            id="npuint-float-lossy",
        ),
        pytest.param([1, "a", None, 2], r"inferred as 'mixed-integer'", id="mixed"),
        pytest.param([None, None], r"inferred as 'empty'", id="all-missing"),
        pytest.param(["a", None], r"inferred as 'string'", id="str-missing"),
        pytest.param([b"a", None], r"inferred as 'bytes'", id="bytes-missing"),
    ],
)
def test_write_object_column_unconvertible(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: list[object],
    match: str,
) -> None:
    """Object columns that cannot be written raise an error naming the column.

    This also covers zarr, which used to silently write `str(value)` for such columns.
    """
    col = pd.Series(values, index=[f"c{i}" for i in range(len(values))], dtype=object)
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}))
    assert orig.obs["col"].dtype == object
    errors = list(iter_unwriteable(orig, store_type=store_type_of(diskfmt)))
    assert [(p, type(e)) for p, e in errors] == [
        ("/obs/col", UnsupportedColumnTypeError)
    ]
    with pytest.raises(TypeError, match=rf"column 'col'.*{match}"):
        getattr(orig, f"write_{diskfmt}")(
            diskfmt_store, convert_strings_to_categoricals=False
        )


@pytest.mark.parametrize(
    "values",
    [
        pytest.param([2**63 - 1, -(2**63), None], id="int-bounds"),
        pytest.param([2**63 - 1, -(2**63)], id="int-bounds-full"),
        pytest.param([True, False], id="bool-full"),
        pytest.param([2**53, 0.5, None], id="int-float-exact"),
        pytest.param([2**60, -(2**60), 0.5, None], id="bigint-float-exact"),
        pytest.param([1e300, -(2.0**60), 0.5, None], id="bigfloat-exact"),
        pytest.param([float("inf"), -float("inf"), 0.5, None], id="inf-float"),
        pytest.param([np.int64(2**53), np.float32(0.5), None], id="npint-float-exact"),
        pytest.param(["a", "b", "c"], id="str"),
        pytest.param([b"a", b"b", b"c"], id="bytes"),
        pytest.param(["a", b"b", "c"], id="str-bytes-mixed"),
        pytest.param([], id="empty"),
    ],
)
def test_write_object_column_roundtrip_values(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: list[object],
) -> None:
    """Values at the edge of what is convertible and object string columns round-trip."""
    col = pd.Series(values, index=[f"c{i}" for i in range(len(values))], dtype=object)
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}))
    assert orig.obs["col"].dtype == object
    getattr(orig, f"write_{diskfmt}")(
        diskfmt_store, convert_strings_to_categoricals=False
    )
    curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert isinstance(curr.obs, pd.DataFrame)
    read = curr.obs["col"]
    assert read.isna().tolist() == col.isna().tolist()
    expected = [v.decode() if isinstance(v, bytes) else v for v in col.dropna()]
    assert [v.item() if isinstance(v, np.generic) else v for v in expected] == [
        v.item() if isinstance(v, np.generic) else v for v in read.dropna()
    ]


def _realistic_object_columns() -> dict[str, tuple[pd.Series, str]]:
    """Object columns the way users actually produce them, with the dtype expected on disk."""
    rng = np.random.default_rng(0)
    counts = pd.Series(
        rng.integers(100, 5000, size=8), index=[f"c{i}" for i in range(8)]
    )
    flags = counts > 2000
    ints_nullable = pd.array([1, None, 3, 4, None, 6, 7, 8], dtype="Int64")
    bools_nullable = pd.array(
        [True, None, False, True, None, True, False, True], dtype="boolean"
    )
    partial = pd.DataFrame({"rank": counts.iloc[:4].astype(object)})
    return {
        "np.where-int-none": (
            pd.Series(
                np.where(flags, counts, None),  # type: ignore[call-overload]  # stubs: no None
                index=counts.index,
            ),
            "Int64",
        ),
        "np.where-bool-none": (
            pd.Series(
                np.where(flags, flags, None),  # type: ignore[call-overload]  # stubs: no None
                index=counts.index,
            ),
            "boolean",
        ),
        "Int64-astype-object": (
            pd.Series(ints_nullable, index=counts.index).astype(object),
            "Int64",
        ),
        "boolean-astype-object": (
            pd.Series(bools_nullable, index=counts.index).astype(object),
            "boolean",
        ),
        # a join that leaves rows unmatched fills them with NaN, keeping dtype object
        "join-unmatched-rows": (counts.to_frame().join(partial)["rank"], "Int64"),
        "float-where-astype-object": (
            counts.astype(float).where(flags).astype(object),
            "float64",
        ),
        # scirpy's case (scverse/scirpy#434): a plain list of bools and None
        "list-bool-none": (
            pd.Series(
                [True, None, False, True, None, True, False, True], index=counts.index
            ),
            "boolean",
        ),
    }


@pytest.mark.parametrize("case", list(_realistic_object_columns()))
def test_write_object_column_realistic(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    case: str,
) -> None:
    """Object columns produced by common pandas/numpy idioms round-trip."""
    col, expected_dtype = _realistic_object_columns()[case]
    assert col.dtype == object
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}))
    getattr(orig, f"write_{diskfmt}")(diskfmt_store)
    curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert isinstance(curr.obs, pd.DataFrame)
    read = curr.obs["col"]
    assert read.dtype == expected_dtype
    assert read.isna().tolist() == col.isna().tolist()
    assert read.dropna().tolist() == col.dropna().tolist()


def test_write_object_column_elem_and_empty_view(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Object columns are converted for direct `write_elem` calls and zero-row views."""
    col = pd.Series([1, None, 3], index=list("abc"), dtype=object)
    df = pd.DataFrame({"col": col, "s": pd.Series(list("xyz"), index=list("abc"))})
    store = open_store(diskfmt_store)
    ad.io.write_elem(store, "df", df)
    read = ad.io.read_elem(store["df"])
    assert isinstance(read, pd.DataFrame)
    if isinstance(store, h5py.File):
        store.close()  # the h5ad write below reopens the file
    assert read["col"].dtype == "Int64"
    assert read["col"].tolist()[::2] == [1, 3]
    assert df["col"].dtype == object

    adata = ad.AnnData(obs=df)
    empty_view = adata[:0]
    assert empty_view.obs["col"].dtype == object
    getattr(empty_view, f"write_{diskfmt}")(diskfmt_store)
    curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert curr.n_obs == 0
    assert "col" in curr.obs


class _LongDouble(np.longdouble):
    pass


class _NaNWithInt(float):
    """A NaN that numpy casts to the integer 2."""

    def __int__(self) -> int:
        return 2


@pytest.mark.skipif(
    np.finfo(np.longdouble).nmant <= np.finfo(np.float64).nmant,
    reason="np.longdouble is not wider than float64 on this platform",
)
@pytest.mark.parametrize("with_missing", [True, False], ids=["missing", "full"])
@pytest.mark.parametrize("sign", [1, -1], ids=["pos", "neg"])
@pytest.mark.parametrize("subclass", [False, True], ids=["longdouble", "subclass"])
def test_write_object_column_longdouble_lossy(
    *, with_missing: bool, sign: int, subclass: bool
) -> None:
    """Extended-precision scalars below 2**53 are still checked for exactness."""
    cls = _LongDouble if subclass else np.longdouble
    x = cls(sign * (np.longdouble(1) + np.finfo(np.longdouble).eps))
    tail = [None] if with_missing else []
    values = np.array([x, 0.5, *tail], dtype=object)
    with pytest.raises(TypeError, match=r"exactly as float64"):
        _object_column_to_writable(values, path="col")
    exact = np.array([cls(sign * 0.5), 2.0, *tail], dtype=object)
    np.testing.assert_array_equal(
        _object_column_to_writable(exact, path="col"),
        np.array([sign * 0.5, 2.0, *([np.nan] if with_missing else [])]),
    )


@pytest.mark.parametrize(
    "missing",
    [
        pytest.param(np.complex128(complex(2, np.nan)), id="complex-nan"),
        pytest.param(_NaNWithInt(np.nan), id="float-nan-subclass"),
    ],
)
def test_write_object_column_int_keeps_exotic_missing(missing: object) -> None:
    """Missing markers numpy would happily cast to an int stay missing."""
    values = np.array([1, missing, 3], dtype=object)
    assert pd.api.types.infer_dtype(values, skipna=True) == "integer"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        converted = _object_column_to_writable(values, path="col")
    assert isinstance(converted, pd.arrays.IntegerArray)
    assert converted.isna().tolist() == [False, True, False]
    assert converted[[0, 2]].tolist() == [1, 3]


def test_write_object_column_nan_distinct_from_na() -> None:
    """With pandas distinguishing NaN from NA, an object float NaN is written as NaN."""
    if not hasattr(pd.options.future, "distinguish_nan_and_na"):
        pytest.skip("pandas has no `future.distinguish_nan_and_na` option")
    values = np.array([0.5, np.nan, None], dtype=object)
    with pd.option_context("future.distinguish_nan_and_na", True):  # noqa: FBT003
        np.testing.assert_array_equal(
            _object_column_to_writable(values, path="col"),
            np.array([0.5, np.nan, np.nan]),
        )


@pytest.mark.parametrize("dtype", ["Float64", "Float32"])
@pytest.mark.parametrize("distinguish", [False, True], ids=["collapse", "distinguish"])
def test_nullable_float_passthrough(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    dtype: str,
    distinguish: bool,
) -> None:
    """Nullable floats are never converted: a masked NA and an unmasked NaN stay distinct."""
    if not hasattr(pd.options.future, "distinguish_nan_and_na"):
        pytest.skip("pandas has no `future.distinguish_nan_and_na` option")
    store_type = store_type_of(diskfmt)
    with pd.option_context("future.distinguish_nan_and_na", distinguish):
        values = pd.arrays.FloatingArray(
            np.array([1.5, 0.0, np.nan], dtype=dtype.lower()),
            np.array([False, True, False]),
        )
        before = values.copy()
        assert normalize_column(values, path="col") is values
        pd.testing.assert_extension_array_equal(values, before)
        assert values.isna().tolist() == [False, True, False]

        adata = adata_with_obs(col=pd.Series(values, index=obs_index()))
        errors = list(iter_unwriteable(adata, store_type=store_type))
        assert [path for path, _ in errors] == ["/obs/col"]
        assert isinstance(errors[0][1], IORegistryError)
        assert "no on-disk encoding for nullable floats" in exception_notes(
            errors[0][1]
        )
        pd.testing.assert_extension_array_equal(
            as_dataframe(adata.obs)["col"].array, before
        )
        with (
            ad.settings.override(validate_before_write=False),
            pytest.raises(IORegistryError) as exc_info,
        ):
            write_and_read(adata, diskfmt, diskfmt_store)
        assert "na_value=np.nan" in exception_notes(exc_info.value)


@pytest.mark.parametrize(
    "values",
    [
        pytest.param(np.array([1, None, 3], dtype=object), id="object-int"),
        pytest.param(np.array([1.5, None, 2.0], dtype=object), id="object-float"),
        pytest.param(np.array(["a", "b"], dtype=object), id="object-str"),
        pytest.param(pd.array([1, None], dtype="Int64"), id="Int64"),
        pytest.param(pd.array(["a", None], dtype="string"), id="string"),
        pytest.param(pd.Categorical(["a", "b", "a"]), id="category"),
        pytest.param(np.array([1.5, np.nan]), id="float64"),
    ],
)
def test_normalize_column_idempotent_and_non_mutating(
    values: np.ndarray | ExtensionArray,
) -> None:
    original = pd.Series(values, copy=True)
    once = normalize_column(values, path="col")
    # already writable: returned as is
    assert normalize_column(once, path="col") is once
    pd.testing.assert_series_equal(pd.Series(values), original)


@pytest.mark.parametrize("copier", ["pickle", "copy"])
def test_unsupported_object_column_error_roundtrips(copier: str) -> None:
    """Errors from worker processes arrive pickled, so they must survive it."""
    with pytest.raises(UnsupportedColumnTypeError) as exc_info:
        normalize_column(np.array([1, "a", None], dtype=object), path="/obs/col")
    e = exc_info.value
    e.add_note("Element: /obs/col")
    restored = pickle.loads(pickle.dumps(e)) if copier == "pickle" else copy.copy(e)
    assert type(restored) is UnsupportedColumnTypeError
    assert str(restored) == str(e)
    assert str(e).startswith(
        "Cannot write column '/obs/col' (dtype object, inferred as 'mixed-integer'): "
    )
    assert exception_notes(restored) == "Element: /obs/col"
    assert (restored.path, restored.kind) == ("/obs/col", "mixed-integer")


@pytest.mark.parametrize("storage", ["python", "pyarrow"])
@pytest.mark.parametrize(
    "allow", [False, None, True], ids=["allow=F", "allow=N", "allow=T"]
)
@pytest.mark.parametrize("categorize", [True, False], ids=["s2c", "no-s2c"])
def test_nullable_string_columns(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    storage: str,
    allow: bool | None,
    categorize: bool,
) -> None:
    """Nullable string columns with missing values write and read back.

    With `allow_write_nullable_strings=False` the writer refuses missing values,
    unless the column was converted to categorical first. With `None`, the writer
    follows `pd.options.future.infer_string` (always on in pandas >= 3).
    """
    col = obs_column(["a", None, "a"], f"string[{storage}]")
    adata = adata_with_obs(col=col)
    with ad.settings.override(
        validate_before_write=True, allow_write_nullable_strings=allow
    ):
        if allow is False and not categorize:
            with pytest.raises(ValueError, match=r"with missing values"):
                getattr(adata, f"write_{diskfmt}")(
                    diskfmt_store, convert_strings_to_categoricals=categorize
                )
            return
        if allow is None and not categorize and not pd.options.future.infer_string:
            # pandas < 3 without `infer_string`: the writer wants a choice
            with pytest.raises(RuntimeError, match=r"allow_write_nullable_strings"):
                getattr(adata, f"write_{diskfmt}")(
                    diskfmt_store, convert_strings_to_categoricals=categorize
                )
            return
        getattr(adata, f"write_{diskfmt}")(
            diskfmt_store, convert_strings_to_categoricals=categorize
        )
    read = as_dataframe(getattr(ad, f"read_{diskfmt}")(diskfmt_store).obs)["col"]
    assert read.dtype == ("category" if categorize else "string")
    assert read.isna().tolist() == [False, True, False]
    assert read.dropna().tolist() == ["a", "a"]


@pytest.mark.parametrize("categorize", [True, False], ids=["s2c", "no-s2c"])
def test_categorical_column_set_unchanged(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    categorize: bool,
) -> None:
    """Column normalisation does not change which columns become categorical.

    See https://github.com/scverse/anndata/issues/2610
    """
    adata = adata_with_obs(
        repeated=obs_column(["a", "a", "b"]),
        unique=obs_column(["a", "b", "c"]),
        obj=obs_column(np.array([1, None, 3], dtype=object)),
        cat=obs_column(["x", "y", "x"]).astype("category"),
    )
    expected = {"repeated", "cat"} if categorize else {"cat"}
    getattr(adata, f"write_{diskfmt}")(
        diskfmt_store, convert_strings_to_categoricals=categorize
    )
    read = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    for obs in (adata.obs, read.obs):
        categorical = {c for c in obs if isinstance(obs[c].dtype, pd.CategoricalDtype)}
        assert categorical == expected
    assert read.obs["obj"].dtype == "Int64"


def _write_period(
    f: h5py.Group | zarr.Group,
    k: str,
    v: pd.arrays.PeriodArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any],
) -> None:
    _writer.write_elem(f, k, v.astype(np.int64), dataset_kwargs=dataset_kwargs)


@pytest.fixture
def period_writer() -> Generator[None]:
    """Temporarily register a writer for `PeriodArray`."""
    keys: list[tuple[type, type, frozenset[str]]] = [
        (store_type, pd.arrays.PeriodArray, frozenset())
        for store_type in (h5py.Group, zarr.Group)
    ]
    assert not any(key in _REGISTRY.write for key in keys)

    for store_type, src_type, _ in keys:
        _REGISTRY.register_write(store_type, src_type, IOSpec("period-test", "0.1.0"))(
            _write_period
        )
    try:
        yield
    finally:
        for key in keys:
            _REGISTRY.write.pop(key)
        _REGISTRY.write_specs.pop(pd.arrays.PeriodArray)


@pytest.mark.usefixtures("period_writer")
def test_custom_registered_writer_is_honoured(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Unknown extension arrays are passed to the registry, not rejected up front."""
    adata = adata_with_obs(
        period=obs_column(pd.period_range("2020-01", periods=3, freq="M"))
    )
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        getattr(adata, f"write_{diskfmt}")(
            diskfmt_store, convert_strings_to_categoricals=False
        )


CATEGORY_FAMILIES: dict[str, tuple[Callable[[], pd.Categorical[Any]], bool]] = {
    # name: (factory, accepted by the preflight)
    "str": (lambda: pd.Categorical(["a", "b", "a"]), True),
    "Int64": (lambda: pd.Categorical(pd.array([1, None, 2], dtype="Int64")), True),
    "Float64": (
        lambda: pd.Categorical(pd.array([1.5, 2.0, 1.5], dtype="Float64")),
        True,
    ),
    "arrow-int": (
        lambda: pd.Categorical(pd.array([1, 2, 1], dtype="int64[pyarrow]")),
        True,
    ),
    "datetime": (
        lambda: pd.Categorical(
            pd.to_datetime(["2020-01-01", "2021-01-01", "2020-01-01"])
        ),
        True,
    ),
    "timedelta": (lambda: pd.Categorical(pd.to_timedelta(["1h", "2h", "1h"])), True),
    "interval": (
        lambda: pd.Categorical.from_codes(
            [0, 1, 0],
            dtype=pd.CategoricalDtype(pd.IntervalIndex.from_tuples([(0, 1), (1, 2)])),
        ),
        False,
    ),
    "period": (
        lambda: pd.Categorical(pd.period_range("2020-01", periods=3, freq="M")),
        False,
    ),
    "arrow-str": (
        lambda: pd.Categorical(
            pd.Series(
                ["a", "b", "a"],
                dtype=pd.ArrowDtype(pytest.importorskip("pyarrow").string()),
            ).array
        ),
        True,
    ),
    "datetime-tz": (
        lambda: pd.Categorical(
            pd.to_datetime(["2020-01-01", "2021-01-01", "2020-01-01"]).tz_localize(
                "UTC"
            )
        ),
        False,
    ),
}


@pytest.mark.parametrize("family", list(CATEGORY_FAMILIES))
def test_categorical_categories_agreement(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    family: str,
) -> None:
    """The walker judges categories on what the writer stores (`to_numpy()`).

    Object arrays of non-strings (intervals, periods, tz-aware datetimes) are refused
    on both formats; the writer itself is unchanged: h5py fails, zarr stores strings.
    """
    make, accepted = CATEGORY_FAMILIES[family]
    store_type = store_type_of(diskfmt)
    adata = adata_with_obs(cat=pd.Series(make(), index=obs_index()))
    errors = list(iter_unwriteable(adata, store_type=store_type))
    assert [p for p, _ in errors] == ([] if accepted else ["/obs/cat/categories"])
    assert adata.unwriteable(store_type=store_type) is not accepted
    if not accepted:
        assert isinstance(errors[0][1], TypeError)
        assert "Cannot write categories" in str(errors[0][1])
        with ad.settings.override(validate_before_write=True), pytest.raises(TypeError):
            write_and_read(adata, diskfmt, diskfmt_store)
    if diskfmt == "h5ad" and (not accepted or family in {"datetime", "timedelta"}):
        # h5py refuses object and datetime64/timedelta64 data
        with pytest.raises(TypeError, match=r"convert|conversion path"):
            write_and_read(adata, diskfmt, diskfmt_store)
        return
    read = as_dataframe(write_and_read(adata, diskfmt, diskfmt_store).obs)["cat"]
    assert read.dtype == "category"
    assert read.isna().tolist() == adata.obs["cat"].isna().tolist()
    if not accepted:  # zarr stores what it cannot encode as strings
        assert pd.api.types.is_string_dtype(read.cat.categories.dtype)
