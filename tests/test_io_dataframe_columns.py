"""Which dataframe columns can be written, how, and how failures are reported.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import contextlib
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
from anndata.acc import A
from anndata.tests.helpers import open_store

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping
    from pathlib import Path
    from typing import Any, Literal

    from pandas.api.extensions import ExtensionArray
    from zarr.storage import MemoryStore

    from anndata._io.specs.registry import Writer


def _adata_with_obs(**columns: object) -> ad.AnnData:
    """An :class:`~anndata.AnnData` without `X` whose three-row `obs` holds `columns`."""
    return ad.AnnData(obs=pd.DataFrame(columns, index=[f"c{i}" for i in range(3)]))


@pytest.mark.parametrize(
    ("values", "hint"),
    [
        pytest.param(
            pd.array([1.5, None, 2.0], dtype="Float64"),
            r"no on-disk encoding for nullable floats.*na_value=np\.nan",
            id="Float64",
        ),
        pytest.param(
            pd.to_datetime(["2020-01-01", "NaT", "2021-01-01"]),
            r"datetime/timedelta columns are not supported yet, see scverse/anndata#2238",
            id="datetime64",
        ),
        pytest.param(
            pd.array([1, None, 3], dtype="int64[pyarrow]"),
            r"pyarrow-backed columns are not supported.*numpy_nullable",
            id="int64-pyarrow",
        ),
        pytest.param(
            pd.period_range("2020-01", periods=3, freq="M"),
            r"Convert the column to a numpy, nullable integer/boolean, string, or categorical",
            id="period",
        ),
    ],
)
def test_unsupported_column_note(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: object,
    hint: str,
) -> None:
    """Columns without a writer raise the registry error, noting the column and a way out."""
    adata = _adata_with_obs(col=values)
    with pytest.raises(IORegistryError) as exc_info:
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    notes = "\n".join(exc_info.value.__notes__)
    assert "Column 'col' (dtype " in notes
    assert re.search(hint, notes, re.DOTALL)


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
        pytest.param([2**63 - 1, -(2**63), None], "Int64", id="int-bounds"),
        pytest.param([2**63 - 1, -(2**63)], "int64", id="int-bounds-full"),
        pytest.param([np.inf, -np.inf, 0.5, None], "float64", id="inf-float"),
        pytest.param(["a", "b", "c"], None, id="str"),
        pytest.param([b"a", b"b", b"c"], None, id="bytes"),
        pytest.param(["a", b"b", "c"], None, id="str-bytes-mixed"),
        pytest.param([], None, id="empty"),
    ],
)
def test_write_object_column_roundtrip(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: list[object] | np.ndarray,
    expected_dtype: str | None,
) -> None:
    """Convertible object columns and object string columns round-trip."""
    col = pd.Series(values, index=[f"c{i}" for i in range(len(values))], dtype=object)
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}))
    assert not orig.unwriteable()
    getattr(orig, f"write_{diskfmt}")(
        diskfmt_store, convert_strings_to_categoricals=False
    )
    assert orig.obs["col"].dtype == object
    curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert isinstance(curr.obs, pd.DataFrame)
    read = curr.obs["col"]
    if expected_dtype is not None:
        assert read.dtype == expected_dtype
    assert read.isna().tolist() == col.isna().tolist()
    expected = [v.decode() if isinstance(v, bytes) else v for v in col.dropna()]
    assert read.dropna().tolist() == expected


@pytest.mark.parametrize(
    ("values", "match"),
    [
        pytest.param(
            [2**63, 1, None], r"inferred as 'integer'.*'Int64' failed", id="int-too-big"
        ),
        pytest.param([-(2**63) - 1, None], r"'Int64' failed", id="int-too-small"),
        pytest.param(
            [2**63, 1], r"inferred as 'integer'.*'Int64' failed", id="int-too-big-full"
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
    """Object columns that cannot be written raise an error naming the column."""
    col = pd.Series(values, index=[f"c{i}" for i in range(len(values))], dtype=object)
    orig = ad.AnnData(obs=pd.DataFrame({"col": col}))
    assert orig.obs["col"].dtype == object
    with pytest.raises(TypeError, match=rf"column 'col'.*{match}"):
        getattr(orig, f"write_{diskfmt}")(
            diskfmt_store, convert_strings_to_categoricals=False
        )


def _realistic_object_columns() -> dict[str, tuple[pd.Series, str]]:
    """Object columns the way users actually produce them, with the dtype expected on disk."""
    rng = np.random.default_rng(0)
    counts = pd.Series(
        rng.integers(100, 5000, size=8), index=[f"c{i}" for i in range(8)]
    )
    flags = counts > 2000
    partial = pd.DataFrame({"rank": counts.iloc[:4].astype(object)})
    return {
        "np.where-int-none": (
            pd.Series(
                np.where(flags, counts, None),  # type: ignore[call-overload]
                index=counts.index,
            ),
            "Int64",
        ),
        "np.where-bool-none": (
            pd.Series(
                np.where(flags, flags, None),  # type: ignore[call-overload]
                index=counts.index,
            ),
            "boolean",
        ),
        "join-unmatched-rows": (counts.to_frame().join(partial)["rank"], "Int64"),
        "float-where-astype-object": (
            counts.astype(float).where(flags).astype(object),
            "float64",
        ),
        "list-bool-none": (
            pd.Series(
                [True, None, False, True, None, True, False, True], index=counts.index
            ),
            "boolean",
        ),
    }


@pytest.mark.parametrize("distinguish", [False, True], ids=["collapse", "distinguish"])
@pytest.mark.parametrize("case", list(_realistic_object_columns()))
def test_write_object_column_realistic(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    case: str,
    distinguish: bool,
) -> None:
    """Object columns produced by common pandas/numpy idioms round-trip.

    pandas >= 3 can keep NaN apart from NA, so an integer column with NaN
    from a join must still be written as `Int64`.
    """
    has_option = hasattr(pd.options.future, "distinguish_nan_and_na")
    if distinguish and not has_option:
        pytest.skip("pandas has no `future.distinguish_nan_and_na` option")
    with (
        pd.option_context("future.distinguish_nan_and_na", distinguish)
        if has_option
        else contextlib.nullcontext()
    ):
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
        store.close()
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


class _NaNWithInt(float):
    """A NaN that numpy casts to the integer 2."""

    def __int__(self) -> int:
        return 2


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
        converted = _object_column_to_writable(values, name="col")
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
            _object_column_to_writable(values, name="col"),
            np.array([0.5, np.nan, np.nan]),
        )


def test_write_object_column_large_int_float_rounds() -> None:
    """Integers above 2**53 mixed with floats are rounded to float64, not refused."""
    values = np.array([2**53 + 1, 0.5, None], dtype=object)
    np.testing.assert_array_equal(
        _object_column_to_writable(values, name="col"), np.array([2.0**53, 0.5, np.nan])
    )


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
    once = normalize_column(values, name="col")
    assert normalize_column(once, name="col") is once
    pd.testing.assert_series_equal(pd.Series(values), original)


@pytest.mark.parametrize("categorize", [True, False], ids=["s2c", "no-s2c"])
def test_categorical_column_set_unchanged(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    categorize: bool,
) -> None:
    """Column normalisation does not change which columns become categorical."""
    adata = _adata_with_obs(
        repeated=["a", "a", "b"],
        unique=["a", "b", "c"],
        obj=np.array([1, None, 3], dtype=object),
        cat=pd.Categorical(["x", "y", "x"]),
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
def period_writer(diskfmt: Literal["h5ad", "zarr"]) -> Generator[None]:
    """Temporarily register a writer for `PeriodArray` with the store type under test."""
    store_type = h5py.Group if diskfmt == "h5ad" else zarr.Group
    key = (store_type, pd.arrays.PeriodArray, frozenset[str]())
    assert key not in _REGISTRY.write
    spec = IOSpec("period-test", "0.1.0")
    _REGISTRY.register_write(store_type, pd.arrays.PeriodArray, spec)(_write_period)
    try:
        yield
    finally:
        _REGISTRY.write.pop(key)
        _REGISTRY.write_specs.pop(pd.arrays.PeriodArray)


@pytest.mark.usefixtures("period_writer")
def test_custom_registered_writer_is_honoured(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Unknown extension arrays are passed to the registry, not rejected up front."""
    adata = _adata_with_obs(period=pd.period_range("2020-01", periods=3, freq="M"))
    assert not adata.unwriteable()
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)


def test_use_object_column_converted_like_full_write(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """`use=A.obs[col]` writes an object column the way a full write does."""
    adata = _adata_with_obs(col=[1, 2, 3])
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    adata.obs["o"] = np.array([1, None, 3], dtype=object)
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, use=A.obs["o"])
    read: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert isinstance(read.obs, pd.DataFrame)
    assert read.obs["o"].dtype == "Int64"
    assert read.obs["o"].isna().tolist() == [False, True, False]
    assert adata.obs["o"].dtype == object


def test_use_unconvertible_object_column_writes_nothing(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """All parts of `use=` are converted before the first one is written."""
    adata = _adata_with_obs(col=[1, 2, 3])
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    adata.obsm["good"] = np.ones((3, 2))
    adata.obs["mixed"] = np.array([1, "a", None], dtype=object)
    with pytest.raises(UnsupportedColumnTypeError, match=r"column 'mixed'"):
        getattr(adata, f"write_{diskfmt}")(
            diskfmt_store, use=[A.obsm["good"], A.obs["mixed"]]
        )
    read: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert "good" not in read.obsm
    assert "mixed" not in read.obs


@pytest.mark.parametrize(
    ("values", "unwriteable"),
    [
        pytest.param([1, None, 3], False, id="int-none"),
        pytest.param([True, None, False], False, id="bool-none"),
        pytest.param(["a", None, "b"], False, id="str-none"),
        pytest.param([1, "a", None], True, id="mixed"),
        pytest.param([b"a", None, b"b"], True, id="bytes-none"),
        pytest.param([None, None, None], True, id="all-missing"),
    ],
)
@pytest.mark.parametrize("where", ["obs", "raw.var"])
def test_unwriteable_object_column(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    values: list[object],
    unwriteable: bool,
    where: str,
) -> None:
    """`unwriteable()` is true exactly when a default write raises for an object column."""
    col = np.array(values, dtype=object)
    if where == "obs":
        adata = _adata_with_obs(col=col)
    else:
        var = pd.DataFrame({"col": col}, index=[f"g{i}" for i in range(3)])
        adata = ad.AnnData(np.zeros((2, 3)), var=var)
        adata.raw = adata
        adata.var = pd.DataFrame(index=adata.var_names)
    assert adata.unwriteable() is unwriteable
    write = getattr(adata, f"write_{diskfmt}")
    if unwriteable:
        with pytest.raises(UnsupportedColumnTypeError, match=r"column 'col'"):
            write(diskfmt_store)
        return
    write(diskfmt_store)
    read: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    if where == "obs":
        frame = read.obs
    else:
        assert read.raw is not None
        frame = read.raw.var
    assert isinstance(frame, pd.DataFrame)
    assert frame["col"].isna().tolist() == pd.isna(col).tolist()


def test_unwriteable_outer_concat_strings(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Strings with missing values from an outer join (object on pandas < 3) are writable."""
    adata = ad.concat(
        [
            ad.AnnData(obs=pd.DataFrame({"s": ["x", "y"]}, index=["a0", "a1"])),
            ad.AnnData(obs=pd.DataFrame({"t": ["z", "z"]}, index=["b0", "b1"])),
        ],
        join="outer",
    )
    assert isinstance(adata.obs, pd.DataFrame)
    assert adata.obs["s"].isna().tolist() == [False, False, True, True]
    assert not adata.unwriteable()
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    read: ad.AnnData = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert isinstance(read.obs, pd.DataFrame)
    assert read.obs["s"].dtype == "category"
