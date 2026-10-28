"""Which dataframe columns can be written, how, and how failures are reported.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import copy
import pickle
import re
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr

import anndata as ad
from anndata._io.specs.dataframe_columns import (
    UnsupportedColumnTypeError,
    normalize_column,
)
from anndata._io.specs.registry import _REGISTRY, IORegistryError, IOSpec
from anndata._io.specs.validate import iter_unwriteable

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Mapping
    from pathlib import Path
    from typing import Any, Literal

    from pandas.api.extensions import ExtensionArray
    from zarr.storage import MemoryStore

    from anndata._io.specs.registry import Writer


def _idx() -> pd.Index:
    """A fresh index for every frame: pandas < 3 shares the Index object between frames."""
    return pd.Index([f"c{i}" for i in range(3)])


def _series(values: object, dtype: str | None = None) -> pd.Series:
    return pd.Series(values, index=_idx(), dtype=dtype)


def _adata(**obs_columns: pd.Series) -> ad.AnnData:
    return ad.AnnData(obs=pd.DataFrame(obs_columns, index=_idx()))


def _write(adata: ad.AnnData, diskfmt: str, store: object, **kwargs) -> ad.AnnData:
    getattr(adata, f"write_{diskfmt}")(
        store, convert_strings_to_categoricals=False, **kwargs
    )
    return getattr(ad, f"read_{diskfmt}")(store)


def _store_type(diskfmt: Literal["h5ad", "zarr"]) -> Literal["h5", "zarr"]:
    return "h5" if diskfmt == "h5ad" else "zarr"


def _df(elem: object) -> pd.DataFrame:
    """Narrow an element (`adata.obs`, a read element, …) to a dataframe."""
    assert isinstance(elem, pd.DataFrame)
    return elem


def _notes(e: BaseException) -> str:
    return "\n".join(getattr(e, "__notes__", []))


# families with a writer and their on-disk dtype; factories give fresh input
SUPPORTED: dict[str, tuple[Callable[[], pd.Series], str]] = {
    "int64": (lambda: _series([1, 2, 3]), "int64"),
    "float64-nan": (lambda: _series([1.5, np.nan, 2.0]), "float64"),
    "bool": (lambda: _series([True, False, True]), "bool"),
    "Int64": (lambda: _series([1, None, 3], "Int64"), "Int64"),
    "boolean": (lambda: _series([True, None, False], "boolean"), "boolean"),
    "string": (lambda: _series(["a", None, "a"], "string"), "string"),
    "category": (lambda: _series(["a", "b", None]).astype("category"), "category"),
    "category-ordered": (
        lambda: pd.Series(
            pd.Categorical(["lo", "hi", "lo"], categories=["lo", "hi"], ordered=True),
            index=_idx(),
        ),
        "category",
    ),
    "object-int-none": (
        lambda: _series(np.array([1, None, 3], dtype=object)),
        "Int64",
    ),
    "object-float-none": (
        lambda: _series(np.array([1.5, None, 2.0], dtype=object)),
        "float64",
    ),
    "object-bool-none": (
        lambda: _series(np.array([True, None, False], dtype=object)),
        "boolean",
    ),
}

# families without a writer: error type and the hint in its notes
UNSUPPORTED: dict[str, tuple[Callable[[], pd.Series], type[Exception], str]] = {
    "Float64": (
        lambda: _series([1.5, None, 2.0], "Float64"),
        IORegistryError,
        r"no on-disk encoding for nullable floats.*na_value=np\.nan",
    ),
    "Float32": (
        lambda: _series([1.5, None, 2.0], "Float32"),
        IORegistryError,
        r"no on-disk encoding for nullable floats",
    ),
    "datetime64": (
        lambda: _series(pd.to_datetime(["2020-01-01", "NaT", "2021-01-01"])),
        IORegistryError,
        r"datetime/timedelta columns are not supported yet, see scverse/anndata#2238",
    ),
    "datetime64-tz": (
        lambda: _series(
            pd.to_datetime(["2020-01-01", "NaT", "2021-01-01"]).tz_localize("UTC")
        ),
        IORegistryError,
        r"scverse/anndata#2238",
    ),
    "timedelta64": (
        lambda: _series(pd.to_timedelta(["1h", None, "2h"])),
        IORegistryError,
        r"scverse/anndata#2238",
    ),
    "int64-pyarrow": (
        lambda: _series([1, None, 3], "int64[pyarrow]"),
        IORegistryError,
        r"pyarrow-backed columns are not supported.*numpy_nullable",
    ),
    "period": (
        lambda: _series(pd.period_range("2020-01", periods=3, freq="M")),
        IORegistryError,
        r"Convert the column to a numpy, nullable integer/boolean, string, or categorical",
    ),
    "interval": (
        lambda: _series(pd.interval_range(0, 3)),
        IORegistryError,
        r"Convert the column to a numpy",
    ),
    "sparse": (
        lambda: _series(pd.arrays.SparseArray([0.0, 1.0, 0.0])),
        IORegistryError,
        r"Convert the column to a numpy",
    ),
    "object-mixed": (
        lambda: _series(np.array([1, "a", None], dtype=object)),
        UnsupportedColumnTypeError,
        r"",
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
    adata = _adata(col=col)
    assert list(iter_unwriteable(adata)) == []
    assert not adata.unwriteable()

    with ad.settings.override(
        validate_before_write=True, allow_write_nullable_strings=True
    ):
        read = _df(_write(adata, diskfmt, diskfmt_store).obs)["col"]
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
    store_type = _store_type(diskfmt)

    adata = _adata(col=make())
    errors = list(iter_unwriteable(adata, store_type=store_type))
    assert [path for path, _ in errors] == ["/obs/col"]
    assert isinstance(errors[0][1], exc_type)
    assert re.search(hint, _notes(errors[0][1]) + str(errors[0][1]), re.DOTALL)
    assert adata.unwriteable(store_type=store_type)

    # the writer raises the same type
    adata = _adata(col=make())
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(exc_type) as exc_info,
    ):
        _write(adata, diskfmt, diskfmt_store)
    assert re.search(hint, _notes(exc_info.value) + str(exc_info.value), re.DOTALL)
    assert re.search(r"'col'", str(exc_info.value) + _notes(exc_info.value))

    # one problem is raised as is, not as a group
    adata = _adata(col=make())
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(exc_type) as exc_info,
    ):
        _write(adata, diskfmt, diskfmt_store)
    assert "Element: /obs/col" in _notes(exc_info.value)


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
    store_type = _store_type(diskfmt)
    with pd.option_context("future.distinguish_nan_and_na", distinguish):
        values = pd.arrays.FloatingArray(
            np.array([1.5, 0.0, np.nan], dtype=dtype.lower()),
            np.array([False, True, False]),
        )
        before = values.copy()
        assert normalize_column(values, path="col") is values
        pd.testing.assert_extension_array_equal(values, before)
        assert values.isna().tolist() == [False, True, False]

        adata = _adata(col=pd.Series(values, index=_idx()))
        errors = list(iter_unwriteable(adata, store_type=store_type))
        assert [path for path, _ in errors] == ["/obs/col"]
        assert isinstance(errors[0][1], IORegistryError)
        assert "no on-disk encoding for nullable floats" in _notes(errors[0][1])
        pd.testing.assert_extension_array_equal(_df(adata.obs)["col"].array, before)
        with (
            ad.settings.override(validate_before_write=False),
            pytest.raises(IORegistryError) as exc_info,
        ):
            _write(adata, diskfmt, diskfmt_store)
        assert "na_value=np.nan" in _notes(exc_info.value)


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
    assert _notes(restored) == "Element: /obs/col"
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
    col = _series(["a", None, "a"], f"string[{storage}]")
    adata = _adata(col=col)
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
    read = _df(getattr(ad, f"read_{diskfmt}")(diskfmt_store).obs)["col"]
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
    adata = _adata(
        repeated=_series(["a", "a", "b"]),
        unique=_series(["a", "b", "c"]),
        obj=_series(np.array([1, None, 3], dtype=object)),
        cat=_series(["x", "y", "x"]).astype("category"),
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
    adata = _adata(period=_series(pd.period_range("2020-01", periods=3, freq="M")))
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
    store_type = _store_type(diskfmt)
    adata = _adata(cat=pd.Series(make(), index=_idx()))
    errors = list(iter_unwriteable(adata, store_type=store_type))
    assert [p for p, _ in errors] == ([] if accepted else ["/obs/cat/categories"])
    assert adata.unwriteable(store_type=store_type) is not accepted
    if not accepted:
        assert isinstance(errors[0][1], TypeError)
        assert "Cannot write categories" in str(errors[0][1])
        with ad.settings.override(validate_before_write=True), pytest.raises(TypeError):
            _write(adata, diskfmt, diskfmt_store)
    if diskfmt == "h5ad" and (not accepted or family in {"datetime", "timedelta"}):
        # h5py refuses object and datetime64/timedelta64 data
        with pytest.raises(TypeError, match=r"convert|conversion path"):
            _write(adata, diskfmt, diskfmt_store)
        return
    read = _df(_write(adata, diskfmt, diskfmt_store).obs)["cat"]
    assert read.dtype == "category"
    assert read.isna().tolist() == adata.obs["cat"].isna().tolist()
    if not accepted:  # zarr stores what it cannot encode as strings
        assert pd.api.types.is_string_dtype(read.cat.categories.dtype)
