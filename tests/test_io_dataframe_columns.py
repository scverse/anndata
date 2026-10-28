"""Which dataframe columns can be written, how, and how failures are reported.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import contextlib
import functools
import re
from contextlib import contextmanager
from typing import TYPE_CHECKING

import dask
import dask.array as da
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr
from zarr.storage import MemoryStore

import anndata as ad
from anndata._io.specs.dataframe_columns import (
    UnsupportedColumnTypeError,
    normalize_column,
)
from anndata._io.specs.registry import _REGISTRY, IORegistryError, IOSpec
from anndata._io.specs.validate import iter_unwriteable, raise_if_unwriteable
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, gen_adata

if TYPE_CHECKING:
    from collections.abc import Callable, Generator
    from pathlib import Path
    from typing import Literal


IDX = pd.Index([f"c{i}" for i in range(3)])


def _series(values: object, dtype: str | None = None) -> pd.Series:
    return pd.Series(values, index=IDX, dtype=dtype)


def _adata(**obs_columns: pd.Series) -> ad.AnnData:
    return ad.AnnData(obs=pd.DataFrame(obs_columns, index=IDX))


def _write(adata: ad.AnnData, diskfmt: str, store: object, **kwargs) -> ad.AnnData:
    getattr(adata, f"write_{diskfmt}")(
        store, convert_strings_to_categoricals=False, **kwargs
    )
    return getattr(ad, f"read_{diskfmt}")(store)


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
            index=IDX,
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
        lambda: _series(pd.to_datetime(["2020-01-01", None, "2021-01-01"])),
        IORegistryError,
        r"datetime/timedelta columns are not supported yet, see scverse/anndata#2238",
    ),
    "datetime64-tz": (
        lambda: _series(
            pd.to_datetime(["2020-01-01", None, "2021-01-01"]).tz_localize("UTC")
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
        read = _write(adata, diskfmt, diskfmt_store).obs["col"]
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
    store_type = "h5" if diskfmt == "h5ad" else "zarr"

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
    store_type = "h5" if diskfmt == "h5ad" else "zarr"
    with pd.option_context("future.distinguish_nan_and_na", distinguish):
        values = pd.arrays.FloatingArray(
            np.array([1.5, 0.0, np.nan], dtype=dtype.lower()),
            np.array([False, True, False]),
        )
        data_before, mask_before = values._data.copy(), values._mask.copy()
        assert normalize_column(values, path="col") is values
        np.testing.assert_array_equal(values._data, data_before)
        np.testing.assert_array_equal(values._mask, mask_before)

        adata = _adata(col=pd.Series(values, index=IDX))
        errors = list(iter_unwriteable(adata, store_type=store_type))
        assert [path for path, _ in errors] == ["/obs/col"]
        assert isinstance(errors[0][1], IORegistryError)
        assert "no on-disk encoding for nullable floats" in _notes(errors[0][1])
        np.testing.assert_array_equal(adata.obs["col"]._values._data, data_before)
        np.testing.assert_array_equal(adata.obs["col"]._values._mask, mask_before)
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
def test_normalize_column_idempotent_and_non_mutating(values: object) -> None:
    original = values.copy()
    once = normalize_column(values, path="col")
    twice = normalize_column(once, path="col")
    assert type(twice) is type(once)
    assert pd.array(twice).equals(pd.array(once))
    assert type(values) is type(original)
    assert pd.array(values).equals(pd.array(original))


def test_unsupported_object_column_error_attributes() -> None:
    with pytest.raises(UnsupportedColumnTypeError) as exc_info:
        normalize_column(np.array([1, "a", None], dtype=object), path="/obs/col")
    e = exc_info.value
    assert isinstance(e, TypeError)
    assert (e.path, e.dtype, e.kind) == ("/obs/col", "object", "mixed-integer")
    assert str(e).startswith(
        "Cannot write column '/obs/col' (dtype object, inferred as 'mixed-integer'): "
    )
    assert e.hint in str(e)


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
    read = getattr(ad, f"read_{diskfmt}")(diskfmt_store).obs["col"]
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


# Preflight validation


def _bad_everywhere() -> ad.AnnData:
    adata = _adata(
        good=_series([1, 2, 3]),
        period=_series(pd.period_range("2020-01", periods=3, freq="M")),
    )
    adata.var["dt"] = pd.to_datetime([])
    adata.obsm["df"] = pd.DataFrame(
        {"mixed": np.array([1, "a", None], dtype=object)}, index=IDX
    )
    adata.uns["nested"] = {"ok": 1, "bad": object()}
    return adata


def test_preflight_reports_all_problems(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    """All problems are raised together, and nothing is written."""
    path = tmp_path / f"adata.{diskfmt}"
    adata = _bad_everywhere()
    paths = ["/obs/period", "/var/dt", "/uns/nested/bad", "/obsm/df/mixed"]

    assert [p for p, _ in iter_unwriteable(adata)] == paths
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(ExceptionGroup) as exc_info,
    ):
        getattr(adata, f"write_{diskfmt}")(path)
    group = exc_info.value
    assert group.message == f"4 elements cannot be written: {', '.join(paths)}"
    assert len(group.exceptions) == 4
    assert all(
        f"Element: {p}" in _notes(e)
        for p, e in zip(paths, group.exceptions, strict=True)
    )
    assert not path.exists()

    # leaves keep their types for `except*`
    registry_misses = group.subgroup(IORegistryError)
    conversion_failures = group.subgroup(TypeError)
    assert registry_misses is not None
    assert len(registry_misses.exceptions) == 3
    assert conversion_failures is not None
    assert len(conversion_failures.exceptions) == 1


def test_preflight_off_fails_at_first_problem(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    adata = _bad_everywhere()
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(IORegistryError, match=r"PeriodArray"),
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata.{diskfmt}")


def test_preflight_leaves_existing_destination_untouched(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    good = _adata(col=_series([1, 2, 3]))
    getattr(good, f"write_{diskfmt}")(path)
    bad = _bad_everywhere()
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(ExceptionGroup),
    ):
        getattr(bad, f"write_{diskfmt}")(path)
    read = getattr(ad, f"read_{diskfmt}")(path)
    assert list(read.obs.columns) == ["col"]
    assert read.obs["col"].tolist() == [1, 2, 3]


def test_preflight_passes_clean_adata(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    adata = gen_adata((10, 8), **GEN_ADATA_NO_XARRAY_ARGS)
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    read = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert read.shape == adata.shape


def test_write_elem_bypasses_preflight(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """`write_elem` is not validated up front, but still normalises columns."""
    if diskfmt == "h5ad":
        store = h5py.File(diskfmt_store, "w")
    else:
        store = zarr.open_group(diskfmt_store, mode="w")
    df = pd.DataFrame({"obj": np.array([1, None, 3], dtype=object)}, index=IDX)
    with ad.settings.override(validate_before_write=True):
        ad.io.write_elem(store, "df", df)
        assert ad.io.read_elem(store["df"])["obj"].dtype == "Int64"
        # no aggregation: the first problem is raised
        with pytest.raises(IORegistryError) as exc_info:
            ad.io.write_elem(store, "uns", {"a": object(), "b": object()})
    assert not isinstance(exc_info.value, ExceptionGroup)


def test_preflight_raw(*, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path) -> None:
    adata = ad.AnnData(np.zeros((3, 2)), obs=pd.DataFrame(index=IDX))
    adata.raw = adata
    adata.raw.var["period"] = pd.period_range("2020-01", periods=2, freq="M")
    assert [p for p, _ in iter_unwriteable(adata)] == ["/raw/var/period"]
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(IORegistryError, match=r"PeriodArray"),
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata.{diskfmt}")


def test_preflight_dataframe_structure() -> None:
    import pyarrow as pa

    adata = _adata(col=_series([1, 2, 3]))
    adata.obs["_index"] = [1, 2, 3]
    assert [(p, type(e)) for p, e in iter_unwriteable(adata)] == [("/obs", ValueError)]

    adata = _adata(col=_series([1, 2, 3]))
    adata.obs = pd.DataFrame(np.ones((3, 2)), columns=["dup", "dup"], index=IDX)
    errors = list(iter_unwriteable(adata))
    assert [p for p, _ in errors] == ["/obs"]
    assert "repeated column names" in str(errors[0][1])

    # a value comparison: writer-time only
    adata = _adata(col=_series([1, 2, 3]))
    adata.obs.index.name = "col"
    assert list(iter_unwriteable(adata)) == []
    with pytest.raises(ValueError, match=r"index.name"):
        adata.write_zarr(MemoryStore())

    with ad.settings.override(restrict_index_types=False):
        adata = ad.AnnData(
            np.zeros((3, 1)),
            obs=pd.DataFrame(
                {"x": [1, 2, 3]},
                index=pd.Index(
                    pd.arrays.ArrowExtensionArray(
                        pa.array([1, 2, 3], type=pa.timestamp("ns"))
                    )
                ),
            ),
        )
    errors = list(iter_unwriteable(adata, store_type="zarr"))
    assert [(p, type(e)) for p, e in errors] == [("/obs/_index", IORegistryError)]
    assert adata.unwriteable()


def _write_period(f, k, v, *, _writer, dataset_kwargs):
    _writer.write_elem(f, k, np.asarray(v.asi8), dataset_kwargs=dataset_kwargs)


@pytest.fixture
def period_writer() -> Generator[None]:
    """Temporarily register a writer for `PeriodArray`."""
    keys = [
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


def _fail_on_compute(dsk, keys, **kwargs):
    pytest.fail("the preflight check computed a dask array")


def _fail_on_packing(*args, **kwargs):
    pytest.fail("the preflight check packed an awkward array")


def test_preflight_does_not_compute_dask() -> None:
    adata = ad.AnnData(
        da.zeros((4, 3), chunks=(2, 3)), obs=pd.DataFrame(index=list("abcd"))
    )
    adata.layers["L"] = da.ones((4, 3), chunks=(1, 1))
    adata.obsm["m"] = da.zeros((4, 2), chunks=(2, 2))

    with dask.config.set(scheduler=_fail_on_compute):
        assert list(iter_unwriteable(adata)) == []
        raise_if_unwriteable(adata)


@pytest.mark.parametrize("sparse_format", [None, "csr", "csc"])
def test_preflight_does_not_read_backed(
    *, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sparse_format: str | None
) -> None:
    """Backed datasets (X, layers, raw.X; also a Raw nested in uns) are never read."""
    import scipy.sparse as sp

    path = tmp_path / "backed.h5ad"
    X = np.arange(12.0).reshape(4, 3)
    if sparse_format is not None:
        X = getattr(sp, f"{sparse_format}_matrix")(X)
    adata = ad.AnnData(X, obs=pd.DataFrame(index=list("abcd")))
    adata.raw = adata
    adata.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    reads: list[str] = []
    getitem = h5py.Dataset.__getitem__
    monkeypatch.setattr(
        h5py.Dataset,
        "__getitem__",
        lambda self, *args, **kw: (reads.append(self.name), getitem(self, *args, **kw))[
            1
        ],
    )
    view = backed[:2]
    assert view.is_view
    nested = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
    nested.uns["raw"] = view.raw
    for elem in (backed, view, nested):
        assert list(iter_unwriteable(elem, store_type="h5")) == []
        raise_if_unwriteable(elem, store_type="h5")
    assert reads == []


def test_preflight_does_not_pack_awkward(*, monkeypatch: pytest.MonkeyPatch) -> None:
    ak = pytest.importorskip("awkward")
    adata = ad.AnnData(np.zeros((3, 1)), obs=pd.DataFrame(index=IDX))
    adata.obsm["ragged"] = ak.Array([[1, 2], [], [3]])
    adata.uns["ragged"] = ak.Array([[1.0], [2.0, 3.0]])

    monkeypatch.setattr(ak, "to_packed", _fail_on_packing)
    monkeypatch.setattr(ak, "to_buffers", _fail_on_packing)
    assert list(iter_unwriteable(adata)) == []
    assert not adata.unwriteable()


def test_preflight_view_uses_own_metadata(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """A view is validated on what it would write, not on its parent's columns."""
    parent = _adata(col=_series(np.array([1, "a", None], dtype=object)))
    parent.layers["L"] = np.zeros((3, 0))
    assert parent.unwriteable()
    view = parent[:1]
    assert view.is_view
    assert view.obs["col"].tolist() == [1]
    assert list(iter_unwriteable(view)) == []
    with ad.settings.override(validate_before_write=True):
        read = _write(view, diskfmt, diskfmt_store)
    assert read.obs["col"].tolist() == [1]


def test_preflight_multiindex_aggregates() -> None:
    with ad.settings.override(restrict_index_types=False):
        adata = ad.AnnData(
            np.zeros((2, 1)),
            obs=pd.DataFrame(
                {"period": pd.period_range("2020-01", periods=2, freq="M")},
                index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2)]),
            ),
        )
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [
        ("/obs/_index", ValueError),
        ("/obs/period", IORegistryError),
    ]
    assert "MultiIndex" in str(errors[0][1])
    assert "Element: /obs/_index" in _notes(errors[0][1])
    assert adata.unwriteable()
    with pytest.raises(ExceptionGroup) as exc_info:
        raise_if_unwriteable(adata)
    assert len(exc_info.value.exceptions) == 2


CATEGORY_FAMILIES: dict[str, tuple[Callable[[], pd.Categorical], bool]] = {
    # name: (factory, has a faithful on-disk representation)
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
            pd.array(
                ["a", "b", "a"],
                dtype=pd.ArrowDtype(pytest.importorskip("pyarrow").string()),
            )
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
    """The walker judges categories the way the writer stores them (`to_numpy()`)."""
    make, registered = CATEGORY_FAMILIES[family]
    store_type = "h5" if diskfmt == "h5ad" else "zarr"
    adata = _adata(cat=pd.Series(make(), index=IDX))
    errors = list(iter_unwriteable(adata, store_type=store_type))
    assert [p for p, _ in errors] == ([] if registered else ["/obs/cat/categories"])
    assert adata.unwriteable(store_type=store_type) is not registered
    if not registered:
        # the writer refuses them too, on both formats
        assert isinstance(errors[0][1], TypeError)
        assert "Cannot write categories" in str(errors[0][1])
        with pytest.raises(TypeError, match=r"Cannot write categories"):
            _write(adata, diskfmt, diskfmt_store)
        return
    if family == "datetime" and diskfmt == "h5ad":
        # h5py rejects the values; not a type problem
        with pytest.raises(TypeError, match=r"No conversion path"):
            _write(adata, diskfmt, diskfmt_store)
        return
    read = _write(adata, diskfmt, diskfmt_store).obs["cat"]
    assert read.dtype == "category"
    assert read.isna().tolist() == adata.obs["cat"].isna().tolist()


class _ArrayApiNoDlpack:
    """An array-API object without DLPack: `to_writeable` cannot convert it."""

    shape, ndim, size, dtype, device = (2,), 1, 2, np.dtype("f8"), "cpu"

    def __array_namespace__(self, *, api_version=None):
        return np

    def to_device(self, device, /, *, stream=None):
        return self

    def __getitem__(self, k, /):
        return self


def test_preflight_matches_writer_dispatch_for_leaves() -> None:
    """Standalone Series/Index, non-dict mappings and non-DLPack array-API objects miss."""
    from collections import UserDict

    adata = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
    adata.uns["series"] = pd.Series([1, 2])
    adata.uns["index"] = pd.Index(["a", "b"])
    adata.uns["userdict"] = UserDict({"a": 1})
    adata.uns["no_dlpack"] = _ArrayApiNoDlpack()
    adata.uns["ok"] = {"nested": np.zeros(2), "none": None, "list": [1, 2]}
    errors = list(iter_unwriteable(adata))
    assert [p for p, _ in errors] == [
        "/uns/series",
        "/uns/index",
        "/uns/userdict",
        "/uns/no_dlpack",
    ]
    assert all(isinstance(e, IORegistryError) for _, e in errors)
    for key in ("series", "index", "userdict", "no_dlpack"):
        with pytest.raises(IORegistryError):
            ad.io.write_elem(
                zarr.open_group(MemoryStore(), mode="w"), key, adata.uns[key]
            )


def test_preflight_invalid_column_name_keeps_destination(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    good = _adata(col=_series([1, 2, 3]))
    getattr(good, f"write_{diskfmt}")(path)

    bad = _adata(col=_series([1, 2, 3]))
    bad.obs[0] = [4, 5, 6]
    errors = list(iter_unwriteable(bad))
    assert [(p, type(e)) for p, e in errors] == [("/obs", TypeError)]
    assert "invalid key" in str(errors[0][1])
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(TypeError, match=r"invalid key"),
    ):
        getattr(bad, f"write_{diskfmt}")(path)
    read = getattr(ad, f"read_{diskfmt}")(path)
    assert list(read.obs.columns) == ["col"]
    assert read.obs["col"].tolist() == [1, 2, 3]


def test_preflight_aggregates_non_2d_with_column(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    adata = _adata(period=_series(pd.period_range("2020-01", periods=3, freq="M")))
    with pytest.warns(UserWarning, match=r"must be 2-dimensional"):
        adata.X = np.zeros((3, 0, 2))
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(ExceptionGroup) as exc_info,
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata.{diskfmt}")
    # writer order: obs before X
    assert [type(e) for e in exc_info.value.exceptions] == [IORegistryError, ValueError]
    assert exc_info.value.message.endswith("/obs/period, /X")
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(ValueError, match=r"must be 2-dimensional"),
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata2.{diskfmt}")


def test_preflight_does_not_read_aligned_mappings_of_views(
    *, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Aligned mappings of a view are inspected via the parent, never subset."""
    with h5py.File(tmp_path / "arrays.h5", "w") as f:
        f["m"] = np.zeros((4, 2))
        f["p"] = np.zeros((4, 4))
        adata = ad.AnnData(
            np.zeros((4, 4)),
            obs=pd.DataFrame(index=list("abcd")),
            var=pd.DataFrame(index=list("wxyz")),
        )
        adata.obsm["m"], adata.varm["m"] = f["m"], f["m"]
        adata.obsp["p"], adata.varp["p"] = f["p"], f["p"]
        adata.obsm["df"] = pd.DataFrame({"x": [1, 2, 3, 4]}, index=adata.obs_names)
        reads: list[str] = []
        getitem = h5py.Dataset.__getitem__
        monkeypatch.setattr(
            h5py.Dataset,
            "__getitem__",
            lambda self, *a, **kw: (reads.append(self.name), getitem(self, *a, **kw))[
                1
            ],
        )
        view = adata[[0, 2], [0, 2]]
        assert view.is_view
        assert list(iter_unwriteable(view, store_type="h5")) == []
        raise_if_unwriteable(view, store_type="h5")
        assert reads == []


def test_preflight_does_not_materialize_range_categories() -> None:
    cat = pd.Categorical.from_codes([0, 1, 0], categories=pd.RangeIndex(100_000))
    adata = _adata(cat=pd.Series(cat, index=IDX))
    assert list(iter_unwriteable(adata)) == []
    assert not adata.unwriteable()
    assert "_data" not in adata.obs["cat"].cat.categories._cache
    assert "_data" not in cat.categories._cache


def test_preflight_multiindex_with_colliding_name() -> None:
    with ad.settings.override(restrict_index_types=False):
        adata = ad.AnnData(
            np.zeros((2, 1)),
            obs=pd.DataFrame(
                {
                    "col": [1, 2],
                    "period": pd.period_range("2020-01", periods=2, freq="M"),
                },
                index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2)]),
            ),
        )
    adata.obs.index.name = "col"
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [
        ("/obs/col", ValueError),
        ("/obs/period", IORegistryError),
    ]
    assert "MultiIndex" in str(errors[0][1])
    assert adata.unwriteable()


@pytest.mark.parametrize("where", ["layers", "uns", "nested-uns"])
def test_preflight_slash_keys_keep_destination(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path, where: str
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    good = _adata(col=_series([1, 2, 3]))
    getattr(good, f"write_{diskfmt}")(path)

    bad = _adata(col=_series([1, 2, 3]))
    if where == "layers":
        bad.layers["bad/key"] = np.zeros((3, 0))
    elif where == "uns":
        bad.uns["bad/key"] = 1
    else:
        bad.uns["outer"] = {"bad/key": 1}
    store_type = "h5" if diskfmt == "h5ad" else "zarr"
    errors = list(iter_unwriteable(bad, store_type=store_type))
    assert [type(e) for _, e in errors] == [ValueError]
    assert "Forward slashes" in str(errors[0][1])
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(ValueError, match=r"Forward slashes"),
    ):
        getattr(bad, f"write_{diskfmt}")(path)
    read = getattr(ad, f"read_{diskfmt}")(path)
    assert list(read.obs.columns) == ["col"]
    assert read.obs["col"].tolist() == [1, 2, 3]


def test_preflight_slash_keys_follow_h5ad_setting(*, tmp_path: Path) -> None:
    """With the h5ad slash setting off, `uns` keys with slashes are allowed (with a warning)."""
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["bad/key"] = 1
    with ad.settings.override(disallow_forward_slash_in_h5ad=False):
        assert list(iter_unwriteable(adata, store_type="h5")) == []
        assert [p for p, _ in iter_unwriteable(adata, store_type="zarr")] == [
            "/uns/bad/key"
        ]
        with (
            ad.settings.override(validate_before_write=True),
            pytest.warns(FutureWarning, match=r"Forward slashes"),
        ):
            adata.write_h5ad(tmp_path / "adata.h5ad")
    # layers keys: refused regardless of the setting
    adata.layers["bad/key"] = np.zeros((3, 0))
    with ad.settings.override(disallow_forward_slash_in_h5ad=False):
        assert [p for p, _ in iter_unwriteable(adata, store_type="h5")] == ["/layers"]


@pytest.mark.parametrize("kind", ["OrderedDict", "defaultdict", "Counter"])
def test_preflight_accepts_dict_subclasses(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore, kind: str
) -> None:
    import collections

    if kind == "OrderedDict":
        value = collections.OrderedDict(a=1, b=np.zeros(2))
    elif kind == "defaultdict":
        value = collections.defaultdict(int, a=1)
    else:
        value = collections.Counter({"a": 2})
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["sub"] = value
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        read = _write(adata, diskfmt, diskfmt_store)
    assert read.uns["sub"]["a"] == value["a"]


# Instrumented "never materialise" harness


class _ForbiddenAccess:
    """Stand-in for a method the preflight check must never call."""

    def __init__(self, what: str) -> None:
        self.what = what

    def __call__(self, *args, **kwargs):
        pytest.fail(f"the preflight check touched payload data: {self.what}")


class _ForbiddenOnNonEmpty:
    """Stand-in for an index method that is only allowed on empty indexes."""

    def __init__(self, what: str, original) -> None:
        self.what = what
        self.original = original

    def __get__(self, index, owner=None):
        # bind like a method when patched onto the class
        return self if index is None else functools.partial(self, index)

    def __call__(self, index, *args, **kwargs):
        if len(index) > 0:
            pytest.fail(f"the preflight check materialised an index: {self.what}")
        return self.original(index, *args, **kwargs)


@contextmanager
def forbid_payload_access() -> Generator[None]:
    """Fail on anything that would subset, copy, materialise or read payload data.

    Patched: dataframe item access, slicing, `take`, `copy` and `to_numpy`; `.array`,
    `.to_numpy()` and `__array__` on any non-empty index, and `RangeIndex._data`;
    anndata's `_subset`/`as_view`; `h5py.Dataset` item reads and `read_direct`;
    `zarr.Array` item reads; dask `compute`/`persist`.
    `np.ndarray.copy`/`__array__` are C-level and cannot be patched.
    """
    from unittest import mock

    import dask.base
    import zarr.core.array as zarr_array
    from pandas.core.indexing import _iLocIndexer, _LocIndexer

    import anndata._core.anndata as core_anndata
    from anndata._core import aligned_mapping, views

    patches = [
        mock.patch.object(
            pd.DataFrame, "__getitem__", _ForbiddenAccess("DataFrame.__getitem__")
        ),
        mock.patch.object(pd.DataFrame, "_slice", _ForbiddenAccess("DataFrame._slice")),
        mock.patch.object(pd.DataFrame, "take", _ForbiddenAccess("DataFrame.take")),
        mock.patch.object(pd.DataFrame, "copy", _ForbiddenAccess("DataFrame.copy")),
        mock.patch.object(
            _iLocIndexer, "__getitem__", _ForbiddenAccess("DataFrame.iloc")
        ),
        mock.patch.object(
            _LocIndexer, "__getitem__", _ForbiddenAccess("DataFrame.loc")
        ),
        mock.patch.object(
            pd.Index,
            "array",
            property(
                _ForbiddenOnNonEmpty(
                    "Index.array",
                    pd.Index.array.__get__,
                )
            ),
        ),
        mock.patch.object(
            pd.Index,
            "to_numpy",
            _ForbiddenOnNonEmpty("Index.to_numpy", pd.Index.to_numpy),
        ),
        mock.patch.object(
            pd.Index,
            "__array__",
            _ForbiddenOnNonEmpty("Index.__array__", pd.Index.__array__),
        ),
        mock.patch.object(
            pd.RangeIndex,
            "_data",
            property(
                _ForbiddenOnNonEmpty("RangeIndex._data", pd.RangeIndex._data.__get__)
            ),
        ),
        mock.patch.object(
            aligned_mapping, "_subset", _ForbiddenAccess("_subset (aligned mapping)")
        ),
        mock.patch.object(
            core_anndata, "_subset", _ForbiddenAccess("_subset (AnnData)")
        ),
        mock.patch.object(views, "as_view", _ForbiddenAccess("as_view")),
        mock.patch.object(
            aligned_mapping, "as_view", _ForbiddenAccess("as_view (aligned mapping)")
        ),
        mock.patch.object(
            pd.DataFrame, "to_numpy", _ForbiddenAccess("DataFrame.to_numpy")
        ),
        mock.patch.object(
            h5py.Dataset, "__getitem__", _ForbiddenAccess("h5py.Dataset.__getitem__")
        ),
        mock.patch.object(
            h5py.Dataset, "read_direct", _ForbiddenAccess("h5py.Dataset.read_direct")
        ),
        mock.patch.object(
            zarr_array.Array, "__getitem__", _ForbiddenAccess("zarr.Array.__getitem__")
        ),
        mock.patch.object(dask.base, "compute", _ForbiddenAccess("dask.compute")),
        mock.patch.object(
            dask.base.DaskMethodsMixin,
            "compute",
            _ForbiddenAccess("dask Array.compute"),
        ),
        mock.patch.object(
            dask.base.DaskMethodsMixin,
            "persist",
            _ForbiddenAccess("dask Array.persist"),
        ),
    ]
    # `np.ndarray.copy`/`__array__` belong to a C type and cannot be patched
    with contextlib.ExitStack() as stack:
        for p in patches:
            stack.enter_context(p)
        yield


N_PAYLOAD = 5_000


def _payload_adata(*, dlpack: bool = True) -> ad.AnnData:
    """A richly populated AnnData; `dlpack=False` leaves out the non-exportable array-API subject."""
    rng = np.random.default_rng(0)
    n = N_PAYLOAD
    index = pd.Index([f"cell{i}" for i in range(n)])
    obj = np.where(rng.random(n) < 0.1, None, rng.integers(0, 9, n))
    obs = pd.DataFrame(
        {
            "cat": pd.Categorical(rng.choice(["a", "b", "c"], n)),
            "num": rng.random(n),
            "obj": pd.Series(obj, dtype=object, index=index),
            "Int64": pd.array(rng.integers(0, 9, n), dtype="Int64"),
            "str": pd.array([f"s{i}" for i in range(n)], dtype="string"),
        },
        index=index,
    )
    adata = ad.AnnData(sp.csr_matrix(rng.random((n, 20))), obs=obs)
    adata.layers["dense"] = np.zeros((n, 20), dtype=np.float32)
    adata.obsm["numeric_df"] = pd.DataFrame(
        rng.random((n, 10)), index=adata.obs_names, columns=[f"c{i}" for i in range(10)]
    )
    adata.obsm["object_df"] = pd.DataFrame(
        {"o": pd.Series(obj, dtype=object, index=adata.obs_names)},
        index=adata.obs_names,
    )
    adata.obsm["arr"] = np.zeros((n, 3))
    adata.varm["arr"] = np.zeros((20, 3))
    adata.obsp["p"] = sp.csr_matrix((n, n))
    adata.varp["p"] = np.zeros((20, 20))
    adata.uns["nested"] = {
        "arr": np.zeros(n),
        "scalar": 1,
        "none": None,
        "list": [1, 2],
    }
    adata.uns["table"] = pd.DataFrame({"number": np.zeros(n)})  # RangeIndex
    adata.uns["cat_index_table"] = pd.DataFrame(
        {"x": np.zeros(3)}, index=pd.CategoricalIndex(["a", "b", "a"])
    )
    adata.uns["range_cat"] = pd.DataFrame(
        {"c": pd.Categorical.from_codes([0, 1, 0], categories=pd.RangeIndex(n))},
        index=list("abc"),
    )
    # the writer compares a same-named column with the index; the preflight must not
    adata.uns["named_collision"] = pd.DataFrame(
        {"id": np.arange(n)}, index=pd.RangeIndex(n, name="id")
    )
    adata.uns["named_index"] = pd.DataFrame(
        {"v": np.zeros(n)}, index=pd.RangeIndex(n, name="row")
    )
    if dlpack:  # cannot be written, so left out of subjects that get written to disk
        adata.uns["dlpack"] = _ArrayApiNoExport(np.zeros(n))
    adata.raw = adata
    return adata


class _ArrayApiNoExport:
    """An array-API array with DLPack support that must not be exported by the walk."""

    def __init__(self, data: np.ndarray) -> None:
        self._data = data
        self.shape, self.ndim, self.size = data.shape, data.ndim, data.size
        self.dtype, self.device = data.dtype, "cpu"

    def __array_namespace__(self, *, api_version=None):
        return np

    def to_device(self, device, /, *, stream=None):
        return self

    def __getitem__(self, k, /):
        return self

    def __dlpack_device__(self):
        return self._data.__dlpack_device__()

    def __dlpack__(self, **kwargs):
        pytest.fail("the preflight check exported an array-API array via DLPack")


def _subject(kind: str, tmp_path: Path) -> object:  # noqa: PLR0911
    adata = _payload_adata()
    if kind == "memory":
        return adata
    if kind == "view":
        return adata[np.arange(0, N_PAYLOAD, 2), [0, 2, 4]]
    if kind == "clean-view-of-dirty-parent":
        dirty = np.array(["bad", *([1] * (N_PAYLOAD - 1))], dtype=object)
        adata.obsm["dirty"] = pd.DataFrame({"o": dirty}, index=adata.obs_names)
        adata.varm["dirty"] = pd.DataFrame(
            {"o": np.array(["bad", *([1] * 19)], dtype=object)}, index=adata.var_names
        )
        assert adata.unwriteable()
        return adata[np.arange(1, N_PAYLOAD), np.arange(1, 20)]
    if kind.startswith("backed"):
        adata = _payload_adata(dlpack=False)
        fmt = kind.split("-")[1]
        if fmt == "dense":
            adata.X = adata.X.toarray()
        elif fmt == "csc":
            adata.X = sp.csc_matrix(adata.X)
        path = tmp_path / f"{kind}.h5ad"
        adata.write_h5ad(path)
        backed = ad.read_h5ad(path, backed="r")
        return backed[[0, 2, 4], [0, 2]]
    if kind == "dask":
        adata.X = da.from_array(adata.X.toarray(), chunks=(1000, 20))
        adata.layers["dense"] = da.zeros((N_PAYLOAD, 20), chunks=(1000, 20))
        return adata
    if kind == "awkward":
        ak = pytest.importorskip("awkward")
        adata.obsm["ragged"] = ak.Array([[i] * (i % 3) for i in range(N_PAYLOAD)])
        return adata
    if kind == "nested-raw":
        outer = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
        outer.uns["raw"] = adata.raw
        outer.uns["adata"] = adata
        outer.uns["view"] = adata[::2]
        return outer
    msg = f"unknown subject {kind}"
    raise ValueError(msg)


SUBJECTS = [
    "memory",
    "view",
    "clean-view-of-dirty-parent",
    "backed-dense",
    "backed-csr",
    "backed-csc",
    "dask",
    "awkward",
    "nested-raw",
]


@pytest.mark.parametrize("kind", SUBJECTS)
@pytest.mark.parametrize("store_type", ["h5", "zarr"])
def test_preflight_never_materialises(
    *, kind: str, store_type: Literal["h5", "zarr"], tmp_path: Path
) -> None:
    subject = _subject(kind, tmp_path)
    with forbid_payload_access():
        assert list(iter_unwriteable(subject, store_type=store_type)) == []
        if isinstance(subject, ad.AnnData):
            assert not subject.unwriteable(store_type=store_type)
            raise_if_unwriteable(subject, store_type=store_type)


def test_harness_catches_payload_access(*, tmp_path: Path) -> None:
    """The harness itself: each forbidden operation is detected."""
    adata = _payload_adata(dlpack=False)
    adata.uns["dlpack"] = _ArrayApiNoExport(np.zeros(3))
    view = adata[::2]
    backed_path = tmp_path / "h.h5ad"
    _payload_adata(dlpack=False).write_h5ad(backed_path)
    backed = ad.read_h5ad(backed_path, backed="r")
    probes = {
        "getitem": lambda: adata.obs["num"],
        "iloc": lambda: adata.obs.iloc[:, 0],
        "take": lambda: adata.obs.take([0]),
        "copy": lambda: adata.obs.copy(),  # noqa: PLW0108
        "index.array": lambda: adata.obs.index.array,
        "index.to_numpy": lambda: adata.obs.index.to_numpy(),  # noqa: PLW0108
        "rangeindex": lambda: adata.uns["table"].index._data,
        "view obsm": lambda: view.obsm["numeric_df"],
        "backed read": lambda: backed.file["obs/num"][()],
        "to_numpy": lambda: adata.obs.to_numpy(),  # noqa: PLW0108
        "read_direct": lambda: backed.file["obs/num"].read_direct(np.empty(N_PAYLOAD)),
        "dask compute": lambda: da.zeros(3).compute(),
        "dlpack": lambda: np.from_dlpack(adata.uns["dlpack"]),
    }
    for probe in probes.values():
        with forbid_payload_access(), pytest.raises(pytest.fail.Exception):
            probe()
    with forbid_payload_access():
        assert len(adata.obs.index[:0].to_numpy()) == 0  # empty indexes are allowed


def test_preflight_reports_cycles() -> None:
    """A container that contains itself is reported, not recursed into forever."""
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["self"] = adata.uns
    adata.uns["adata"] = adata
    adata.uns["shared"] = shared = {"x": 1}
    adata.uns["shared_again"] = shared  # ordinary shared references are fine
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [
        ("/uns/self", ValueError),
        ("/uns/adata", ValueError),
    ]
    assert all("reference cycle" in str(e) for _, e in errors)
    assert adata.unwriteable()
    with pytest.raises(ExceptionGroup):
        raise_if_unwriteable(adata)


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        pytest.param(7, TypeError, id="int"),
        pytest.param(None, TypeError, id="none"),
        pytest.param(("a", "b"), TypeError, id="tuple"),
        pytest.param("/uns/foo", None, id="absolute-inside"),
        pytest.param("/else/key", ValueError, id="absolute-outside"),
        pytest.param("plain", None, id="plain"),
    ],
)
def test_preflight_mapping_key_rules(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    tmp_path: Path,
    key: object,
    expected: type | None,
) -> None:
    """Mapping keys are judged the way `Writer.write_elem` treats them."""
    path = tmp_path / f"adata.{diskfmt}"
    getattr(_adata(col=_series([1, 2, 3])), f"write_{diskfmt}")(path)
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns[key] = 1
    store_type = "h5" if diskfmt == "h5ad" else "zarr"
    with ad.settings.override(disallow_forward_slash_in_h5ad=False):
        errors = list(iter_unwriteable(adata, store_type=store_type))
        if expected is None:
            assert errors == []
            with ad.settings.override(validate_before_write=True):
                getattr(adata, f"write_{diskfmt}")(path)
            read = getattr(ad, f"read_{diskfmt}")(path)
            assert read.uns[str(key).removeprefix("/uns/")] == 1
            return
        assert [type(e) for _, e in errors] == [expected]
        with (
            ad.settings.override(validate_before_write=True),
            pytest.raises(expected),
        ):
            getattr(adata, f"write_{diskfmt}")(path)
    read = getattr(ad, f"read_{diskfmt}")(path)
    assert list(read.obs.columns) == ["col"]
    assert "uns" not in read.uns


@pytest.mark.parametrize(
    "categories",
    [
        pytest.param(pd.period_range("2020-01", periods=2, freq="M"), id="period"),
        pytest.param(pd.IntervalIndex.from_tuples([(0, 1), (1, 2)]), id="interval"),
        pytest.param(
            pd.to_datetime(["2020-01-01", "2021-01-01"]).tz_localize("UTC"), id="tz"
        ),
    ],
)
def test_preflight_categorical_index_policy(categories: pd.Index) -> None:
    """A categorical index goes through the same category policy as a column."""
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["table"] = pd.DataFrame(
        {"x": [1.0, 2.0, 3.0]},
        index=pd.CategoricalIndex(
            pd.Categorical.from_codes([0, 1, 0], categories=categories)
        ),
    )
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [
        ("/uns/table/_index/categories", TypeError)
    ]
    assert "Cannot write categories" in str(errors[0][1])
    assert adata.unwriteable()
    good = pd.DataFrame(
        {"x": [1.0, 2.0]}, index=pd.CategoricalIndex(["a", "b"], categories=["a", "b"])
    )
    adata.uns["table"] = good
    assert list(iter_unwriteable(adata)) == []


def test_write_entry_points_do_not_read_view_payloads_before_preflight(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The `no_write_dataset_2d` guard inspects a view through its parent too."""
    with h5py.File(tmp_path / "arrays.h5", "w") as f:
        f["m"] = np.zeros((3, 2))
        adata = _adata(
            col=_series([1, 2, 3]),
            period=_series(pd.period_range("2020-01", periods=3, freq="M")),
        )
        adata.obsm["m"] = f["m"]
        adata.varm["m"] = np.zeros((0, 1))
        view = adata[[0, 2]]
        reads: list[str] = []
        getitem = h5py.Dataset.__getitem__
        monkeypatch.setattr(
            h5py.Dataset,
            "__getitem__",
            lambda self, *a, **kw: (reads.append(self.name), getitem(self, *a, **kw))[
                1
            ],
        )
        with (
            ad.settings.override(validate_before_write=True),
            pytest.raises(IORegistryError),
        ):
            getattr(view, f"write_{diskfmt}")(
                tmp_path / f"out.{diskfmt}", convert_strings_to_categoricals=False
            )
        assert reads == []


def test_preflight_rejects_dataset2d() -> None:
    xr = pytest.importorskip("xarray")
    from anndata.experimental.backed import Dataset2D

    adata = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
    empty = Dataset2D(
        xr.Dataset(coords={"obs_names": ("obs_names", np.array(list("ab")))})
    )
    adata.obs = empty
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [("/obs", IORegistryError)]
    assert adata.unwriteable()


def test_preflight_view_object_columns_in_mappings(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Object columns of a view's obsm/varm dataframes are judged on the view's rows."""
    adata = ad.AnnData(
        np.zeros((3, 3)),
        obs=pd.DataFrame(index=IDX),
        var=pd.DataFrame(index=list("xyz")),
    )
    adata.obsm["t"] = pd.DataFrame(
        {"o": np.array([1, "bad", None], dtype=object)}, index=adata.obs_names
    )
    adata.varm["t"] = pd.DataFrame(
        {"o": np.array(["bad", 2, None], dtype=object)}, index=adata.var_names
    )
    assert [p for p, _ in iter_unwriteable(adata)] == ["/obsm/t/o", "/varm/t/o"]
    view = adata[[0], [1, 2]]
    assert list(iter_unwriteable(view)) == []
    with ad.settings.override(validate_before_write=True):
        read = _write(view, diskfmt, diskfmt_store)
    assert read.obsm["t"]["o"].tolist() == [1]
    assert read.varm["t"]["o"].isna().tolist() == [False, True]


@pytest.mark.parametrize("mapping", ["obsm", "varm", "obsp", "varp"])
def test_preflight_none_key_only_for_layers(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path, mapping: str
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    getattr(_adata(col=_series([1, 2, 3])), f"write_{diskfmt}")(path)
    adata = _adata(col=_series([1, 2, 3]))
    n = 3 if mapping.startswith("obs") else 0
    getattr(adata, mapping)[None] = np.zeros((n, n if mapping.endswith("p") else 2))
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [(f"/{mapping}/None", TypeError)]
    assert adata.unwriteable()
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(TypeError, match=r"invalid key"),
    ):
        getattr(adata, f"write_{diskfmt}")(path)
    read = getattr(ad, f"read_{diskfmt}")(path)
    assert list(read.obs.columns) == ["col"]


def test_preflight_absolute_keys_normalised_for_children(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["/uns/foo"] = {"/uns/foo/bar": 1}
    with ad.settings.override(disallow_forward_slash_in_h5ad=False):
        assert list(iter_unwriteable(adata)) == []
        with ad.settings.override(validate_before_write=True):
            read = _write(adata, diskfmt, diskfmt_store)
    assert read.uns["foo"]["bar"] == 1


def test_preflight_reports_sequence_cycles() -> None:
    adata = _adata(col=_series([1, 2, 3]))
    cycle: list = []
    cycle.append(cycle)
    adata.uns["cycle"] = cycle
    indirect: dict = {}
    indirect["seq"] = [1, [2, indirect]]
    adata.uns["indirect"] = indirect
    shared = [1, 2]
    adata.uns["shared"] = {"a": shared, "b": [shared, shared]}
    via_dict: list = [{"back": None}]
    via_dict[0]["back"] = via_dict
    adata.uns["via_dict"] = via_dict
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [
        ("/uns/cycle/0", ValueError),
        ("/uns/indirect/seq/1/1", ValueError),
        ("/uns/via_dict/0/back", ValueError),
    ]
    assert adata.unwriteable()


def test_preflight_accepts_sequences_of_writable_members(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """Sequence members are not validated individually: the list writer stores the whole list."""
    adata = _adata(col=_series([1, 2, 3]))
    adata.uns["pairs"] = [(1, 2), (3, 4)]
    adata.uns["nested"] = [[1, 2], [3, 4]]
    adata.uns["strings"] = ["a", "b"]
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        read = _write(adata, diskfmt, diskfmt_store)
    np.testing.assert_array_equal(read.uns["pairs"], [[1, 2], [3, 4]])
    np.testing.assert_array_equal(read.uns["nested"], [[1, 2], [3, 4]])
    assert list(read.uns["strings"]) == ["a", "b"]
    # a standalone tuple has no writer
    adata.uns["tuple"] = (1, 2)
    assert [(p, type(e)) for p, e in iter_unwriteable(adata)] == [
        ("/uns/tuple", IORegistryError)
    ]
