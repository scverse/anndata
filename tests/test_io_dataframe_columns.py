"""Which dataframe columns can be written, how, and how failures are reported.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import dask
import dask.array as da
import h5py
import numpy as np
import pandas as pd
import pytest
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


# Column families with a writer: how they come back from disk.
# Every entry is a factory so each test gets fresh, unmodified input.
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

# Column families without a writer: the error type and the hint in its notes.
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

    # writing without preflight raises the same type from inside `write_dataframe`
    adata = _adata(col=make())
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(exc_type) as exc_info,
    ):
        _write(adata, diskfmt, diskfmt_store)
    assert re.search(hint, _notes(exc_info.value) + str(exc_info.value), re.DOTALL)
    assert re.search(r"'col'", str(exc_info.value) + _notes(exc_info.value))

    # a single problem found by the preflight is raised as is, not as a group
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
            # pandas < 3 without `future.infer_string`: the writer asks for an explicit choice
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

    # the leaves keep their types, so `except*` works
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
        # two problems, but no aggregation: the first one is raised
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

    adata = _adata(col=_series([1, 2, 3]))
    adata.obs.index.name = "col"
    errors = list(iter_unwriteable(adata))
    assert [p for p, _ in errors] == ["/obs"]
    assert "index.name" in str(errors[0][1])

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


@pytest.fixture
def period_writer() -> Generator[None]:
    """Temporarily register a writer for `PeriodArray`."""
    keys = [
        (store_type, pd.arrays.PeriodArray, frozenset())
        for store_type in (h5py.Group, zarr.Group)
    ]
    assert not any(key in _REGISTRY.write for key in keys)

    def write_period(f, k, v, *, _writer, dataset_kwargs):
        _writer.write_elem(f, k, np.asarray(v.asi8), dataset_kwargs=dataset_kwargs)

    for store_type, src_type, _ in keys:
        _REGISTRY.register_write(store_type, src_type, IOSpec("period-test", "0.1.0"))(
            write_period
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


def test_preflight_does_not_compute_dask() -> None:
    adata = ad.AnnData(
        da.zeros((4, 3), chunks=(2, 3)), obs=pd.DataFrame(index=list("abcd"))
    )
    adata.layers["L"] = da.ones((4, 3), chunks=(1, 1))
    adata.obsm["m"] = da.zeros((4, 2), chunks=(2, 2))

    def no_scheduler(dsk, keys, **kwargs):
        pytest.fail("the preflight check computed a dask array")

    with dask.config.set(scheduler=no_scheduler):
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

    def no_packing(*args, **kwargs):
        pytest.fail("the preflight check packed an awkward array")

    monkeypatch.setattr(ak, "to_packed", no_packing)
    monkeypatch.setattr(ak, "to_buffers", no_packing)
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
    # name: (factory, writable on both formats per the registry)
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
    assert [p for p, _ in errors] == ([] if registered else ["/obs/cat"])
    assert adata.unwriteable(store_type=store_type) is not registered
    if not registered:
        # The writer hands these categories to the string writer as arbitrary objects:
        # h5py refuses them, zarr silently stores `str(value)`. The preflight (like
        # `unwriteable()` before it) rejects both rather than mirror the corruption.
        if diskfmt == "h5ad":
            with pytest.raises(TypeError):
                _write(adata, diskfmt, diskfmt_store)
        else:
            read = _write(adata, diskfmt, diskfmt_store).obs["cat"]
            assert read.cat.categories.dtype != adata.obs["cat"].cat.categories.dtype
        return
    if family == "datetime" and diskfmt == "h5ad":
        # a value-level failure in h5py's writer, outside the preflight's scope
        with pytest.raises(TypeError, match=r"No conversion path"):
            _write(adata, diskfmt, diskfmt_store)
        return
    read = _write(adata, diskfmt, diskfmt_store).obs["cat"]
    assert read.dtype == "category"
    assert read.isna().tolist() == adata.obs["cat"].isna().tolist()


def test_preflight_matches_writer_dispatch_for_leaves() -> None:
    """Standalone Series/Index, non-dict mappings and non-DLPack array-API objects miss."""
    from collections import UserDict

    class ArrayApiNoDlpack:
        shape, ndim, size, dtype, device = (2,), 1, 2, np.dtype("f8"), "cpu"

        def __array_namespace__(self, *, api_version=None):
            return np

        def to_device(self, device, /, *, stream=None):
            return self

        def __getitem__(self, k, /):
            return self

    adata = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
    adata.uns["series"] = pd.Series([1, 2])
    adata.uns["index"] = pd.Index(["a", "b"])
    adata.uns["userdict"] = UserDict({"a": 1})
    adata.uns["no_dlpack"] = ArrayApiNoDlpack()
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
    # traversal order is that of the writers: obs before X
    assert [type(e) for e in exc_info.value.exceptions] == [IORegistryError, ValueError]
    assert exc_info.value.message.endswith("/obs/period, /X")
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(ValueError, match=r"must be 2-dimensional"),
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata2.{diskfmt}")
