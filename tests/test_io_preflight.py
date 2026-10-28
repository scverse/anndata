"""The preflight check behind `settings.validate_before_write` and `AnnData.unwriteable`.

See https://github.com/scverse/anndata/issues/1068
"""

from __future__ import annotations

import contextlib
import functools
import inspect
from contextlib import contextmanager
from typing import TYPE_CHECKING

import dask.array as da
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr
from zarr.storage import MemoryStore

import anndata as ad
from anndata._io.specs.dataframe_columns import UnsupportedColumnTypeError
from anndata._io.specs.registry import IORegistryError
from anndata._io.specs.validate import (
    iter_unwriteable,
    raise_if_unwriteable,
    raise_if_unwriteable_parts,
)
from anndata._io.utils import _resolve_use
from anndata.acc import A
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    adata_with_obs,
    as_dataframe,
    exception_notes,
    gen_adata,
    obs_column,
    obs_index,
    store_type_of,
    write_and_read,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Generator
    from pathlib import Path
    from typing import Any, Literal

    from anndata.acc import AdRef, MapAcc, RefAcc


def _bad_everywhere() -> ad.AnnData:
    adata = adata_with_obs(
        good=obs_column([1, 2, 3]),
        period=obs_column(pd.period_range("2020-01", periods=3, freq="M")),
    )
    adata.var["dt"] = pd.to_datetime([])
    adata.obsm["df"] = pd.DataFrame(
        {"mixed": np.array([1, "a", None], dtype=object)}, index=obs_index()
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
        f"Element: {p}" in exception_notes(e)
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
    good = adata_with_obs(col=obs_column([1, 2, 3]))
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
    store: h5py.File | zarr.Group
    if diskfmt == "h5ad":
        store = h5py.File(diskfmt_store, "w")
    else:
        store = zarr.open_group(diskfmt_store, mode="w")
    df = pd.DataFrame({"obj": np.array([1, None, 3], dtype=object)}, index=obs_index())
    with ad.settings.override(validate_before_write=True):
        ad.io.write_elem(store, "df", df)
        assert as_dataframe(ad.io.read_elem(store["df"]))["obj"].dtype == "Int64"
        # no aggregation: the first problem is raised
        with pytest.raises(IORegistryError) as exc_info:
            ad.io.write_elem(store, "uns", {"a": object(), "b": object()})
    assert not isinstance(exc_info.value, ExceptionGroup)


def test_preflight_raw(*, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path) -> None:
    adata = ad.AnnData(np.zeros((3, 2)), obs=pd.DataFrame(index=obs_index()))
    adata.raw = adata
    assert adata.raw is not None
    adata.raw.var["period"] = pd.period_range("2020-01", periods=2, freq="M")
    assert [p for p, _ in iter_unwriteable(adata)] == ["/raw/var/period"]
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(IORegistryError, match=r"PeriodArray"),
    ):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata.{diskfmt}")


def _strings_with_missing(source: str) -> ad.AnnData:
    if source == "concat":
        # an outer join leaves the other object's rows missing (object columns on pandas < 3)
        return ad.concat(
            [
                ad.AnnData(obs=pd.DataFrame({"s": ["x", "y"]}, index=["a0", "a1"])),
                ad.AnnData(obs=pd.DataFrame({"t": ["z", "z"]}, index=["b0", "b1"])),
            ],
            join="outer",
        )
    # pandas >= 3 infers `str` unless told otherwise
    adata = ad.AnnData(
        np.zeros((3, 2)),
        obs=pd.DataFrame({"s": obs_column(["a", None, "b"], "object")}),
        var=pd.DataFrame({"s": pd.Series(["a", None], index=["x", "y"], dtype=object)}),
    )
    adata.raw = adata
    return adata


@pytest.mark.parametrize("validate", [False, True], ids=["validate=F", "validate=T"])
@pytest.mark.parametrize("source", ["concat", "object"])
def test_unwriteable_judges_strings_as_the_writers_convert_them(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    source: str,
    validate: bool,
) -> None:
    """Strings with missing values in `obs`, `var` and `raw.var` are written as categoricals."""
    adata = _strings_with_missing(source)
    assert all(
        pd.api.types.infer_dtype(col) == "string" and col.hasnans
        for _, col in as_dataframe(adata.obs).items()
    )
    store_type = store_type_of(diskfmt)
    assert not adata.unwriteable()
    assert not adata.unwriteable(store_type=store_type)
    if source == "object":
        # without the conversion, the writers refuse them
        paths = ["/obs/s", "/var/s", "/raw/var/s"]
        errors = list(iter_unwriteable(adata, store_type=store_type))
        assert [p for p, _ in errors] == paths
        expected = ExceptionGroup if validate else UnsupportedColumnTypeError
        with (
            ad.settings.override(validate_before_write=validate),
            pytest.raises(expected),
        ):
            write_and_read(adata, diskfmt, diskfmt_store)
    with ad.settings.override(validate_before_write=validate):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    read = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert all(
        isinstance(col.dtype, pd.CategoricalDtype)
        for _, col in as_dataframe(read.obs).items()
    )


def test_unwriteable_does_not_convert_nested_or_aligned_strings() -> None:
    """The writers convert only the top-level `obs`, `var` and `raw.var`."""
    strings = obs_column(["a", None, "a"], "object")
    inner = adata_with_obs(s=strings)
    assert not inner.unwriteable()
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.uns["inner"] = inner
    adata.obsm["df"] = pd.DataFrame({"s": strings})
    errors = list(iter_unwriteable(adata, strings_to_categoricals=True))
    assert [(p, type(e)) for p, e in errors] == [
        ("/uns/inner/obs/s", UnsupportedColumnTypeError),
        ("/obsm/df/s", UnsupportedColumnTypeError),
    ]
    assert adata.unwriteable()
    with pytest.raises(UnsupportedColumnTypeError):
        adata.write_zarr(MemoryStore())


def test_preflight_dataframe_structure() -> None:
    import pyarrow as pa

    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.obs["_index"] = [1, 2, 3]
    assert [(p, type(e)) for p, e in iter_unwriteable(adata)] == [("/obs", ValueError)]

    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.obs = pd.DataFrame(np.ones((3, 2)), columns=["dup", "dup"], index=obs_index())
    errors = list(iter_unwriteable(adata))
    assert [p for p, _ in errors] == ["/obs"]
    assert "repeated column names" in str(errors[0][1])

    # a value comparison: writer-time only
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.obs = as_dataframe(adata.obs).rename_axis("col")
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


def test_preflight_view_uses_own_metadata(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """A view is validated on what it would write, not on its parent's columns."""
    parent = adata_with_obs(col=obs_column(np.array([1, "a", None], dtype=object)))
    parent.layers["L"] = np.zeros((3, 0))
    assert parent.unwriteable()
    view = parent[:1]
    assert view.is_view
    assert view.obs["col"].tolist() == [1]
    assert list(iter_unwriteable(view)) == []
    with ad.settings.override(validate_before_write=True):
        read = write_and_read(view, diskfmt, diskfmt_store)
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
    assert "Element: /obs/_index" in exception_notes(errors[0][1])
    assert adata.unwriteable()
    with pytest.raises(ExceptionGroup) as exc_info:
        raise_if_unwriteable(adata)
    assert len(exc_info.value.exceptions) == 2


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


class _DataFrameSubclass(pd.DataFrame):
    @property
    def _constructor(self) -> type[_DataFrameSubclass]:
        return _DataFrameSubclass


class _CategoricalSubclass(pd.Categorical):
    pass


class _ArraySubclass(np.ndarray):
    pass


@pytest.mark.parametrize(
    "make",
    [
        pytest.param(lambda: _DataFrameSubclass({"a": [1, 2]}), id="DataFrame"),
        pytest.param(lambda: _CategoricalSubclass(["a", "b"]), id="Categorical"),
        pytest.param(lambda: np.zeros(3).view(_ArraySubclass), id="ndarray"),
    ],
)
def test_preflight_dispatches_subclasses_like_writer(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    make: Callable[[], object],
) -> None:
    """Containers are only walked if their exact type has a writer, as in `write_elem`."""
    adata = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame(index=list("ab")))
    adata.uns["x"] = make()
    errors = list(iter_unwriteable(adata))
    assert [(p, type(e)) for p, e in errors] == [("/uns/x", IORegistryError)]
    with pytest.raises(IORegistryError):
        write_and_read(adata, diskfmt, diskfmt_store)


def test_preflight_invalid_column_name_keeps_destination(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    good = adata_with_obs(col=obs_column([1, 2, 3]))
    getattr(good, f"write_{diskfmt}")(path)

    bad = adata_with_obs(col=obs_column([1, 2, 3]))
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
    adata = adata_with_obs(
        period=obs_column(pd.period_range("2020-01", periods=3, freq="M"))
    )
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


def test_preflight_does_not_read_aligned_mappings_of_views(*, tmp_path: Path) -> None:
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
        view = adata[[0, 2], [0, 2]]
        assert view.is_view
        with forbid_payload_access():
            assert list(iter_unwriteable(view, store_type="h5")) == []
            raise_if_unwriteable(view, store_type="h5")


def test_preflight_does_not_materialize_range_categories() -> None:
    cat = pd.Categorical.from_codes([0, 1, 0], categories=pd.RangeIndex(100_000))
    adata = adata_with_obs(cat=pd.Series(cat, index=obs_index()))
    with forbid_payload_access():  # forbids `RangeIndex._data`
        assert list(iter_unwriteable(adata)) == []
        assert not adata.unwriteable()


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
    good = adata_with_obs(col=obs_column([1, 2, 3]))
    getattr(good, f"write_{diskfmt}")(path)

    bad = adata_with_obs(col=obs_column([1, 2, 3]))
    if where == "layers":
        bad.layers["bad/key"] = np.zeros((3, 0))
    elif where == "uns":
        bad.uns["bad/key"] = 1
    else:
        bad.uns["outer"] = {"bad/key": 1}
    store_type = store_type_of(diskfmt)
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
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
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

    value: dict[str, Any]
    if kind == "OrderedDict":
        value = collections.OrderedDict(a=1, b=np.zeros(2))
    elif kind == "defaultdict":
        value = collections.defaultdict(int, a=1)
    else:
        value = collections.Counter({"a": 2})
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.uns["sub"] = value
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        read = write_and_read(adata, diskfmt, diskfmt_store)
    assert read.uns["sub"]["a"] == value["a"]


# Instrumented "never materialise" harness


class _ForbiddenAccess:
    """Stand-in for a method the preflight check must never call."""

    def __init__(self, what: str) -> None:
        self.what = what

    def __call__(self, *args, **kwargs):
        pytest.fail(f"the preflight check touched payload data: {self.what}")


class _ForbiddenAbove:
    """Stand-in for a method only allowed on objects of at most `limit` elements."""

    def __init__(
        self,
        what: str,
        original: Callable[..., object],
        *,
        limit: int = 0,
        size: Callable[[Any], int] = len,
    ) -> None:
        self.what, self.original, self.limit, self.size = what, original, limit, size

    def __get__(self, obj: object, owner: type | None = None) -> object:
        # bind like a method when patched onto the class
        return self if obj is None else functools.partial(self, obj)

    def __call__(self, obj: object, *args: object, **kwargs: object) -> object:
        if self.size(obj) > self.limit:
            pytest.fail(f"the preflight check materialised payload data: {self.what}")
        return self.original(obj, *args, **kwargs)


# column labels of the payload frames are shorter, all their row indexes longer
_MAX_LABELS = 16


def _n_elements(matrix: sp.csr_matrix | sp.csr_array) -> int:
    return matrix.shape[0] * matrix.shape[1]


@contextmanager
def forbid_payload_access() -> Generator[None]:
    """Fail on anything that would subset, copy, materialise or read payload data.

    Patched: dataframe item access, slicing, `take`, `copy` and `to_numpy`; `.array`,
    `.to_numpy()` and `__array__` on any non-empty index, and `RangeIndex._data`;
    `.values` and iteration of indexes, and `__array__` of extension arrays, longer
    than `_MAX_LABELS`; `toarray`/`todense` of non-empty scipy sparse matrices;
    anndata's `_subset`/`as_view`; `h5py.Dataset` item reads and `read_direct`;
    `zarr.Array` item reads; dask `compute`/`persist`.
    Not patched: `np.ndarray.copy`/`__array__` (C-level), and `Index._values`, which
    pandas reads without copying to build every column Series, and which only
    materialises a `RangeIndex` (via the patched `_data`).
    """
    from unittest import mock

    import dask.base
    import zarr.core.array as zarr_array
    from pandas.core.arrays.masked import BaseMaskedArray
    from pandas.core.indexing import _iLocIndexer, _LocIndexer

    import anndata._core.anndata as core_anndata
    from anndata._core import aligned_mapping, views

    patches: list[contextlib.AbstractContextManager[object]] = [
        mock.patch.object(
            pd.Index,
            "values",
            property(
                _ForbiddenAbove(
                    "Index.values",
                    inspect.getattr_static(pd.Index, "values").fget,
                    limit=_MAX_LABELS,
                )
            ),
        ),
        *(
            mock.patch.object(
                cls,
                "__iter__",
                _ForbiddenAbove(
                    f"iter({cls.__name__})", cls.__iter__, limit=_MAX_LABELS
                ),
            )
            for cls in (pd.Index, pd.RangeIndex)
        ),
        *(
            mock.patch.object(
                cls,
                "__array__",
                _ForbiddenAbove(
                    f"np.asarray({cls.__name__})",
                    inspect.getattr_static(cls, "__array__"),
                    limit=_MAX_LABELS,
                ),
            )
            for cls in (
                pd.Categorical,
                BaseMaskedArray,
                pd.arrays.NumpyExtensionArray,
                pd.arrays.ArrowExtensionArray,
            )
        ),
        *(
            mock.patch.object(
                cls,
                name,
                _ForbiddenAbove(
                    f"{cls.__name__}.{name}", getattr(cls, name), size=_n_elements
                ),
            )
            for cls in (sp.csr_matrix, sp.csc_matrix, sp.csr_array, sp.csc_array)
            for name in ("toarray", "todense")
        ),
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
                _ForbiddenAbove(
                    "Index.array", inspect.getattr_static(pd.Index, "array").__get__
                )
            ),
        ),
        mock.patch.object(
            pd.Index,
            "to_numpy",
            _ForbiddenAbove("Index.to_numpy", pd.Index.to_numpy),
        ),
        mock.patch.object(
            pd.Index,
            "__array__",
            _ForbiddenAbove("Index.__array__", pd.Index.__array__),
        ),
        mock.patch.object(
            pd.RangeIndex,
            "_data",
            property(
                _ForbiddenAbove(
                    "RangeIndex._data",
                    inspect.getattr_static(pd.RangeIndex, "_data").__get__,
                )
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
    obj = rng.integers(0, 9, n).astype(object)
    obj[rng.random(n) < 0.1] = None
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
        assert isinstance(adata.X, sp.csr_matrix)
        fmt = kind.split("-")[1]
        if fmt == "dense":
            adata.X = adata.X.toarray()
        elif fmt == "csc":
            adata.X = sp.csc_matrix(adata.X)
        path = tmp_path / f"{kind}.h5ad"
        with ad.settings.override(allow_write_nullable_strings=True):
            adata.write_h5ad(path)
        backed = ad.read_h5ad(path, backed="r")
        return backed[[0, 2, 4], [0, 2]]
    if kind == "dask":
        assert isinstance(adata.X, sp.csr_matrix)
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


def _caught(probe: Callable[[], object]) -> bool:
    with forbid_payload_access():
        try:
            probe()
        except pytest.fail.Exception:
            return True
    return False


def test_harness_catches_payload_access(*, tmp_path: Path) -> None:
    """The harness itself: each forbidden operation is detected."""
    adata = _payload_adata(dlpack=False)
    adata.uns["dlpack"] = _ArrayApiNoExport(np.zeros(3))
    view = adata[::2]
    backed_path = tmp_path / "h.h5ad"
    with ad.settings.override(allow_write_nullable_strings=True):
        _payload_adata(dlpack=False).write_h5ad(backed_path)
    backed = ad.read_h5ad(backed_path, backed="r")
    obs, var = as_dataframe(adata.obs), as_dataframe(adata.var)
    table = as_dataframe(adata.uns["table"])
    columns = {name: obs[name].array for name in ("cat", "Int64", "str")}
    X = adata.X
    assert isinstance(X, sp.csr_matrix)
    X_array = sp.csr_array(X)
    backed_num = backed.file["obs/num"]
    assert isinstance(backed_num, h5py.Dataset)
    probes: dict[str, Callable[[], object]] = {
        "getitem": lambda: obs["num"],
        "iloc": lambda: obs.iloc[:, 0],
        "take": lambda: obs.take([0]),
        "copy": lambda: obs.copy(),  # noqa: PLW0108
        "index.array": lambda: obs.index.array,
        "index.to_numpy": lambda: obs.index.to_numpy(),  # noqa: PLW0108
        "index.values": lambda: obs.index.values,  # noqa: PD011
        "var index.values": lambda: var.index.values,  # noqa: PD011
        "iter(index)": lambda: list(obs.index),
        "iter(RangeIndex)": lambda: list(table.index),
        "rangeindex": lambda: adata.uns["table"].index._data,
        **{
            f"np.asarray({name})": functools.partial(np.asarray, values)
            for name, values in columns.items()
        },
        "toarray": lambda: X.toarray(),  # noqa: PLW0108
        "todense": lambda: X.todense(),  # noqa: PLW0108
        "csr_array.toarray": lambda: X_array.toarray(),  # noqa: PLW0108
        "view obsm": lambda: view.obsm["numeric_df"],
        "backed read": lambda: backed_num[()],
        "to_numpy": lambda: obs.to_numpy(),  # noqa: PLW0108
        "read_direct": lambda: backed_num.read_direct(np.empty(N_PAYLOAD)),
        "dask compute": lambda: da.zeros(3).compute(),
        "dlpack": lambda: np.from_dlpack(adata.uns["dlpack"]),
    }
    assert [name for name, probe in probes.items() if not _caught(probe)] == []
    # metadata-sized access is allowed
    with forbid_payload_access():
        assert len(obs.index[:0].to_numpy()) == 0
        assert list(obs.columns) == ["cat", "num", "obj", "Int64", "str"]
        assert "num" in obs.columns


def test_preflight_reports_cycles() -> None:
    """A container that contains itself is reported, not recursed into forever."""
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
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
    getattr(adata_with_obs(col=obs_column([1, 2, 3])), f"write_{diskfmt}")(path)
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.uns[key] = 1
    store_type = store_type_of(diskfmt)
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
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
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
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path
) -> None:
    """The `no_write_dataset_2d` guard inspects a view through its parent too."""
    with h5py.File(tmp_path / "arrays.h5", "w") as f:
        f["m"] = np.zeros((3, 2))
        adata = adata_with_obs(
            col=obs_column([1, 2, 3]),
            period=obs_column(pd.period_range("2020-01", periods=3, freq="M")),
        )
        adata.obsm["m"] = f["m"]
        adata.varm["m"] = np.zeros((0, 1))
        view = adata[[0, 2]]
        with (
            forbid_payload_access(),
            ad.settings.override(validate_before_write=True),
            pytest.raises(IORegistryError),
        ):
            getattr(view, f"write_{diskfmt}")(
                tmp_path / f"out.{diskfmt}", convert_strings_to_categoricals=False
            )


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
        obs=pd.DataFrame(index=obs_index()),
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
        read = write_and_read(view, diskfmt, diskfmt_store)
    assert as_dataframe(read.obsm["t"])["o"].tolist() == [1]
    assert as_dataframe(read.varm["t"])["o"].isna().tolist() == [False, True]


@pytest.mark.parametrize("mapping", ["obsm", "varm", "obsp", "varp"])
def test_preflight_none_key_only_for_layers(
    *, diskfmt: Literal["h5ad", "zarr"], tmp_path: Path, mapping: str
) -> None:
    path = tmp_path / f"adata.{diskfmt}"
    getattr(adata_with_obs(col=obs_column([1, 2, 3])), f"write_{diskfmt}")(path)
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
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
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.uns["/uns/foo"] = {"/uns/foo/bar": 1}
    with ad.settings.override(disallow_forward_slash_in_h5ad=False):
        assert list(iter_unwriteable(adata)) == []
        with ad.settings.override(validate_before_write=True):
            read = write_and_read(adata, diskfmt, diskfmt_store)
    assert read.uns["foo"]["bar"] == 1


def test_preflight_reports_sequence_cycles() -> None:
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
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
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    adata.uns["pairs"] = [(1, 2), (3, 4)]
    adata.uns["nested"] = [[1, 2], [3, 4]]
    adata.uns["strings"] = ["a", "b"]
    assert list(iter_unwriteable(adata)) == []
    with ad.settings.override(validate_before_write=True):
        read = write_and_read(adata, diskfmt, diskfmt_store)
    np.testing.assert_array_equal(read.uns["pairs"], [[1, 2], [3, 4]])
    np.testing.assert_array_equal(read.uns["nested"], [[1, 2], [3, 4]])
    assert list(read.uns["strings"]) == ["a", "b"]
    # a standalone tuple has no writer
    adata.uns["tuple"] = (1, 2)
    assert [(p, type(e)) for p, e in iter_unwriteable(adata)] == [
        ("/uns/tuple", IORegistryError)
    ]


def test_preflight_use_checks_only_used_elements(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore
) -> None:
    """A partial write is judged on what it writes, not on the rest of the object."""
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    write_and_read(adata, diskfmt, diskfmt_store)
    adata.obs["period"] = pd.period_range("2020-01", periods=3, freq="M")
    adata.obsm["clean"] = np.ones((3, 2))
    assert adata.unwriteable()
    with ad.settings.override(validate_before_write=True):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store, use=A.obsm["clean"])
    read = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    np.testing.assert_array_equal(read.obsm["clean"], np.ones((3, 2)))
    assert list(as_dataframe(read.obs).columns) == ["col"]


@pytest.mark.parametrize("bad", ["column", "mapping"])
def test_preflight_use_refuses_unwritable_part(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore, bad: str
) -> None:
    """An unwritable part fails the partial write before any part is written."""
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    write_and_read(adata, diskfmt, diskfmt_store)
    period = pd.period_range("2020-01", periods=3, freq="M")
    adata.obsm["good"] = np.ones((3, 2))
    if bad == "column":
        adata.obs["period"] = period
        use, path = [A.obsm["good"], A.obs["period"]], "/obs/period"
    else:
        adata.obsm["bad"] = pd.DataFrame({"p": period}, index=adata.obs_names)
        use, path = [A.obsm["good"], A.obsm["bad"]], "/obsm/bad/p"
    with (
        ad.settings.override(validate_before_write=True),
        pytest.raises(IORegistryError) as exc_info,
    ):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store, use=use)
    assert f"Element: {path}" in exception_notes(exc_info.value)
    assert _keys_on_disk(diskfmt_store) == set()

    # without the preflight, the parts before the bad one are written
    with (
        ad.settings.override(validate_before_write=False),
        pytest.raises(IORegistryError),
    ):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store, use=use)
    assert "good" in _keys_on_disk(diskfmt_store)


def _keys_on_disk(store: Path | MemoryStore, group_name: str = "obsm") -> set[str]:
    """Read the store directly: a failed partial write may leave it unreadable or stale."""
    if isinstance(store, MemoryStore):
        root = zarr.open_group(store, mode="r", use_consolidated=False)
        if group_name not in root:
            return set()
        group = root[group_name]
        assert isinstance(group, zarr.Group)
        return set(group.keys())
    with h5py.File(store, "r") as f:
        return set(f[group_name].keys()) if group_name in f else set()


@pytest.mark.parametrize("validate", [False, True], ids=["validate=F", "validate=T"])
def test_use_object_column_converted_like_full_write(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    tmp_path: Path,
    validate: bool,
) -> None:
    """`use=A.obs[col]` writes an object column the way a full write does."""
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    write_and_read(adata, diskfmt, diskfmt_store)
    adata.obs["o"] = np.array([1, None, 3], dtype=object)
    full_store = tmp_path / "full.h5ad" if diskfmt == "h5ad" else MemoryStore()
    full = as_dataframe(write_and_read(adata, diskfmt, full_store).obs)["o"]
    with ad.settings.override(validate_before_write=validate):
        read = as_dataframe(
            write_and_read(adata, diskfmt, diskfmt_store, use=A.obs["o"]).obs
        )["o"]
    assert read.dtype == full.dtype == "Int64"
    assert read.isna().tolist() == [False, True, False]
    assert as_dataframe(adata.obs)["o"].dtype == object


@pytest.mark.parametrize("validate", [False, True], ids=["validate=F", "validate=T"])
def test_use_unconvertible_object_column_writes_nothing(
    *,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    validate: bool,
) -> None:
    adata = adata_with_obs(col=obs_column([1, 2, 3]))
    write_and_read(adata, diskfmt, diskfmt_store)
    adata.obsm["good"] = np.ones((3, 2))
    adata.obs["mixed"] = np.array([1, "a", None], dtype=object)
    with (
        ad.settings.override(validate_before_write=validate),
        pytest.raises(UnsupportedColumnTypeError, match=r"'mixed'"),
    ):
        write_and_read(
            adata, diskfmt, diskfmt_store, use=[A.obsm["good"], A.obs["mixed"]]
        )
    assert _keys_on_disk(diskfmt_store) == set()
    assert "mixed" not in _keys_on_disk(diskfmt_store, "obs")


@pytest.mark.parametrize("part", ["X", "layer"])
def test_use_array_parts(
    *, diskfmt: Literal["h5ad", "zarr"], diskfmt_store: Path | MemoryStore, part: str
) -> None:
    adata = ad.AnnData(
        np.zeros((3, 2)),
        obs=pd.DataFrame(index=obs_index()),
        var=pd.DataFrame(index=["x", "y"]),
    )
    adata.layers["L"] = np.zeros((3, 2))
    write_and_read(adata, diskfmt, diskfmt_store)
    use = A.X if part == "X" else A.layers["L"]
    adata.obsm["bad"] = pd.DataFrame(
        {"p": pd.period_range("2020-01", periods=3, freq="M")}, index=adata.obs_names
    )
    new = np.arange(6.0).reshape(3, 2)
    if part == "X":
        adata.X = new
    else:
        adata.layers["L"] = new
    with ad.settings.override(validate_before_write=True):
        # an unwritable part alongside: nothing is written
        with pytest.raises(IORegistryError):
            write_and_read(adata, diskfmt, diskfmt_store, use=[use, A.obsm["bad"]])
        read = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
        assert not read.X.any()
        assert not read.layers["L"].any()
        read = write_and_read(adata, diskfmt, diskfmt_store, use=use)
    np.testing.assert_array_equal(read.X if part == "X" else read.layers["L"], new)


@pytest.mark.parametrize("source", ["backed", "view"])
def test_use_from_backed_or_view(*, tmp_path: Path, source: str) -> None:
    """Parts of a backed object or a view are checked without reading them."""
    rng = np.random.default_rng(0)
    base = ad.AnnData(
        rng.random((4, 3)),
        obs=pd.DataFrame(index=list("abcd")),
        var=pd.DataFrame(index=list("xyz")),
        layers={"L": sp.random(4, 3, density=0.5, format="csr", rng=rng)},
        obsm={"m": rng.random((4, 2))},
    )
    period = pd.period_range("2020-01", periods=4, freq="M")
    subject: ad.AnnData
    if source == "backed":
        base.write_h5ad(tmp_path / "src.h5ad")
        subject = ad.read_h5ad(tmp_path / "src.h5ad", backed="r")
        subject.obsm["bad"] = pd.DataFrame({"p": period}, index=subject.obs_names)
    else:
        base.obsm["bad"] = pd.DataFrame({"p": period}, index=base.obs_names)
        subject = base[:2]
        assert subject.is_view
    dest = tmp_path / "dest.h5ad"
    blank = ad.AnnData(
        np.zeros(subject.shape),
        obs=pd.DataFrame(index=subject.obs_names),
        var=pd.DataFrame(index=subject.var_names),
    )
    blank.write_h5ad(dest)
    uses: list[AdRef | RefAcc | MapAcc] = [A.X, A.layers["L"], A.obsm["m"]]

    with h5py.File(dest, "r+") as f:
        parts = [_resolve_use(f, subject, use) for use in uses]
        with forbid_payload_access():
            raise_if_unwriteable_parts(parts, store=f)
    with ad.settings.override(validate_before_write=True):
        with pytest.raises(IORegistryError):
            subject.write_h5ad(dest, use=[*uses, A.obsm["bad"]])
        assert _keys_on_disk(dest) == set()
        assert not _dense(ad.read_h5ad(dest).X).any()
        subject.write_h5ad(dest, use=uses)
    read = ad.read_h5ad(dest)
    expected = subject.to_memory()
    np.testing.assert_array_equal(_dense(read.X), _dense(expected.X))
    np.testing.assert_array_equal(
        _dense(read.layers["L"]), _dense(expected.layers["L"])
    )
    np.testing.assert_array_equal(_dense(read.obsm["m"]), _dense(expected.obsm["m"]))


def _dense(x: object) -> np.ndarray:
    if isinstance(x, sp.csr_matrix | sp.csc_matrix | sp.csr_array | sp.csc_array):
        return x.toarray()
    assert isinstance(x, np.ndarray)
    return x


def test_unwriteable_leaves_closed_backed_file_closed(*, tmp_path: Path) -> None:
    """Stopping at the first problem still closes a backed file the walk reopened."""
    adata = ad.AnnData(np.ones((3, 2)), obs=pd.DataFrame(index=obs_index()))
    adata.write_h5ad(tmp_path / "a.h5ad")
    backed = ad.read_h5ad(tmp_path / "a.h5ad", backed="r")
    # `layers` holds `X`, so walking it reopens the file
    backed.layers["bad"] = np.zeros((3, 2)).view(_ArraySubclass)
    backed.file.close()
    assert backed.unwriteable()
    assert not backed.file.is_open
