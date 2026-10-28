"""\
Best-effort preflight check of an :class:`~anndata.AnnData` object before writing.

The walker in this module visits every element the writers would visit,
inspects its *type* (plus cheap metadata such as dtypes and shapes) and asks
the IO registry whether a writer exists, without writing anything.
It backs :meth:`anndata.AnnData.unwriteable` and, when
:attr:`anndata.settings.validate_before_write` is set, runs at the start of
:func:`~anndata.io.write_h5ad` and :func:`~anndata.io.write_zarr`.

What it catches: elements and dataframe columns with no registered writer,
object columns that cannot be converted, non-2D ``X``/layers, and dataframe
structure problems (reserved or invalid names, duplicate columns, index-name
collisions, unsupported index types).
What it does not catch: backend, codec and IO failures, and writers that
reject particular *values* rather than types (e.g. a categorical with
datetime categories on h5ad). Those can still fail mid-write.

What it never does: load data from a backed file, compute a dask array,
densify, subset a view's payload, or export an array-API array. The one
allocation it makes is the conversion of object-dtype dataframe columns,
which has to look at the values to find out whether they are writable.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import PurePosixPath
from types import MappingProxyType
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr
from pandas.api.types import is_string_dtype

from anndata._core.aligned_mapping import _on_disk_x
from anndata._core.anndata import AnnData
from anndata._core.raw import Raw
from anndata._core.storage import _non_2d_message
from anndata._settings import settings
from anndata.compat import XDataset, has_xp
from anndata.utils import iter_outer

from ..utils import check_key
from .dataframe_columns import add_column_note, normalize_column
from .registry import _REGISTRY, IORegistryError, _iter_patterns

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal

    from pandas.api.extensions import ExtensionArray

    from ..._types import _GroupStorageType


__all__ = ["iter_unwriteable", "raise_if_unwriteable"]


_STORE_TYPES: Mapping[str, tuple[type[_GroupStorageType], ...]] = MappingProxyType({
    "h5": (h5py.Group,),
    "zarr": (zarr.Group,),
})


def _dest_types(store_type: Literal["h5", "zarr"] | None) -> tuple[type, ...]:
    if store_type is None:
        return tuple(t for ts in _STORE_TYPES.values() for t in ts)
    return _STORE_TYPES[store_type]


def iter_unwriteable(
    elem: object,
    *,
    store_type: Literal["h5", "zarr"] | None = None,
    path: str = "/",
) -> Generator[tuple[str, Exception]]:
    """\
    Yield ``(path, exception)`` for every part of `elem` that cannot be written.

    `elem` is an :class:`~anndata.AnnData`, a :class:`~anndata.Raw`,
    a dataframe, a dict, or a single element.
    Paths are slash-separated like the keys in a file (``/obs/cluster``),
    and are also attached to each exception as a note.
    Each exception is what writing that part would (roughly) have raised.
    See the module docstring for what is and is not caught.

    Parameters
    ----------
    store_type
        Which backing store to check against; `None` means the element only
        needs a writer for one of them.
    """
    dests = _dest_types(store_type)
    yield from _walk(elem, dests, path)


def raise_if_unwriteable(
    adata: AnnData, *, store_type: Literal["h5", "zarr"] | None = None
) -> None:
    """\
    Raise if any part of `adata` cannot be written, before anything is written.

    A single problem raises the exception for that part
    (an :exc:`~anndata._io.specs.registry.IORegistryError` for a missing writer,
    a :exc:`TypeError` for an unconvertible column or an invalid key, a
    :exc:`ValueError` for a structural problem). Several problems are raised
    together as an :exc:`ExceptionGroup` whose message lists their paths and
    whose members keep their types; catch those with ``except*``.
    """
    errors = list(iter_unwriteable(adata, store_type=store_type))
    if not errors:
        return
    if len(errors) == 1:
        raise errors[0][1]
    paths = ", ".join(path for path, _ in errors)
    msg = f"{len(errors)} elements cannot be written: {paths}"
    raise ExceptionGroup(msg, [e for _, e in errors])


def _join(path: str, key: object) -> str:
    return f"{path.rstrip('/')}/{key}"


def _with_path[E: Exception](e: E, path: str) -> E:
    e.add_note(f"Element: {path}")
    return e


def _registered(elem: object, dests: tuple[type, ...]) -> bool:
    return any(
        _REGISTRY.has_write(dest, pattern, frozenset())
        for dest in dests
        for pattern in _iter_patterns(elem)
    )


def _registry_miss(elem: object, dests: tuple[type, ...], path: str) -> Exception:
    e = IORegistryError._from_write_parts(dests[0], type(elem), frozenset())
    return _with_path(e, path)


def _slash_forbidden(dests: tuple[type, ...]) -> bool:
    """Whether `write_elem` would refuse a key containing a slash for every destination."""
    return all(
        dest is zarr.Group or settings.disallow_forward_slash_in_h5ad for dest in dests
    )


def _walk(
    elem: object,
    dests: tuple[type, ...],
    path: str,
    ancestors: frozenset[int] = frozenset(),
) -> Generator[tuple[str, Exception]]:
    # The branches mirror the writers' dispatch: whatever `write_elem` would do
    # with `elem` decides what is inspected here. `ancestors` holds the ids of
    # the containers on the current branch, so a cycle is reported once and the
    # branch stopped (the writers would recurse until `RecursionError`).
    if elem is None:
        return
    if isinstance(elem, AnnData | Raw | XDataset | Mapping):
        if id(elem) in ancestors:
            msg = f"{type(elem).__name__} at {path!r} contains itself (reference cycle)"
            yield path, _with_path(ValueError(msg), path)
            return
        ancestors = ancestors | {id(elem)}
    if isinstance(elem, AnnData):
        yield from _walk_anndata(elem, dests, path, ancestors)
    elif isinstance(elem, Raw):
        yield from _walk_raw(elem, dests, path, ancestors)
    elif isinstance(elem, pd.DataFrame):
        yield from _walk_dataframe(elem, dests, path)
    elif isinstance(elem, pd.Categorical):
        yield from _walk_categorical(elem, dests, path)
    elif isinstance(elem, dict) or (
        isinstance(elem, XDataset | Mapping) and _registered(elem, dests)
    ):
        # `normalize_nested` turns any dict subclass into a plain dict, which
        # has a writer that writes each item; other mappings need their own
        for key, value in elem.items():
            key_path = _join(path, key)
            if (key_error := _check_mapping_key(key, path, dests)) is not None:
                yield key_path, _with_path(key_error, key_path)
                continue
            yield from _walk(value, dests, key_path, ancestors)
    else:
        yield from _walk_leaf(elem, dests, path)


def _check_mapping_key(
    key: object, path: str, dests: tuple[type, ...]
) -> Exception | None:
    """\
    Mirror what `Writer.write_elem` does with a key written into the group at `path`.

    See `registry.py`: a non-string key fails at `k.startswith("/")`; an absolute key
    is made relative to the group (`PurePosixPath(k).relative_to(store.name)`, which
    rejects keys outside it); a remaining slash is refused for zarr and, depending on
    `settings.disallow_forward_slash_in_h5ad`, for h5ad.
    """
    if not isinstance(key, str):
        msg = f"{key!r} of type {type(key)} is an invalid key. Should be str."
        return TypeError(msg)
    if key.startswith("/"):
        try:
            key = str(PurePosixPath(key).relative_to(path, walk_up=False))
        except ValueError:
            msg = f"{key!r} is not in the subpath of {path!r}"
            return ValueError(msg)
    if "/" in key and _slash_forbidden(dests):
        msg = f"Forward slashes are not allowed in keys in {dests[0]}"
        return ValueError(msg)
    return None


def _walk_leaf(
    elem: object, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # `to_writeable` turns DLPack-capable array-API arrays into numpy arrays at
    # write time; never export them here, just accept them
    if has_xp(elem) and not isinstance(elem, np.ndarray):
        return
    if not _registered(elem, dests):
        yield path, _registry_miss(elem, dests, path)


_ALIGNED_MAPPINGS = frozenset({"layers", "obsm", "varm", "obsp", "varp"})


def _walk_anndata(
    adata: AnnData, dests: tuple[type, ...], path: str, ancestors: frozenset[int]
) -> Generator[tuple[str, Exception]]:
    # For a view, indexing `X`, `layers`, `obsm`, … or `raw.X` subsets the parent's
    # payload, so those are inspected on the parent (same keys, same types; a
    # dataframe's columns are judged on the parent's values). The view's own
    # `obs`/`var` (materialised when the view was made) and `uns` are inspected
    # directly.
    payload_source = adata._adata_ref if adata.is_view else adata
    assert payload_source is not None
    for attr_name, attr in iter_outer(adata):
        attr_path = _join(path, attr_name)
        if attr_name == "uns":
            assert isinstance(attr, Mapping)
            # the writer copies `uns` into a plain dict; remember the original too
            yield from _walk(dict(attr), dests, attr_path, ancestors | {id(attr)})
        elif attr_name in _ALIGNED_MAPPINGS:
            stored: Mapping[str | None, object] = getattr(payload_source, attr_name)
            # the writers refuse slashes in these keys up front (`_check_has_no_slash_key`)
            if any(isinstance(k, str) and "/" in k for k in stored):
                slash_msg = f"Forward slashes are not allowed in keys in {attr_name}"
                yield attr_path, _with_path(ValueError(slash_msg), attr_path)
            for key, value in stored.items():
                if key is None:  # `.X`, stored under `layers[None]`
                    name, elem_path = "X", _join(path, "X")
                else:
                    name, elem_path = f"Layer {key!r}", _join(attr_path, key)
                    if not isinstance(key, str):  # slashes were reported above
                        key_error = _check_mapping_key(key, attr_path, dests)
                        assert key_error is not None
                        yield elem_path, _with_path(key_error, elem_path)
                        continue
                if attr_name == "layers" and (
                    (msg := _non_2d_message(value, name=name)) is not None
                ):
                    yield elem_path, _with_path(ValueError(msg), elem_path)
                    continue
                yield from _walk(value, dests, elem_path, ancestors)
        elif isinstance(attr, Raw):
            yield from _walk_raw(attr, dests, attr_path, ancestors)
        else:  # obs, var: a dataframe, or an unsupported container such as Dataset2D
            yield from _walk(attr, dests, attr_path, ancestors)


def _walk_raw(
    raw: Raw, dests: tuple[type, ...], path: str, ancestors: frozenset[int]
) -> Generator[tuple[str, Exception]]:
    # `raw.X` subsets a backed view's dataset; inspect the stored object instead
    X: object
    if raw._adata.isbacked:
        file = raw._adata.file
        if not file.is_open:
            file.open()
        key = next((k for k in ("raw/X", "raw.X") if k in file), None)
        X = None if key is None else _on_disk_x(file[key])
    else:
        X = raw._X
    yield from _walk(X, dests, _join(path, "X"), ancestors)
    yield from _walk(raw.var, dests, _join(path, "var"), ancestors)
    yield from _walk(dict(raw.varm), dests, _join(path, "varm"), ancestors)


def _walk_categorical(
    cat: pd.Categorical, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # The writer stores `categories.to_numpy()`; an empty slice has the same dtype
    # as the whole index, so neither the categories nor (for a `RangeIndex`) their
    # materialisation is ever touched.
    categories_path = _join(path, "categories")
    categories = cat.categories[:0]
    representative = categories.to_numpy()
    if representative.dtype == object and not is_string_dtype(categories.dtype):
        # e.g. intervals, periods, tz-aware datetimes: `to_numpy()` makes them
        # arbitrary Python objects, which the string writer would stringify (zarr)
        # or refuse (h5py); say so instead of letting either happen
        msg = (
            f"Categories of dtype {categories.dtype} are not written: "
            "zarr would store them as strings and h5py refuses them. "
            "Convert the categories first, e.g. with `.cat.rename_categories(str)`."
        )
        yield categories_path, _with_path(TypeError(msg), categories_path)
    else:
        yield from _walk_leaf(representative, dests, categories_path)


def _walk_dataframe(  # noqa: PLR0912, PLR0915
    df: pd.DataFrame, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # The structural checks mirror `write_dataframe`. Columns are reached by
    # position so that neither the column index nor the row index (possibly a
    # `RangeIndex`) is materialised, and no column is boxed into a Series.
    columns = df.columns
    if "_index" in columns:
        msg = "'_index' is a reserved name for dataframe columns."
        yield path, _with_path(ValueError(msg), path)
    if not columns.is_unique:
        duplicates = list(columns[columns.duplicated()])
        msg = f"Found repeated column names: {duplicates}. Column names must be unique."
        yield path, _with_path(ValueError(msg), path)
    index_name = df.index.name if df.index.name is not None else "_index"
    n_columns = df.shape[1]
    for key in (*(columns[i] for i in range(n_columns)), index_name):
        try:
            key = check_key(key)
        except TypeError as e:
            yield path, _with_path(e, path)
        else:
            if "/" in key and _slash_forbidden(dests):
                msg = f"Forward slashes are not allowed in keys in {dests[0]}"
                yield _join(path, key), _with_path(ValueError(msg), _join(path, key))
    is_multi = isinstance(df.index, pd.MultiIndex)
    if (
        not is_multi  # a MultiIndex cannot be compared this way; reported below
        and df.index.name is not None
        and df.index.name in columns
        and not pd.Series(df.index, index=df.index).equals(df[df.index.name])
    ):
        msg = (
            f"DataFrame.index.name ({df.index.name!r}) is also used by a column "
            "whose values are different. This is not supported. Please make sure "
            "the values are the same, or use a different name."
        )
        yield path, _with_path(ValueError(msg), path)
    index_path = _join(path, index_name)
    if is_multi:
        msg = "Cannot write a DataFrame with a MultiIndex. Use a single-level index."
        yield index_path, _with_path(ValueError(msg), index_path)
    else:
        # The writer stores `index.to_numpy()` for numpy-backed indexes and
        # `index.array` otherwise (e.g. a `Categorical`); an empty slice has the
        # same dtype and array type without materialising the index.
        empty_index = df.index[:0]
        index_values = (
            empty_index.to_numpy()
            if isinstance(empty_index.array, pd.arrays.NumpyExtensionArray)
            else empty_index.array
        )
        yield from _walk(index_values, dests, index_path)
    dtypes = list(df.dtypes)
    numpy_kind_registered: dict[str, bool] = {}
    for i in range(n_columns):
        colname, dtype = columns[i], dtypes[i]
        col_path = _join(path, colname)
        if isinstance(dtype, np.dtype) and dtype.kind not in "OMm":
            # A plain numpy column is written as is and its registry pattern depends
            # only on the dtype kind, so avoid boxing the column into a Series
            # (wide dataframes in `uns` can have tens of thousands of columns).
            if dtype.kind not in numpy_kind_registered:
                numpy_kind_registered[dtype.kind] = _registered(
                    np.empty(0, dtype=dtype), dests
                )
            if not numpy_kind_registered[dtype.kind]:
                stub = np.empty(0, dtype=dtype)
                miss = _registry_miss(stub, dests, col_path)
                add_column_note(miss, stub, col_path)
                yield col_path, miss
            continue
        try:
            values = normalize_column(_column_array(df, i), path=col_path)
        except TypeError as conversion_error:
            yield col_path, _with_path(conversion_error, col_path)
            continue
        if isinstance(values, pd.Categorical):
            yield from _walk_categorical(values, dests, col_path)
        elif not _registered(values, dests):
            miss = _registry_miss(values, dests, col_path)
            add_column_note(miss, values, col_path)
            yield col_path, miss


def _column_array(df: pd.DataFrame, i: int) -> np.ndarray | ExtensionArray:
    """The backing array of column `i`, like `series._values`, without boxing a Series."""
    return df._get_column_array(i)  # type: ignore[operator]
