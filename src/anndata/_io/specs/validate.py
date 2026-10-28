"""\
Preflight check of an :class:`~anndata.AnnData` object before writing.

The walk visits every element the writers would visit, inspects its type and
cheap metadata (dtypes, shapes), and asks the IO registry for a writer.
It backs :meth:`anndata.AnnData.unwriteable` and, with
:attr:`anndata.settings.validate_before_write`, runs at the start of
:func:`~anndata.io.write_h5ad` and :func:`~anndata.io.write_zarr`.

Caught: missing writers, unconvertible object columns, refused categories,
non-2D ``X``/layers, invalid or duplicate names, unsupported index types.
Not caught: backend, codec and IO failures, and writers that reject values
rather than types.

The walk never materialises data: no reads from backed files, no dask
computation, no densifying, no subsetting of a view's payload, no dataframe
copies, no index materialisation, no awkward packing, no DLPack export.
Only object-dtype columns are scanned (for a view, along its selection),
an allocation bounded by the column. Whether a column named like the index
holds the index's values is therefore left to the writer.
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

from anndata._core.aligned_mapping import _on_disk_x
from anndata._core.anndata import AnnData
from anndata._core.raw import Raw
from anndata._core.storage import _non_2d_message
from anndata._settings import settings
from anndata.compat import XDataset, has_xp
from anndata.utils import iter_outer

from ..utils import check_key
from .dataframe_columns import (
    add_column_note,
    check_categories_writable,
    normalize_column,
)
from .registry import _REGISTRY, IORegistryError, _iter_patterns

if TYPE_CHECKING:
    from collections.abc import Generator, Hashable
    from typing import Literal

    from pandas.api.extensions import ExtensionArray

    from ..._types import _GroupStorageType

    type Problems = Generator[tuple[str, Exception]]
    type GroupTypes = tuple[type[_GroupStorageType], ...]


__all__ = ["iter_unwriteable", "raise_if_unwriteable"]


_GROUP_TYPES: Mapping[str, GroupTypes] = MappingProxyType({
    "h5": (h5py.Group,),
    "zarr": (zarr.Group,),
})
_ALIGNED_MAPPINGS = frozenset({"layers", "obsm", "varm", "obsp", "varp"})


def iter_unwriteable(
    elem: object,
    *,
    store_type: Literal["h5", "zarr"] | None = None,
    path: str = "/",
) -> Problems:
    """\
    Yield ``(path, exception)`` for every part of `elem` that cannot be written.

    Parameters
    ----------
    elem
        An :class:`~anndata.AnnData`, :class:`~anndata.Raw`, dataframe, dict,
        or single element.
    store_type
        Which backing store to check against;
        `None` means a writer for either store suffices.
    path
        Slash-separated location of `elem`, used in the yielded paths
        and attached to each exception as a note.
    """
    yield from _walk(elem, _group_types(store_type), path)


def raise_if_unwriteable(
    adata: AnnData, *, store_type: Literal["h5", "zarr"] | None = None
) -> None:
    """\
    Raise if any part of `adata` cannot be written, before anything is written.

    One problem raises its own exception
    (:exc:`~anndata._io.specs.registry.IORegistryError` for a missing writer,
    :exc:`TypeError` for an unconvertible column or invalid key,
    :exc:`ValueError` for a structural problem).
    Several problems raise an :exc:`ExceptionGroup` listing their paths;
    catch those with ``except*``.
    """
    errors = list(iter_unwriteable(adata, store_type=store_type))
    if not errors:
        return
    if len(errors) == 1:
        raise errors[0][1]
    paths = ", ".join(path for path, _ in errors)
    msg = f"{len(errors)} elements cannot be written: {paths}"
    raise ExceptionGroup(msg, [e for _, e in errors])


def _group_types(store_type: Literal["h5", "zarr"] | None) -> GroupTypes:
    if store_type is None:
        return tuple(t for ts in _GROUP_TYPES.values() for t in ts)
    return _GROUP_TYPES[store_type]


def _child_path(path: str, key: object) -> str:
    return f"{path.rstrip('/')}/{key}"


def _note_path[E: Exception](e: E, path: str) -> E:
    e.add_note(f"Element: {path}")
    return e


def _has_writer(elem: object, group_types: GroupTypes) -> bool:
    return any(
        _REGISTRY.has_write(group_type, pattern, frozenset())
        for group_type in group_types
        for pattern in _iter_patterns(elem)
    )


def _missing_writer_error(
    elem: object, group_types: GroupTypes, path: str
) -> IORegistryError:
    e = IORegistryError._from_write_parts(group_types[0], type(elem), frozenset())
    return _note_path(e, path)


def _slash_forbidden(group_types: GroupTypes) -> bool:
    return all(
        group_type is zarr.Group or settings.disallow_forward_slash_in_h5ad
        for group_type in group_types
    )


def _stored_key(key: object, path: str, group_types: GroupTypes) -> str:
    """Return `key` as `Writer.write_elem` would store it under `path`, or raise."""
    if not isinstance(key, str):
        msg = f"{key!r} of type {type(key)} is an invalid key. Should be str."
        raise TypeError(msg)
    if key.startswith("/"):
        try:
            key = str(PurePosixPath(key).relative_to(path, walk_up=False))
        except ValueError:
            msg = f"{key!r} is not in the subpath of {path!r}"
            raise ValueError(msg) from None
    if "/" in key and _slash_forbidden(group_types):
        msg = f"Forward slashes are not allowed in keys in {group_types[0]}"
        raise ValueError(msg)
    return key


def _walk(
    elem: object,
    group_types: GroupTypes,
    path: str,
    ancestors: frozenset[int] = frozenset(),
) -> Problems:
    """Dispatch like `write_elem`; `ancestors` are the containers on this branch."""
    if elem is None:
        return
    if isinstance(elem, AnnData | Raw | XDataset | Mapping | list | tuple):
        if id(elem) in ancestors:
            msg = f"{type(elem).__name__} at {path!r} contains itself (reference cycle)"
            yield path, _note_path(ValueError(msg), path)
            return
        ancestors = ancestors | {id(elem)}
    if isinstance(elem, list | tuple):
        # the list writer stores the whole sequence; only cycles matter here
        yield from _sequence_cycles(elem, path, ancestors)
    if isinstance(elem, AnnData):
        yield from _walk_anndata(elem, group_types, path, ancestors)
    elif isinstance(elem, Raw):
        yield from _walk_raw(elem, group_types, path, ancestors)
    elif isinstance(elem, pd.DataFrame):
        yield from _walk_dataframe(elem, group_types, path)
    elif isinstance(elem, pd.Categorical):
        yield from _walk_categorical(elem, group_types, path)
    elif isinstance(elem, dict) or (
        isinstance(elem, XDataset | Mapping) and _has_writer(elem, group_types)
    ):
        # `normalize_nested` turns dict subclasses into dicts
        for key, value in elem.items():
            try:
                stored = _stored_key(key, path, group_types)
            except (TypeError, ValueError) as e:
                yield _child_path(path, key), _note_path(e, _child_path(path, key))
                continue
            yield from _walk(value, group_types, _child_path(path, stored), ancestors)
    else:
        yield from _walk_leaf(elem, group_types, path)


def _sequence_cycles(
    seq: list | tuple | Mapping, path: str, ancestors: frozenset[int]
) -> Problems:
    """Report containers inside `seq` that are already on the current branch."""
    members = seq.items() if isinstance(seq, Mapping) else enumerate(seq)
    for key, item in members:
        if not isinstance(item, Mapping | list | tuple):
            continue
        item_path = _child_path(path, key)
        if id(item) in ancestors:
            msg = f"{type(item).__name__} at {item_path!r} contains itself (reference cycle)"
            yield item_path, _note_path(ValueError(msg), item_path)
            continue
        yield from _sequence_cycles(item, item_path, ancestors | {id(item)})


def _walk_leaf(elem: object, group_types: GroupTypes, path: str) -> Problems:
    # `to_writeable` converts these via DLPack at write time
    if has_xp(elem) and not isinstance(elem, np.ndarray):
        return
    if not _has_writer(elem, group_types):
        yield path, _missing_writer_error(elem, group_types, path)


def _walk_anndata(
    adata: AnnData, group_types: GroupTypes, path: str, ancestors: frozenset[int]
) -> Problems:
    # indexing a view's X/layers/obsm/… subsets the parent's payload
    payload_source = adata._adata_ref if adata.is_view else adata
    assert payload_source is not None
    for attr_name, attr in iter_outer(adata):
        attr_path = _child_path(path, attr_name)
        if attr_name == "uns":
            assert isinstance(attr, Mapping)
            yield from _walk(dict(attr), group_types, attr_path, ancestors | {id(attr)})
        elif attr_name in _ALIGNED_MAPPINGS:
            stored = getattr(payload_source, attr_name)
            yield from _walk_aligned_mapping(
                adata,
                attr_name,
                stored=stored,
                group_types=group_types,
                adata_path=path,
                ancestors=ancestors,
            )
        elif isinstance(attr, Raw):
            yield from _walk_raw(attr, group_types, attr_path, ancestors)
        else:
            yield from _walk(attr, group_types, attr_path, ancestors)


def _walk_aligned_mapping(
    adata: AnnData,
    attr_name: str,
    *,
    stored: Mapping[str | None, object],
    group_types: GroupTypes,
    adata_path: str,
    ancestors: frozenset[int],
) -> Problems:
    path = _child_path(adata_path, attr_name)
    if any(isinstance(k, str) and "/" in k for k in stored):
        msg = f"Forward slashes are not allowed in keys in {attr_name}"
        yield path, _note_path(ValueError(msg), path)
    for key, value in stored.items():
        if key is None and attr_name == "layers":
            # the AnnData writers store `layers[None]` as `X`
            name, elem_path = "X", _child_path(adata_path, "X")
        else:
            name, elem_path = f"Layer {key!r}", _child_path(path, key)
            if not isinstance(key, str):
                try:
                    _stored_key(key, path, group_types)
                except TypeError as e:
                    yield elem_path, _note_path(e, elem_path)
                    continue
        if attr_name == "layers" and (
            (shape_msg := _non_2d_message(value, name=name)) is not None
        ):
            yield elem_path, _note_path(ValueError(shape_msg), elem_path)
        elif isinstance(value, pd.DataFrame) and adata.is_view:
            # the view writes a row (obsm/obsp) or column (varm/varp) selection
            rows = adata._oidx if attr_name.startswith("obs") else adata._vidx
            yield from _walk_dataframe(value, group_types, elem_path, rows=rows)
        else:
            yield from _walk(value, group_types, elem_path, ancestors)


def _walk_raw(
    raw: Raw, group_types: GroupTypes, path: str, ancestors: frozenset[int]
) -> Problems:
    # `raw.X` subsets a backed view's dataset
    X: object
    if raw._adata.isbacked:
        file = raw._adata.file
        if not file.is_open:
            file.open()
        key = next((k for k in ("raw/X", "raw.X") if k in file), None)
        X = None if key is None else _on_disk_x(file[key])
    else:
        X = raw._X
    yield from _walk(X, group_types, _child_path(path, "X"), ancestors)
    yield from _walk(raw.var, group_types, _child_path(path, "var"), ancestors)
    yield from _walk(dict(raw.varm), group_types, _child_path(path, "varm"), ancestors)


def _walk_categorical(
    cat: pd.Categorical, group_types: GroupTypes, path: str
) -> Problems:
    categories_path = _child_path(path, "categories")
    # an empty slice has the dtype of the whole index
    categories = cat.categories[:0]
    try:
        check_categories_writable(categories, path=categories_path)
    except TypeError as e:
        yield categories_path, _note_path(e, categories_path)
        return
    yield from _walk_leaf(categories.to_numpy(), group_types, categories_path)


def _walk_dataframe(
    df: pd.DataFrame,
    group_types: GroupTypes,
    path: str,
    *,
    rows: object = None,
) -> Problems:
    """Check a dataframe like `write_dataframe` would, without touching its data.

    `rows` is a view's selection along the dataframe's axis; it is applied to
    object columns only.
    """
    yield from _check_dataframe_names(df, group_types, path)
    yield from _walk_index(df, group_types, path)
    yield from _walk_columns(df, group_types, path, rows=rows)


def _check_dataframe_names(
    df: pd.DataFrame, group_types: GroupTypes, path: str
) -> Problems:
    columns = df.columns
    if "_index" in columns:
        msg = "'_index' is a reserved name for dataframe columns."
        yield path, _note_path(ValueError(msg), path)
    if not columns.is_unique:
        duplicates = list(columns[columns.duplicated()])
        msg = f"Found repeated column names: {duplicates}. Column names must be unique."
        yield path, _note_path(ValueError(msg), path)
    # positional access keeps a `RangeIndex` of columns unmaterialised
    names: list[Hashable] = [columns[i] for i in range(df.shape[1])]
    names.append(_index_name(df))
    for name in names:
        try:
            name = check_key(name)
        except TypeError as e:
            yield path, _note_path(e, path)
            continue
        if "/" in name and _slash_forbidden(group_types):
            msg = f"Forward slashes are not allowed in keys in {group_types[0]}"
            yield (
                _child_path(path, name),
                _note_path(ValueError(msg), _child_path(path, name)),
            )


def _index_name(df: pd.DataFrame) -> Hashable:
    return df.index.name if df.index.name is not None else "_index"


def _walk_index(df: pd.DataFrame, group_types: GroupTypes, path: str) -> Problems:
    index_path = _child_path(path, _index_name(df))
    if isinstance(df.index, pd.MultiIndex):
        msg = "Cannot write a DataFrame with a MultiIndex. Use a single-level index."
        yield index_path, _note_path(ValueError(msg), index_path)
        return
    # an empty slice has the dtype and array type of the whole index
    index = df.index[:0]
    index_values = (
        index.to_numpy()
        if isinstance(index.array, pd.arrays.NumpyExtensionArray)
        else index.array
    )
    yield from _walk(index_values, group_types, index_path)


def _walk_columns(
    df: pd.DataFrame, group_types: GroupTypes, path: str, *, rows: object
) -> Problems:
    columns = df.columns
    dtypes = list(df.dtypes)
    numpy_kind_has_writer: dict[str, bool] = {}
    for i in range(df.shape[1]):
        col_path = _child_path(path, columns[i])
        dtype = dtypes[i]
        if isinstance(dtype, np.dtype) and dtype.kind not in "OMm":
            # the writer pattern of a numpy column depends on the dtype kind only
            stub = np.empty(0, dtype=dtype)
            if dtype.kind not in numpy_kind_has_writer:
                numpy_kind_has_writer[dtype.kind] = _has_writer(stub, group_types)
            if not numpy_kind_has_writer[dtype.kind]:
                missing = _missing_writer_error(stub, group_types, col_path)
                add_column_note(missing, stub, col_path)
                yield col_path, missing
            continue
        values = _column_values(df, i)
        if rows is not None and isinstance(values, np.ndarray) and dtype.kind == "O":
            values = values[rows if isinstance(rows, slice) else np.asarray(rows)]
        try:
            values = normalize_column(values, path=col_path)
        except TypeError as e:
            yield col_path, _note_path(e, col_path)
            continue
        if isinstance(values, pd.Categorical):
            yield from _walk_categorical(values, group_types, col_path)
        elif not _has_writer(values, group_types):
            missing = _missing_writer_error(values, group_types, col_path)
            add_column_note(missing, values, col_path)
            yield col_path, missing


def _column_values(df: pd.DataFrame, i: int) -> np.ndarray | ExtensionArray:
    """Backing array of column `i`, like `series._values`, without a Series."""
    return df._get_column_array(i)  # type: ignore[operator]
