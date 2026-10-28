"""\
Check whether an :class:`~anndata.AnnData` object can be written, without writing it.

The walk visits what the writers would visit and judges types, dtypes and shapes
with the writers' own checks, without touching payload data;
see :ref:`type-mapping` for what it catches.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr
from pandas.api.types import is_string_dtype

from anndata._core.anndata import AnnData
from anndata._core.raw import Raw
from anndata._core.storage import _non_2d_message
from anndata.compat import XDataset, has_xp
from anndata.utils import iter_outer

from ..utils import _check_has_no_slash_key, check_key
from .dataframe_columns import (
    UnsupportedColumnTypeError,
    add_column_note,
    normalize_column,
)
from .methods import (
    _check_dataframe_columns,
    _dataframe_index_name,
    _dataframe_index_values,
)
from .registry import _REGISTRY, IORegistryError, _iter_patterns, _normalize_write_key

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable
    from typing import Literal

    from ..._types import _GroupStorageType

    type Problems = Generator[tuple[str, Exception]]
    type GroupTypes = tuple[type[_GroupStorageType], ...]


__all__ = ["iter_unwriteable", "raise_if_unwriteable", "raise_if_unwriteable_parts"]


_GROUP_TYPES: Mapping[str, GroupTypes] = MappingProxyType({
    "h5": (h5py.Group,),
    "zarr": (zarr.Group,),
})
_ALIGNED_MAPPINGS = frozenset({"layers", "obsm", "varm", "obsp", "varp"})
_MAX_PATHS_IN_MESSAGE = 10


def iter_unwriteable(
    elem: object,
    *,
    store_type: Literal["h5", "zarr"] | None = None,
    path: str = "/",
) -> Problems:
    """\
    Yield ``(path, exception)`` for every part of `elem` that cannot be written.

    `store_type` `None` means a writer for either store suffices.
    """
    group_types = _group_types(store_type)
    if isinstance(elem, AnnData):
        # the write entry points take any AnnData; nested ones dispatch on their type
        yield from _walk_anndata(elem, group_types, path, frozenset({id(elem)}))
    else:
        yield from _walk(elem, group_types, path)


def raise_if_unwriteable(
    adata: AnnData, *, store_type: Literal["h5", "zarr"] | None = None
) -> None:
    """\
    Raise if any part of `adata` cannot be written, before anything is written.

    One problem raises its own exception: `IORegistryError` for a missing writer,
    `TypeError` for an unconvertible column or invalid key, `ValueError` for
    a structural problem (the writers raise `AttributeError` for a non-str key
    and `RecursionError` for a cycle). Several raise an `ExceptionGroup`.
    """
    _raise_problems(list(iter_unwriteable(adata, store_type=store_type)))


def raise_if_unwriteable_parts(
    parts: Iterable[tuple[str | None, str, object]], *, store: _GroupStorageType
) -> None:
    """\
    Like :func:`raise_if_unwriteable`, for the `(parent, key, element)` triples
    `write_use` writes into `store` (`parent` `None` is the root).
    """
    group_types = (zarr.Group,) if isinstance(store, zarr.Group) else (h5py.Group,)
    errors: list[tuple[str, Exception]] = []
    for parent, key, elem in parts:
        group_path = "/" if parent is None else f"/{parent}"
        try:
            stored = _stored_key(key, group_path, group_types)
        except (TypeError, ValueError) as e:
            path = _child_path(group_path, key)
            errors.append((path, _note_path(e, path)))
            continue
        errors.extend(_walk(elem, group_types, _child_path(group_path, stored)))
    _raise_problems(errors)


def _raise_problems(errors: list[tuple[str, Exception]]) -> None:
    if not errors:
        return
    if len(errors) == 1:
        raise errors[0][1]
    paths = ", ".join(path for path, _ in errors[:_MAX_PATHS_IN_MESSAGE])
    if (n_more := len(errors) - _MAX_PATHS_IN_MESSAGE) > 0:
        paths += f", and {n_more} more"
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


def _stored_key(key: object, path: str, group_types: GroupTypes) -> str:
    """Return `key` as `Writer.write_elem` would store it under `path`, or raise."""
    name: str = check_key(key)
    errors: list[ValueError] = []
    for group_type in group_types:  # with several, one accepting the key suffices
        try:
            return _normalize_write_key(name, path, group_type)
        except ValueError as e:
            errors.append(e)
    raise errors[0]


def _walk(
    elem: object,
    group_types: GroupTypes,
    path: str,
    ancestors: frozenset[int] = frozenset(),
) -> Problems:
    """Dispatch like `write_elem`; `ancestors` are the containers on this branch."""
    if elem is None:
        return
    # `normalize_nested` turns dict subclasses into dicts;
    # everything else needs a writer for its exact type, like `find_write_func`
    if not isinstance(elem, dict) and not _has_writer(elem, group_types):
        yield from _walk_leaf(elem, group_types, path)
        return
    if isinstance(elem, AnnData | Raw | XDataset | Mapping | list):
        if id(elem) in ancestors:
            msg = f"{type(elem).__name__} at {path!r} contains itself (reference cycle)"
            yield path, _note_path(ValueError(msg), path)
            return
        ancestors = ancestors | {id(elem)}
    if isinstance(elem, list):
        # the list writer stores the whole sequence; only cycles matter here
        yield from _sequence_cycles(elem, path, ancestors)
    elif isinstance(elem, AnnData):
        yield from _walk_anndata(elem, group_types, path, ancestors)
    elif isinstance(elem, Raw):
        yield from _walk_raw(elem, group_types, path, ancestors)
    elif isinstance(elem, pd.DataFrame):
        yield from _walk_dataframe(elem, group_types, path)
    elif isinstance(elem, pd.Categorical):
        yield from _walk_categorical(elem, group_types, path)
    elif isinstance(elem, XDataset | Mapping):
        yield from _walk_mapping(elem, group_types, path, ancestors)


def _walk_mapping(
    elem: XDataset | Mapping,
    group_types: GroupTypes,
    path: str,
    ancestors: frozenset[int],
) -> Problems:
    for key, value in elem.items():
        try:
            stored = _stored_key(key, path, group_types)
        except (TypeError, ValueError) as e:
            yield _child_path(path, key), _note_path(e, _child_path(path, key))
            continue
        yield from _walk(value, group_types, _child_path(path, stored), ancestors)


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
    try:
        _check_has_no_slash_key(attr_name, stored)
    except ValueError as e:
        yield path, _note_path(e, path)
    except TypeError:
        pass  # a non-str key, reported below
    for key, value in stored.items():
        if key is None and attr_name == "layers":
            # the AnnData writers store `layers[None]` as `X`
            name, elem_path = "X", _child_path(adata_path, "X")
        else:
            name, elem_path = f"Layer {key!r}", _child_path(path, key)
            try:
                check_key(key)
            except TypeError as e:
                yield elem_path, _note_path(e, elem_path)
                continue
        if attr_name == "layers" and (
            (shape_msg := _non_2d_message(value, name=name)) is not None
        ):
            yield elem_path, _note_path(ValueError(shape_msg), elem_path)
        elif (
            adata.is_view
            and isinstance(value, pd.DataFrame)
            and _has_writer(value, group_types)
        ):
            # the view writes a row (obsm/obsp) or column (varm/varp) selection
            rows = adata._oidx if attr_name.startswith("obs") else adata._vidx
            yield from _walk_dataframe(value, group_types, elem_path, rows=rows)
        else:
            yield from _walk(value, group_types, elem_path, ancestors)


def _walk_raw(
    raw: Raw, group_types: GroupTypes, path: str, ancestors: frozenset[int]
) -> Problems:
    X_path = _child_path(path, "X")
    X: object = raw._X
    if raw._adata.isbacked:
        try:
            # `raw.X` subsets a backed view's dataset
            X = raw._x_on_disk()
        except AttributeError as e:
            yield X_path, _note_path(e, X_path)
    yield from _walk(X, group_types, X_path, ancestors)
    yield from _walk(raw.var, group_types, _child_path(path, "var"), ancestors)
    yield from _walk(dict(raw.varm), group_types, _child_path(path, "varm"), ancestors)


def _walk_categorical(
    cat: pd.Categorical, group_types: GroupTypes, path: str
) -> Problems:
    categories_path = _child_path(path, "categories")
    # `write_categorical` stores `categories.to_numpy()`; an empty slice has its dtype
    categories = cat.categories[:0]
    stored = categories.to_numpy()
    if stored.dtype == object and not is_string_dtype(categories.dtype):
        # stricter than the writers, which fail (h5py) or store strings (zarr)
        msg = (
            f"Cannot write categories of dtype {categories.dtype}: "
            "h5ad cannot store them and zarr would store them as strings. "
            "Convert them first, e.g. with `.cat.rename_categories(str)`."
        )
        yield categories_path, _note_path(TypeError(msg), categories_path)
        return
    yield from _walk(stored, group_types, categories_path)


def _walk_dataframe(
    df: pd.DataFrame,
    group_types: GroupTypes,
    path: str,
    *,
    rows: object = None,
) -> Problems:
    """Check `df` like `write_dataframe`; `rows`, a view's selection, applies to object columns."""
    yield from _check_dataframe_names(df, group_types, path)
    yield from _walk_index(df, group_types, path)
    yield from _walk_columns(df, group_types, path, rows=rows)


def _check_dataframe_names(
    df: pd.DataFrame, group_types: GroupTypes, path: str
) -> Problems:
    try:
        names = _check_dataframe_columns(df)
        index_name = check_key(_dataframe_index_name(df))
    except (TypeError, ValueError) as e:
        yield path, _note_path(e, path)
        return
    for name in [*names, index_name]:
        try:
            _stored_key(name, path, group_types)
        except ValueError as e:
            yield _child_path(path, name), _note_path(e, _child_path(path, name))


def _walk_index(df: pd.DataFrame, group_types: GroupTypes, path: str) -> Problems:
    index_path = _child_path(path, _dataframe_index_name(df))
    try:
        # an empty slice has the dtype and array type of the whole index
        values = _dataframe_index_values(df.index[:0])
    except ValueError as e:  # a `MultiIndex` has no single backing array
        yield index_path, _note_path(e, index_path)
        return
    yield from _walk(values, group_types, index_path)


def _walk_columns(
    df: pd.DataFrame, group_types: GroupTypes, path: str, *, rows: object
) -> Problems:
    # wide frames repeat a few numpy dtypes, whose writer depends on the dtype only
    numpy_dtype_has_writer: dict[np.dtype, bool] = {}
    for i, colname in enumerate(df.columns):
        col_path = _child_path(path, colname)
        # `series._values` without boxing a Series like `df.items()`: with the cache,
        # 256k columns in 11 `uns` frames take 0.3 s instead of 1.4 s
        values = df._get_column_array(i)  # type: ignore[operator]  # stubs: a column
        if isinstance(values, np.ndarray) and values.dtype != object:
            if (has_writer := numpy_dtype_has_writer.get(values.dtype)) is None:
                has_writer = _has_writer(values, group_types)
                numpy_dtype_has_writer[values.dtype] = has_writer
        else:
            if rows is not None and isinstance(values, np.ndarray):
                values = values[rows if isinstance(rows, slice) else np.asarray(rows)]
            try:
                values = normalize_column(values, path=col_path)
            except UnsupportedColumnTypeError as e:
                yield col_path, _note_path(e, col_path)
                continue
            has_writer = _has_writer(values, group_types)
        if not has_writer:
            missing = _missing_writer_error(values, group_types, col_path)
            add_column_note(missing, values, col_path)
            yield col_path, missing
        elif isinstance(values, pd.Categorical):
            yield from _walk_categorical(values, group_types, col_path)
