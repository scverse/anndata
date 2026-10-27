"""\
Best-effort preflight check of an :class:`~anndata.AnnData` object before writing.

The walker in this module inspects the *types* of every element the writers
would visit and asks the IO registry whether a writer exists, without writing
and without loading, computing, or copying any data. It backs
:meth:`anndata.AnnData.unwriteable` and, when
:attr:`anndata.settings.validate_before_write` is set, runs at the start of
:func:`~anndata.io.write_h5ad` and :func:`~anndata.io.write_zarr`.

It catches unsupported element and column types and dataframe structure
problems before anything is written. Backend, codec and IO failures, and
writers that reject particular *values* (not types), can still fail mid-write.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import zarr

from anndata._core.anndata import AnnData
from anndata._core.raw import Raw
from anndata._core.storage import _non_2d_message
from anndata.compat import AwkArray, XDataset
from anndata.types import SupportsArrayApiBase
from anndata.utils import iter_outer

from .dataframe_columns import add_column_note, normalize_column
from .registry import _REGISTRY, IORegistryError, _iter_patterns

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal

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
    a dataframe, a mapping, or a single element.
    Paths are slash-separated like the keys in a file (``/obs/cluster``),
    and are also attached to each exception as a note.
    Each exception is what writing that part would (roughly) have raised.

    Only types and cheap metadata (dtypes, shapes, dask metadata) are inspected:
    nothing is loaded from disk, computed, densified, subset, or copied.

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
    a :exc:`TypeError` for an unconvertible column, a :exc:`ValueError` for a
    structural problem). Several problems are raised together as an
    :exc:`ExceptionGroup` whose message lists their paths; catch those with
    ``except*``.
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


def _registered(elem: object, dests: tuple[type, ...]) -> bool:
    return any(
        _REGISTRY.has_write(dest, pattern, frozenset())
        for dest in dests
        for pattern in _iter_patterns(elem)
    )


def _registry_miss(elem: object, dests: tuple[type, ...], path: str) -> Exception:
    e = IORegistryError._from_write_parts(dests[0], type(elem), frozenset())
    e.add_note(f"Element: {path}")
    return e


def _walk(  # noqa: PLR0911
    elem: object, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    if elem is None:
        return
    if isinstance(elem, AnnData):
        yield from _walk_anndata(elem, dests, path)
        return
    if isinstance(elem, Raw):
        yield from _walk(elem.X, dests, _join(path, "X"))
        yield from _walk(elem.var, dests, _join(path, "var"))
        yield from _walk(dict(elem.varm), dests, _join(path, "varm"))
        return
    if isinstance(elem, pd.DataFrame):
        yield from _walk_dataframe(elem, dests, path)
        return
    if isinstance(elem, pd.Categorical):
        # the categories are written as their own array
        yield from _walk_leaf(elem.categories.array, dests, _join(path, "categories"))
        return
    if isinstance(elem, pd.Series | pd.Index):
        # matches behavior in methods.py
        yield from _walk_leaf(elem.array, dests, path)
        return
    if isinstance(elem, AwkArray):
        import awkward as ak

        container = ak.to_buffers(ak.to_packed(elem))
        for key, buffer in container[2].items():
            yield from _walk_leaf(buffer, dests, _join(path, key))
        return
    if isinstance(elem, XDataset | Mapping):
        for key, value in elem.items():
            yield from _walk(value, dests, _join(path, key))
        return
    yield from _walk_leaf(elem, dests, path)


def _walk_leaf(
    elem: object, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # array-API arrays are converted via DLPack at write time; never export them here
    if isinstance(elem, SupportsArrayApiBase) and not isinstance(elem, np.ndarray):
        return
    if isinstance(elem, pd.arrays.NumpyExtensionArray):
        elem = elem.to_numpy()
    if not _registered(elem, dests):
        yield path, _registry_miss(elem, dests, path)


def _walk_anndata(
    adata: AnnData, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # Subsetting a view materialises its elements; the parent holds the same types.
    src = adata._adata_ref if adata.is_view else adata
    for attr_name, attr in iter_outer(src):
        attr_path = _join(path, attr_name)
        if attr_name == "layers":
            assert isinstance(attr, Mapping)
            for key, value in attr.items():
                name, elem_path = (
                    ("X", _join(path, "X"))
                    if key is None
                    else (f"Layer {key!r}", _join(attr_path, key))
                )
                if (msg := _non_2d_message(value, name=name)) is not None:
                    e = ValueError(msg)
                    e.add_note(f"Element: {elem_path}")
                    yield elem_path, e
                else:
                    yield from _walk(value, dests, elem_path)
            continue
        yield from _walk(
            dict(attr) if isinstance(attr, Mapping) else attr, dests, attr_path
        )


def _walk_dataframe(
    df: pd.DataFrame, dests: tuple[type, ...], path: str
) -> Generator[tuple[str, Exception]]:
    # the structural checks mirror `write_dataframe`
    if "_index" in df.columns:
        msg = "'_index' is a reserved name for dataframe columns."
        yield path, _structural_error(msg, path)
    if not df.columns.is_unique:
        duplicates = list(df.columns[df.columns.duplicated()])
        msg = f"Found repeated column names: {duplicates}. Column names must be unique."
        yield path, _structural_error(msg, path)
    if (
        df.index.name is not None
        and df.index.name in df.columns
        and not pd.Series(df.index, index=df.index).equals(df[df.index.name])
    ):
        msg = (
            f"DataFrame.index.name ({df.index.name!r}) is also used by a column "
            "whose values are different. This is not supported. Please make sure "
            "the values are the same, or use a different name."
        )
        yield path, _structural_error(msg, path)
    index_path = _join(path, df.index.name if df.index.name is not None else "_index")
    yield from _walk_leaf(df.index.array, dests, index_path)
    numpy_kind_registered: dict[str, bool] = {}
    for i, (colname, dtype) in enumerate(df.dtypes.items()):
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
            values = normalize_column(df.iloc[:, i]._values, path=col_path)
        except TypeError as conversion_error:
            conversion_error.add_note(f"Element: {col_path}")
            yield col_path, conversion_error
            continue
        if isinstance(values, pd.Categorical):
            yield from _walk(values, dests, col_path)
        elif not _registered(values, dests):
            miss = _registry_miss(values, dests, col_path)
            add_column_note(miss, values, col_path)
            yield col_path, miss


def _structural_error(msg: str, path: str) -> ValueError:
    e = ValueError(msg)
    e.add_note(f"Element: {path}")
    return e
