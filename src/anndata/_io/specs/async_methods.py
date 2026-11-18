"""Coroutine readers for zarr, see :func:`anndata.experimental.read_elem_async`.

An encoding without a reader here is an error, see :class:`AsyncReader`.
"""

from __future__ import annotations

import asyncio
from functools import partial
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from zarr import AsyncArray, AsyncGroup

from anndata import AnnData
from anndata._core.sparse_dataset import (
    _check_array,
    _construct,
    _format_from_attrs,
    _shape_from_attrs,
)
from anndata._types import ANNDATA_ELEMS

from .methods import (
    _accessor_from_json,
    _awkward_from_parts,
    _categorical_from_parts,
    _dataframe_from_parts,
    _dataframe_keys,
    _nullable_from_parts,
    _nullable_string_from_parts,
    _recarray_from_raw,
)
from .registry import _ASYNC_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Coroutine, Iterable, Sequence
    from typing import Any, Literal

    from numpy import typing as npt

    from anndata import acc
    from anndata._core.sparse_dataset import DenseType
    from anndata.compat import AwkArray, CSArray, CSMatrix, CupyCSMatrix

    from .methods import _BaseMaskedArray
    from .registry import AsyncReader


async def _gather[T](coros: Iterable[Coroutine[Any, Any, T]]) -> list[T]:
    """Like :func:`asyncio.gather`, but cancel the reads that are still pending when it ends."""
    tasks = [asyncio.ensure_future(coro) for coro in coros]
    try:
        return list(await asyncio.gather(*tasks))
    finally:
        if pending := [task for task in tasks if not task.done()]:
            for task in pending:
                task.cancel()
            await asyncio.gather(*pending, return_exceptions=True)


_MISSING = object()


async def _read_children(
    group: AsyncGroup,
    _reader: AsyncReader,
    keys: Sequence[str] | None = None,
    *,
    missing_ok: bool = False,
) -> dict[str, Any]:
    """Read the children `keys` of `group`, or all of them, each as soon as it is found."""

    async def read(k: str, node: AsyncArray | AsyncGroup | None = None) -> object:
        if node is None:
            try:
                node = await group.getitem(k)
            except KeyError:
                if not missing_ok:
                    raise
                return _MISSING
        return await _reader.read_elem(node)

    if keys is None:
        members = [m async for m in group.members()]
        keys = [k for k, _ in members]
        values = await _gather(read(k, node) for k, node in members)
    else:
        # Like the sync `elem[k]`: look up only these keys, without listing the group
        values = await _gather(map(read, keys))
    return {k: v for k, v in zip(keys, values, strict=True) if v is not _MISSING}


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("array", "0.2.0"))
@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("string-array", "0.2.0"))
async def read_array(elem: AsyncArray, *, _reader: AsyncReader) -> np.ndarray:
    return np.asarray(await elem.getitem(()))


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("numeric-scalar", "0.2.0"))
@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("bytes", "0.2.0"))
async def read_scalar(elem: AsyncArray, *, _reader: AsyncReader) -> np.number:
    return np.asarray(await elem.getitem(())).item()


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("string", "0.2.0"))
async def read_string(elem: AsyncArray, *, _reader: AsyncReader) -> str:
    return str(await elem.getitem(()))


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("rec-array", "0.2.0"))
async def read_recarray(
    elem: AsyncArray, *, _reader: AsyncReader
) -> np.recarray | npt.NDArray:
    return _recarray_from_raw(np.asarray(await elem.getitem(())))


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("accessor", "0.1.0"))
async def read_accessor(
    elem: AsyncArray, *, _reader: AsyncReader
) -> acc.AdRef | acc.LayerAcc | acc.MultiAcc | acc.GraphAcc:
    return _accessor_from_json(str(await elem.getitem(())), elem.attrs)


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("null", "0.1.0"))
async def read_null(elem: AsyncArray, *, _reader: AsyncReader) -> None:
    return None


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("dict", "0.1.0"))
async def read_mapping(elem: AsyncGroup, *, _reader: AsyncReader) -> dict[str, Any]:
    return await _read_children(elem, _reader)


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("awkward-array", "0.1.0"))
async def read_awkward(elem: AsyncGroup, *, _reader: AsyncReader) -> AwkArray:
    return _awkward_from_parts(await _read_children(elem, _reader), elem.attrs)


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("anndata", "0.1.0"))
@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("raw", "0.1.0"))
async def read_anndata(elem: AsyncGroup, *, _reader: AsyncReader) -> AnnData:
    return AnnData(
        **await _read_children(elem, _reader, ANNDATA_ELEMS, missing_ok=True)
    )


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("csr_matrix", "0.1.0"))
@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("csc_matrix", "0.1.0"))
async def read_sparse(
    elem: AsyncGroup, *, _reader: AsyncReader
) -> CSMatrix | CSArray | CupyCSMatrix:
    # Like the sync reader, a legacy `h5sparse_format` wins over the encoding type
    if (format := _format_from_attrs(elem.attrs, elem)) not in {"csr", "csc"}:
        msg = f"Unknown encoding type {format}"
        raise ValueError(msg)

    async def read_dense(k: str) -> DenseType:
        _check_array(arr := await elem.getitem(k), elem, k, AsyncArray)
        return cast("DenseType", await cast("AsyncArray", arr).getitem(...))

    data, indices, indptr = await _gather(
        map(read_dense, ["data", "indices", "indptr"])
    )
    return _construct(
        cast("Literal['csr', 'csc']", format),
        (data, indices, indptr),
        shape=_shape_from_attrs(elem.attrs),
        dtype=data.dtype,
    )


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("dataframe", "0.2.0"))
async def read_dataframe(elem: AsyncGroup, *, _reader: AsyncReader) -> pd.DataFrame:
    columns, idx_key = _dataframe_keys(elem.attrs)
    parts = await _read_children(elem, _reader, [*columns, idx_key])
    return _dataframe_from_parts(
        {k: parts[k] for k in columns}, parts[idx_key], columns, idx_key
    )


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("categorical", "0.2.0"))
async def read_categorical(elem: AsyncGroup, *, _reader: AsyncReader) -> pd.Categorical:
    parts = await _read_children(elem, _reader, ["codes", "categories"])
    return _categorical_from_parts(parts["codes"], parts["categories"], elem.attrs)


async def _read_nullable(
    elem: AsyncGroup, *, _reader: AsyncReader, array_type: _BaseMaskedArray
) -> pd.api.extensions.ExtensionArray:
    parts = await _read_children(elem, _reader, ["values", "mask"])
    return _nullable_from_parts(parts["values"], parts["mask"], array_type)


_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("nullable-integer", "0.1.0"))(
    partial(_read_nullable, array_type=pd.arrays.IntegerArray)
)
_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("nullable-boolean", "0.1.0"))(
    partial(_read_nullable, array_type=pd.arrays.BooleanArray)
)


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("nullable-string-array", "0.1.0"))
async def read_nullable_string(
    elem: AsyncGroup, *, _reader: AsyncReader
) -> pd.api.extensions.ExtensionArray:
    parts = await _read_children(elem, _reader, ["values", "mask"])
    return _nullable_string_from_parts(parts["values"], parts["mask"], elem.attrs)
