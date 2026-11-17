"""Coroutine readers for zarr, see :func:`anndata.experimental.read_elem_async`.

Encodings without a reader here fall back to the synchronous reader, see :class:`AsyncReader`.
"""

from __future__ import annotations

import asyncio
from functools import partial
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from zarr import AsyncArray, AsyncGroup

from anndata import AnnData
from anndata._core.sparse_dataset import _construct, _shape_from_attrs
from anndata._types import AnnDataElem
from anndata.compat import _read_attr
from anndata.utils import get_literal_members

from .methods import (
    _dataframe_from_parts,
    _dataframe_keys,
    _nullable_string_from_parts,
)
from .registry import _ASYNC_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any, Literal

    from anndata._core.sparse_dataset import DenseType
    from anndata.compat import CSArray, CSMatrix, CupyCSMatrix
    from anndata.typing import RWAble

    from .methods import _BaseMaskedArray
    from .registry import AsyncReader


async def _read_children(
    group: AsyncGroup, keys: Iterable[str], _reader: AsyncReader
) -> list[Any]:
    async def read_child(k: str) -> RWAble:
        return await _reader.read_elem(await group.getitem(k))

    return await asyncio.gather(*map(read_child, keys))


async def _read_members(
    group: AsyncGroup, _reader: AsyncReader, keys: Iterable[str] | None = None
) -> dict[str, Any]:
    """Read all children of `group`, or the ones of `keys` that exist."""
    if keys is None:
        members = {k: v async for k, v in group.members()}
    else:
        # Like the sync `k in elem`: look up only these keys, without listing the group
        async def get(k: str) -> AsyncArray | AsyncGroup | None:
            try:
                return await group.getitem(k)
            except KeyError:
                return None

        nodes = await asyncio.gather(*map(get, keys))
        members = {k: v for k, v in zip(keys, nodes, strict=True) if v is not None}
    values = await asyncio.gather(*map(_reader.read_elem, members.values()))
    return dict(zip(members, values, strict=True))


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


@_ASYNC_REGISTRY.register_read(AsyncArray, IOSpec("null", "0.1.0"))
async def read_null(elem: AsyncArray, *, _reader: AsyncReader) -> None:
    return None


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("dict", "0.1.0"))
async def read_mapping(elem: AsyncGroup, *, _reader: AsyncReader) -> dict[str, Any]:
    return await _read_members(elem, _reader)


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("anndata", "0.1.0"))
@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("raw", "0.1.0"))
async def read_anndata(elem: AsyncGroup, *, _reader: AsyncReader) -> AnnData:
    keys: tuple[AnnDataElem, ...] = get_literal_members(AnnDataElem)
    return AnnData(**await _read_members(elem, _reader, keys))


async def _read_sparse(
    elem: AsyncGroup, *, _reader: AsyncReader, format: Literal["csr", "csc"]
) -> CSMatrix | CSArray | CupyCSMatrix:
    async def read_dense(k: str) -> DenseType:
        if not isinstance(arr := await elem.getitem(k), AsyncArray):
            msg = f"Expected {k!r} of {elem} to be an array, got {type(arr)}"
            raise ValueError(msg)
        return cast("DenseType", await arr.getitem(...))

    data, indices, indptr = await asyncio.gather(
        *map(read_dense, ["data", "indices", "indptr"])
    )
    return _construct(
        format,
        (data, indices, indptr),
        shape=_shape_from_attrs(elem.attrs),
        dtype=data.dtype,
    )


_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("csr_matrix", "0.1.0"))(
    partial(_read_sparse, format="csr")
)
_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("csc_matrix", "0.1.0"))(
    partial(_read_sparse, format="csc")
)


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("dataframe", "0.2.0"))
async def read_dataframe(elem: AsyncGroup, *, _reader: AsyncReader) -> pd.DataFrame:
    columns, idx_key = _dataframe_keys(elem.attrs)
    *data, index = await _read_children(elem, [*columns, idx_key], _reader)
    return _dataframe_from_parts(
        dict(zip(columns, data, strict=True)), index, columns, idx_key
    )


@_ASYNC_REGISTRY.register_read(AsyncGroup, IOSpec("categorical", "0.2.0"))
async def read_categorical(elem: AsyncGroup, *, _reader: AsyncReader) -> pd.Categorical:
    codes, categories = await _read_children(elem, ["codes", "categories"], _reader)
    return pd.Categorical.from_codes(
        codes=codes,
        categories=categories,
        ordered=bool(_read_attr(elem.attrs, "ordered")),
    )


async def _read_nullable(
    elem: AsyncGroup, *, _reader: AsyncReader, array_type: _BaseMaskedArray
) -> pd.api.extensions.ExtensionArray:
    values, mask = await _read_children(elem, ["values", "mask"], _reader)
    return array_type(values, mask=mask)


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
    values, mask = await _read_children(elem, ["values", "mask"], _reader)
    return _nullable_string_from_parts(values, mask, elem.attrs)
