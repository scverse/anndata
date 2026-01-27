from __future__ import annotations

from typing import TYPE_CHECKING

from . import AdRef

if TYPE_CHECKING:
    from collections.abc import Sequence

    from . import AdAcc, Idx2D


def as_idx_2d[Idx: int | str](
    data: Sequence[str | int | None], typ: type[Idx]
) -> Idx2D[Idx] | None:
    if (
        len(data) != 2
        or not all(isinstance(i, typ | None) for i in data)
        or sum(i is None for i in data) == 0
    ):
        return None
    return tuple(slice(None) if i is None else i for i in data)  # type: ignore


def parse_json[R: AdRef](a: AdAcc[R], data: Sequence[str | int | None]) -> R:
    match data:
        case ["l", str() | None as layer, *idx] if idx := as_idx_2d(idx, str):
            return a.layers[layer][idx]
        case ["o" | "v" as ax, str(col)]:
            return (a.obs if ax == "o" else a.var)[col]
        case ["om" | "vm" as ax, str(col), int(i)]:
            return (a.obsm if ax == "om" else a.varm)[col][i]
        case ["op" | "vp" as ax, str(k), *idx] if idx := as_idx_2d(idx, str):
            return (a.obsp if ax == "op" else a.varp)[k][idx]
        case _:
            msg = f"Cannot parse {data!r}"
            raise ValueError(msg)
