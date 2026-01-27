from __future__ import annotations

from typing import TYPE_CHECKING

from . import AdRef, GraphAcc, LayerAcc, MetaAcc, MultiAcc

if TYPE_CHECKING:
    from collections.abc import Sequence

    from . import AdAcc, Idx2D


__all__ = ["parse_json", "to_json"]


def _as_idx_2d(data: Sequence[str | int | None]) -> Idx2D | None:
    match data:
        case None, None:
            return (slice(None), slice(None))
        case None, str(i):
            return (slice(None), i)
        case str(i), None:
            return (i, slice(None))
        case _:  # e.g two strings or a longer sequence
            return None


def parse_json[R: AdRef](a: AdAcc[R], data: Sequence[str | int | None]) -> R:
    match data:
        case ["layers", str() | None as l, *idx] if idx := _as_idx_2d(idx):
            layer = a if l is None else a.layers[l]
            return layer[idx]
        case ["obs" | "var" as ax, str() | None as col]:
            acc = a.obs if ax == "obs" else a.var
            return acc.index if col is None else acc[col]
        case ["obsm" | "varm" as ax, str(col), int(i)]:
            return (a.obsm if ax == "obsm" else a.varm)[col][i]
        case ["obsp" | "varp" as ax, str(k), *idx] if idx := _as_idx_2d(idx):
            return (a.obsp if ax == "obsp" else a.varp)[k][idx]
        case _:
            msg = f"Cannot parse {data!r}"
            raise ValueError(msg)


def _idx_2d_ser(idx: Idx2D) -> tuple[str | None, None] | tuple[None, str]:
    return tuple(i if isinstance(i, str) else None for i in idx)  # type: ignore


def to_json(ref: AdRef) -> list[str | int | None]:
    match ref.acc:
        case LayerAcc(k):
            return ["layers", k, *_idx_2d_ser(ref.idx)]
        case MetaAcc(ax):
            return [ax, ref.idx]
        case MultiAcc(ax, k):
            return [f"{ax}m", k, ref.idx]
        case GraphAcc(ax, k):
            return [f"{ax}p", k, *_idx_2d_ser(ref.idx)]
        case _:
            msg = f"Unsupported vector accessor {ref.acc!r}"
            raise AssertionError(msg)
