from __future__ import annotations

from typing import TYPE_CHECKING, cast

from . import AdRef, DFAcc, GraphAcc, ItemAcc, LayerAcc, MetaAcc, MultiAcc
from ._parse_str import _check_vec, _uns_item

if TYPE_CHECKING:
    from collections.abc import Sequence

    from . import AccJson, AdAcc, Idx2D, RefAcc


__all__ = ["parse_json", "to_json"]


def _as_idx_2d(data: Sequence[object]) -> Idx2D | None:
    match data:
        case None, None:
            return (slice(None), slice(None))
        case None, str(i):
            return (slice(None), i)
        case str(i), None:
            return (i, slice(None))
        case _:  # e.g two strings or a longer sequence
            return None


def _item[R: AdRef](a: AdAcc[R], path: Sequence[object]) -> ItemAcc[R]:
    # `__getitem__` validates the keys
    return _uns_item(a.uns, cast("Sequence[str]", path))


def _cast[R: AdRef](
    target: MultiAcc[R] | ItemAcc[R], cast: AccJson
) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]:
    """Parse the trailing cast name and its payload, e.g. `"DF", "c"`."""
    match cast:
        case ["DF"]:
            return target.DF
        case ["DF", str() | None as col]:
            return target.DF[col]
        case ["AD", [*inner]] if isinstance(target, ItemAcc):
            return _parse_json(target.AD, inner)
    msg = f"Cannot cast {target!r} with {list(cast)!r}"
    raise ValueError(msg)


def parse_json[R: AdRef](
    a: AdAcc[R], data: AccJson, *, vec: bool | None = None
) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]:
    parsed = _parse_json(a, data)
    _check_vec(data, vec=vec, actual=isinstance(parsed, AdRef))
    return parsed


def _parse_json[R: AdRef](  # noqa: PLR0911
    a: AdAcc[R], data: AccJson
) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]:
    match data:
        case ["layers", str() | None as l]:
            return a.layers[l]
        case ["layers", str() | None as l, *idx] if idx_2d := _as_idx_2d(idx):
            return a.layers[l][idx_2d]
        case ["obs" | "var" as dim, str() | None as col]:
            acc = a.obs if dim == "obs" else a.var
            return acc.index if col is None else acc[col]
        case ["obsm" | "varm" as dim, str(col)]:
            return (a.obsm if dim == "obsm" else a.varm)[col]
        case ["obsm" | "varm" as dim, str(col), int(i)]:
            return (a.obsm if dim == "obsm" else a.varm)[col][i]
        case ["obsm" | "varm" as dim, str(k), str(), *_]:  # cast name, e.g. "DF"
            return _cast((a.obsm if dim == "obsm" else a.varm)[k], data[2:])
        case ["obsp" | "varp" as dim, str(k)]:
            return (a.obsp if dim == "obsp" else a.varp)[k]
        case ["obsp" | "varp" as dim, str(k), *idx] if idx_2d := _as_idx_2d(idx):
            return (a.obsp if dim == "obsp" else a.varp)[k][idx_2d]
        case ["uns", [*path]]:
            return _item(a, path)
        case ["uns", [*path], str(), *_]:  # cast name, e.g. "AD"
            return _cast(_item(a, path), data[2:])
        case _:
            msg = f"Cannot parse {data!r}"
            raise ValueError(msg)


def _idx_2d_ser(idx: Idx2D) -> tuple[str | None, str | None]:
    i, j = idx
    return (i if isinstance(i, str) else None, j if isinstance(j, str) else None)


def to_json(ref: AdRef | RefAcc) -> list[str | int | AccJson | None]:
    local = _local_json(ref)
    scope = (ref.acc if isinstance(ref, AdRef) else ref).scope
    while scope is not None:  # wrap from the innermost `.AD` outwards
        local = [*_local_json(scope), "AD", local]
        scope = scope.scope
    return local


def _local_json(ref: AdRef | RefAcc) -> list[str | int | AccJson | None]:  # noqa: PLR0911
    """Serialize `ref` relative to the `.AD` cast it’s in, ignoring its scope."""
    match ref:
        case AdRef(LayerAcc() | GraphAcc() as acc, idx):
            return [*_local_json(acc), *_idx_2d_ser(idx)]
        case AdRef(MetaAcc(dim), idx):
            return [dim, idx]
        case AdRef(DFAcc(frame_acc=MultiAcc() | ItemAcc() as frame), col):
            return [*_local_json(frame), "DF", col]
        case AdRef(MultiAcc() as acc, i):
            return [*_local_json(acc), i]
        case LayerAcc(k):
            return ["layers", k]
        case MultiAcc(dim, k):
            return [f"{dim}m", k]
        case GraphAcc(dim, k):
            return [f"{dim}p", k]
        case DFAcc(frame_acc=MultiAcc() | ItemAcc() as frame):
            return [*_local_json(frame), "DF"]
        case ItemAcc(path=path):
            return ["uns", list(path)]
        case _:
            msg = f"Unsupported accessor {ref!r}"
            raise TypeError(msg)
