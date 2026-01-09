from __future__ import annotations

import re
from typing import TYPE_CHECKING, overload

from . import GraphAcc, LayerAcc, MetaVecAcc, MultiAcc

if TYPE_CHECKING:
    from typing import Literal

    from . import AdAcc, AdPath, Idx2D


@overload
def parse[P: AdPath](a: AdAcc[P], spec: str, *, strict: Literal[True] = True) -> P: ...
@overload
def parse[P: AdPath](a: AdAcc[P], spec: str, *, strict: Literal[False]) -> P | None: ...
def parse[P: AdPath](a: AdAcc[P], spec: str, *, strict: bool = True) -> P | None:
    """Create accessor from string."""
    if not strict:
        try:
            parse(a, spec)
        except ValueError:
            return None

    if "." not in spec:
        msg = f"Cannot parse accessor {spec!r}"
        raise ValueError(msg)
    acc, rest = spec.split(".", 1)
    match getattr(a, acc, None):
        # TODO: X
        case LayerAcc() as layers:
            return _parse_path_layer(layers, rest)
        case MetaVecAcc() as meta:
            return meta[rest]
        case MultiAcc() as multi:
            return _parse_path_multi(multi, rest)
        case GraphAcc():
            msg = "TODO"
            raise NotImplementedError(msg)
        case None:  # pragma: no cover
            msg = (
                f"Unknown accessor {spec!r}. "
                f"We support '{a.ATTRS}.*' and `AdPath` instances."
            )
            raise ValueError(msg)
    msg = f"Unhandled accessor {spec!r}. This is a bug!"  # pragma: no cover
    raise AssertionError(msg)  # pragma: no cover


def _parse_path_layer[P: AdPath](layers: LayerAcc[P], spec: str) -> P:
    if not (
        m := re.fullmatch(r"([^\[]+)\[([^,]+),\s?([^\]]+)\]", spec)
    ):  # pragma: no cover
        msg = f"Cannot parse layer accessor {spec!r}: should be `name[i,:]`/`name[:,j]`"
        raise ValueError(msg)
    layer, i, j = m.groups("")  # "" just for typing
    return layers[layer][_parse_idx_2d(i, j, str)]


def _parse_path_multi[P: AdPath](multi: MultiAcc[P], spec: str) -> P:
    if not (m := re.fullmatch(r"([^.]+)\.([\d_]+)", spec)):  # pragma: no cover
        msg = f"Cannot parse multi accessor {spec!r}: should be `name.i`"
        raise ValueError(msg)
    key, i = m.groups("")  # "" just for typing
    return multi[key][int(i)]


def _parse_idx_2d[Idx: int | str](i: str, j: str, cls: type[Idx]) -> Idx2D[Idx]:
    match i, j:
        case _, ":":
            return cls(i), slice(None)
        case ":", _:
            return slice(None), cls(j)
        case _:  # pragma: no cover
            msg = f"Unknown indices {i!r}, {j!r}"
            raise ValueError(msg)
