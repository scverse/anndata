from __future__ import annotations

import re
from typing import TYPE_CHECKING, overload

from . import GraphMapAcc, LayerMapAcc, MetaAcc, MultiMapAcc

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from . import AdAcc, AdRef, GraphAcc, Idx2D, LayerAcc


@overload
def parse[P: AdRef](a: AdAcc[P], spec: str, *, strict: Literal[True] = True) -> P: ...
@overload
def parse[P: AdRef](a: AdAcc[P], spec: str, *, strict: Literal[False]) -> P | None: ...
def parse[P: AdRef](a: AdAcc[P], spec: str, *, strict: bool = True) -> P | None:
    """Create accessor from string."""
    if not strict:
        try:
            parse(a, spec)
        except ValueError:
            return None

    if spec.startswith("X["):
        return _parse_path_2d(lambda _: a, spec)
    if "." not in spec:
        msg = f"Cannot parse accessor {spec!r}"
        raise ValueError(msg)
    acc_name, rest = spec.split(".", 1)
    match getattr(a, acc_name, None):
        case LayerMapAcc() | GraphMapAcc() as acc:
            return _parse_path_2d(acc.__getitem__, rest)
        case MetaAcc() as meta:
            return meta[rest]
        case MultiMapAcc() as multi:
            return _parse_path_multi(multi, rest)
        case None:  # pragma: no cover
            msg = (
                f"Unknown accessor {spec!r}. "
                f"We support '{a.ATTRS}.*' and `AdRef` instances."
            )
            raise ValueError(msg)
    msg = f"Unhandled accessor {spec!r}. This is a bug!"  # pragma: no cover
    raise AssertionError(msg)  # pragma: no cover


def _parse_path_2d[P: AdRef](
    get_vec_acc: Callable[[str], LayerAcc[P] | GraphAcc[P]], spec: str
) -> P:
    if not (
        m := re.fullmatch(r"([^\[]+)\[([^,]+),\s?([^\]]+)\]", spec)
    ):  # pragma: no cover
        msg = f"Cannot parse accessor {spec!r}: should be `name[i,:]`/`name[:,j]`"
        raise ValueError(msg)
    name, i, j = m.groups("")  # "" just for typing
    vec_acc = get_vec_acc(name)
    return vec_acc[_parse_idx_2d(i, j, str)]


def _parse_path_multi[P: AdRef](multi: MultiMapAcc[P], spec: str) -> P:
    if not (m := re.fullmatch(r"([^.]+)\.([\d_]+)", spec)):  # pragma: no cover
        msg = f"Cannot parse multi accessor {spec!r}: should be `name.i`"
        raise ValueError(msg)
    key, i = m.groups("")  # "" just for typing
    return multi[key][int(i)]


def _parse_idx_2d(i: str, j: str, cls: type) -> Idx2D:
    match i, j:
        case ":", ":":
            return slice(None), slice(None)
        case _, ":":
            return cls(i), slice(None)
        case ":", _:
            return slice(None), cls(j)
        case _:  # pragma: no cover
            msg = f"Unknown indices {i!r}, {j!r}"
            raise ValueError(msg)
