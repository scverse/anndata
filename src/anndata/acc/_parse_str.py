from __future__ import annotations

import re
from typing import TYPE_CHECKING, overload

from . import GraphMapAcc, LayerMapAcc, MetaAcc, MultiMapAcc

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from . import AdAcc, AdRef, GraphAcc, Idx2D, LayerAcc, MultiAcc


@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[True] = True, vec: Literal[True]
) -> P: ...
@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[False], vec: Literal[True]
) -> P | None: ...
@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[True] = True, vec: Literal[False]
) -> LayerAcc[P] | MultiAcc[P] | GraphAcc[P]: ...
@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[False], vec: Literal[False]
) -> LayerAcc[P] | MultiAcc[P] | GraphAcc[P] | None: ...
@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[True] = True, vec: None = None
) -> P | LayerAcc[P] | MultiAcc[P] | GraphAcc[P]: ...
@overload
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: Literal[False], vec: None = None
) -> P | LayerAcc[P] | MultiAcc[P] | GraphAcc[P] | None: ...
def parse[P: AdRef](
    a: AdAcc[P], spec: str, *, strict: bool = True, vec: bool | None = None
) -> P | LayerAcc[P] | MultiAcc[P] | GraphAcc[P] | None:
    """Create accessor from string."""
    if not strict:
        try:
            parse(a, spec, vec=vec)
        except ValueError:
            return None

    if spec.startswith("X"):
        _check_vec(spec, vec=vec, actual=(do_vec := "[" in spec))
        return _parse_path_2d(lambda _: a.X, spec) if do_vec else a.X
    if "." not in spec:
        msg = (
            f"Cannot parse accessor {spec!r} that isn’t period-separated, "
            "e.g. `obs.x` for accessing the 1d array at column `x` of `obs`."
        )
        raise ValueError(msg)
    acc_name, rest = spec.split(".", 1)
    acc = getattr(a, acc_name, None)
    if acc is None:  # pragma: no cover
        msg = (
            f"Unknown accessor {spec!r}. "
            f"We support '{{{','.join(a.ATTRS)}}}.*' and `AdRef` instances."
        )
        raise ValueError(msg)
    return _parse_dotted(acc, rest, spec, vec=vec)


def _parse_dotted[P: AdRef](
    acc: LayerMapAcc[P] | GraphMapAcc[P] | MetaAcc[P] | MultiMapAcc[P],
    rest: str,
    spec: str,
    *,
    vec: bool | None,
) -> P | LayerAcc[P] | MultiAcc[P] | GraphAcc[P]:
    match acc:
        case LayerMapAcc() | GraphMapAcc():
            _check_vec(spec, vec=vec, actual=(do_vec := "[" in rest))
            return _parse_path_2d(acc.__getitem__, rest) if do_vec else acc[rest]
        case MetaAcc() as meta:
            _check_vec(spec, vec=vec, actual=True)
            return meta[rest]
        case MultiMapAcc() as multi:
            _check_vec(spec, vec=vec, actual=(do_vec := "." in rest))
            return _parse_path_multi(multi, rest) if do_vec else multi[rest]
    msg = f"Unhandled accessor {spec!r}. This is a bug!"  # pragma: no cover
    raise AssertionError(msg)  # pragma: no cover


def _check_vec(spec: str, *, vec: bool | None, actual: bool) -> None:
    if vec is not None and vec != actual:
        kind = "a vector/`AdRef`" if actual else "a whole container"
        msg = f"Accessor {spec!r} refers to {kind}, but {vec=} was requested"
        raise ValueError(msg)


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
