"""Find the :class:`~anndata.acc.AdAcc` that (de)serializes an accessor, via entry points."""

from __future__ import annotations

from importlib.metadata import entry_points
from typing import TYPE_CHECKING

from packaging.utils import canonicalize_name

from . import AdRef, GraphAcc, LayerAcc, MultiAcc

if TYPE_CHECKING:
    from importlib.metadata import EntryPoint
    from typing import TypeAlias

    from . import AdAcc


__all__ = ["WritableAcc", "load_acc", "to_json"]

GROUP = "anndata.acc"

# Accessor base classes with a registered writer (subclasses dispatch to these).
# `MetaAcc` (e.g. `A.obs`) is missing on purpose: `to_json` cannot serialize it
WritableAcc: TypeAlias = AdRef | LayerAcc | MultiAcc | GraphAcc  # noqa: UP040  # used at runtime


def _project_of(elem: WritableAcc) -> str:
    """Get the name of the entry point pointing into the top-level package that defines `elem`’s classes."""
    types = {type(elem), *([type(elem.acc)] if isinstance(elem, AdRef) else [])}
    mods = {t.__module__.partition(".")[0] for t in types}
    # extensions build on anndata’s classes, so those don’t determine the project
    match list(mods - {"anndata"} or mods):
        case [mod]:
            pass
        case several:
            msg = f"Cannot write {elem!r}: it mixes classes from {sorted(several)}."
            raise TypeError(msg)
    for ep in entry_points(group=GROUP):
        if ep.module.partition(".")[0] == mod:
            return _check_name(ep).name
    msg = (
        f"Cannot write {elem!r}: no {GROUP!r} entry point registered "
        f"for the package defining it ({mod!r})."
    )
    # also hit if the `AdAcc` lives in another top-level package than the classes
    raise TypeError(msg)


def load_acc(project: str) -> AdAcc:
    """Load the :class:`~anndata.acc.AdAcc` registered by `project`."""
    if (ep := next(iter(entry_points(group=GROUP, name=project)), None)) is None:
        msg = f"No {GROUP!r} entry point named {project!r}. Is {project!r} installed?"
        raise LookupError(msg)
    return _check_name(ep).load()


def _check_name(ep: EntryPoint) -> EntryPoint:
    """Make sure `ep` is named after the project registering it, so names are unique and stable."""
    # `dist` is `None` for entry points not from an installed distribution
    if ep.dist is not None and canonicalize_name(ep.dist.name) != canonicalize_name(
        ep.name
    ):
        msg = (
            f"{GROUP!r} entry point {ep.name!r} must be named after "
            f"the project registering it ({ep.dist.name!r})."
        )
        raise ValueError(msg)
    return ep


def to_json(elem: WritableAcc) -> tuple[str, list[str | int | None]]:
    """Serialize `elem` using the :class:`~anndata.acc.AdAcc` of the project defining it.

    Raises :exc:`TypeError` if that doesn’t round-trip, so nothing gets silently lost.
    """
    project = _project_of(elem)
    adacc = load_acc(project)
    data = adacc.to_json(elem)
    if adacc.from_json(data) != elem:
        msg = f"Cannot write {elem!r}: {adacc!r} serializes it lossily as {data!r}."
        raise TypeError(msg)
    return project, data
