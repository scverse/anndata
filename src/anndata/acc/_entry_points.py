"""Find the :class:`~anndata.acc.AdAcc` that (de)serializes an accessor, via entry points."""

from __future__ import annotations

from functools import cache
from importlib.metadata import entry_points
from typing import TYPE_CHECKING

from packaging.utils import canonicalize_name

from . import AdRef

if TYPE_CHECKING:
    from importlib.metadata import EntryPoint

    from . import AdAcc, RefAcc


__all__ = ["load_acc", "to_json"]

GROUP = "anndata.acc"


@cache  # scanning takes milliseconds. Plugins installed while running need a restart
def _entry_points() -> list[tuple[EntryPoint, str | None]]:
    """Get entry points with the project registering them (`dist` is `None` if not installed)."""
    return [(ep, ep.dist.name if ep.dist else None) for ep in entry_points(group=GROUP)]


def _project_of(elem: AdRef | RefAcc) -> str:
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
    for ep, dist in _entry_points():
        if ep.module.partition(".")[0] == mod:
            return _check_name(ep, dist).name
    msg = (
        f"Cannot write {elem!r}: no {GROUP!r} entry point registered "
        f"for the package defining it ({mod!r})."
    )
    # also hit if the `AdAcc` lives in another top-level package than the classes
    raise TypeError(msg)


def load_acc(project: str) -> AdAcc:
    """Load the :class:`~anndata.acc.AdAcc` registered by `project`."""
    for ep, dist in _entry_points():
        if ep.name == project:
            return _check_name(ep, dist).load()
    msg = f"No {GROUP!r} entry point named {project!r}. Is {project!r} installed?"
    raise LookupError(msg)


def _check_name(ep: EntryPoint, dist: str | None) -> EntryPoint:
    """Make sure `ep` is named after the project registering it, so names are unique and stable."""
    if dist is not None and canonicalize_name(dist) != canonicalize_name(ep.name):
        msg = (
            f"{GROUP!r} entry point {ep.name!r} must be named after "
            f"the project registering it ({dist!r})."
        )
        raise ValueError(msg)
    return ep


def to_json(elem: AdRef | RefAcc) -> tuple[str, list[str | int | None]]:
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
