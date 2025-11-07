from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from typing import TypeAlias
    from warnings import _ActionKind

    # action, message, category, module, lineno
    WarningFilter: TypeAlias = tuple[_ActionKind, str, type[Warning], str, int]


F = TypeVar("F", bound=Callable)


def doctest_needs(mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    def decorator(func: F) -> F:
        func._doctest_needs = mod
        return func

    return decorator


def doctest_filterwarnings(
    action: _ActionKind,
    message: str = r"",
    category: type[Warning] = Warning,
    module: str = r"",
    lineno: int = 0,
) -> Callable[[F], F]:
    """Mark function with warning filter."""
    filter: WarningFilter = (action, message, category, module, lineno)

    def decorator(func: F) -> F:
        if not hasattr(func, "_doctest_warning_filter"):
            func._doctest_warning_filter = []
        func._doctest_warning_filter.insert(0, filter)
        return func

    return decorator
