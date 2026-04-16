from __future__ import annotations

from warnings import warn


def __getattr__(key: str):
    match key:
        case "ExtensionNamespace":
            from scverse_misc import ExtensionNamespace

            msg = (
                "Importing ExtensionNamespace from `types` is deprecated. "
                "Please use scverse_misc instead."
            )
            warn(msg, FutureWarning, stacklevel=2)
            return ExtensionNamespace
        case _:
            msg = f"types has no attribute {key!r}"
            raise AttributeError(msg)
