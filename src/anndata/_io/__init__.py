from __future__ import annotations

import warnings

__all__: list[str] = []


def __getattr__(key: str):
    from .. import io

    attr = getattr(io, key)
    msg = (
        f"Importing {key} from `anndata._io` is deprecated. "
        "Please use anndata.io instead."
    )
    warnings.warn(msg, FutureWarning, stacklevel=2)
    return attr
