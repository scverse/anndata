from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import legacy_api_wrap

if TYPE_CHECKING:
    from typing import Any


_FILE_PREFIXES: tuple[str, ...] = (
    str(Path(__file__).parent),
    str(Path(legacy_api_wrap.__file__).parent),
)

__all__ = [
    "ExperimentalFeatureWarning",
    "ImplicitModificationWarning",
    "OldFormatWarning",
    "WriteWarning",
]


class WriteWarning(UserWarning):
    pass


class OldFormatWarning(PendingDeprecationWarning):
    """Raised when a file in an old file format is read."""


class ImplicitModificationWarning(UserWarning):
    """\
    Raised whenever initializing an object or assigning a property changes
    the type of a part of a parameter or the value being assigned.

    Examples
    ========
    >>> import pandas as pd
    >>> adata = AnnData(obs=pd.DataFrame(index=[0, 1, 2]))  # doctest: +SKIP
    ImplicitModificationWarning: Transforming to str index.
    """


class ExperimentalFeatureWarning(Warning):
    """Raised when an unstable experimental feature is used."""


# needs to be defined outside of `utils` because it’s used in `compat`, `_settings`, etc.
def warn(
    msg: str,
    category: type[Warning],
    source: Any = None,
    *,
    skip_file_prefixes: tuple[str, ...] = (),
    more_file_prefixes: tuple[str, ...] = (),
) -> None:
    if not skip_file_prefixes:
        skip_file_prefixes = (*_FILE_PREFIXES, *more_file_prefixes)
    elif more_file_prefixes:
        msg = "Cannot specify both `skip_file_prefixes` and `more_file_prefixes`"
        raise TypeError(msg)
    warnings.warn(msg, category, source=source, skip_file_prefixes=skip_file_prefixes)  # noqa: TID251
