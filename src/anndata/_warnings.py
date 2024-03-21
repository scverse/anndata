from __future__ import annotations


class WriteWarning(UserWarning):
    pass


class OldFormatWarning(PendingDeprecationWarning):
    """Raised when a file in an old file format is read."""

    pass


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

    pass


class ExperimentalFeatureWarning(Warning):
    """Raised when an unstable experimental feature is used."""

    pass
