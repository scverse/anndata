from __future__ import annotations

from enum import Enum, Flag, auto
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .types import WriteCompatStr


class WriteFeatures(Flag):
    """Fine-grained write behavior: one flag per on-disk feature.

    Pass a combination of these to a write function’s `compat` argument to pick features individually.
    Flags you leave out are disabled, so to add one to a profile, combine it with the profile’s
    :attr:`~WriteCompat.features`, e.g. `WriteCompat.DEFAULT.features | WriteFeatures.ACCESSOR`.
    """

    ACCESSOR = auto()
    """Write :mod:`anndata.acc` accessors using the `accessor` encoding (anndata ≥0.14)."""

    INDEXED_KEYS = auto()
    """Store mapping keys and dataframe column names in an array instead of as child names (anndata ≥0.14).

    This allows any string as key, e.g. `a/b`, `a:b`, or `NUL` in :attr:`~anndata.AnnData.uns`,
    and makes dataframes follow the `zarr dataframe convention`_.

    .. _zarr dataframe convention: https://github.com/scverse/zarr-convention-dataframe
    """

    if TYPE_CHECKING:  # https://github.com/python/mypy/issues/16712

        def __init__(
            self, value: WriteFeatures | WriteCompat | WriteCompatStr | int, /
        ) -> None: ...

    @classmethod
    def _missing_(cls, value: object) -> WriteFeatures:
        """Make `WriteFeatures(compat)` accept profiles too, normalizing a write function’s `compat`."""
        if isinstance(value, WriteCompat | str):
            return WriteCompat(value).features
        return super()._missing_(value)


class WriteCompat(Enum):
    """Which anndata version’s write behavior to target.

    Each member is the oldest anndata version whose I/O behavior we promise to be
    compatible with, i.e. encodings it cannot read are disallowed.
    See :doc:`/write-compat` for the whole table.
    """

    V0_13 = "0.13"
    """anndata 0.13 (2026-07-07)."""

    V0_14 = "0.14"
    """anndata 0.14 (unreleased).

    Enables

    - the `accessor` encoding for :mod:`anndata.acc` accessors.
    - storing mapping keys and dataframe column names in an array, so any string works as one,
      e.g. `a/b`, `a:b`, or `NUL`.
    """

    DEFAULT = V0_13
    """Alias for the profile write functions use by default, currently :attr:`V0_13`."""

    _value_: str

    @property
    def features(self) -> WriteFeatures:
        """The features this profile enables."""
        return _PROFILE_FEATURES[self]

    def __str__(self) -> str:
        return self._value_


_PROFILE_FEATURES = {
    WriteCompat.V0_13: WriteFeatures(0),
    WriteCompat.V0_14: WriteFeatures.ACCESSOR | WriteFeatures.INDEXED_KEYS,
}
