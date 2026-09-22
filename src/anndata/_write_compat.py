from __future__ import annotations

from enum import Enum
from functools import total_ordering
from typing import TYPE_CHECKING

from packaging.version import Version

if TYPE_CHECKING:
    from .types import WriteCompatStr


@total_ordering
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

    - the `sequence` encoding, which stores a sequence’s elements individually instead
      of converting it to an array. That allows writing heterogeneous, ragged and
      nested sequences, i.e. arbitrary JSON-like structures in
      :attr:`~anndata.AnnData.uns`.
    - the `accessor` encoding for :mod:`anndata.acc` accessors.
    - escaping keys that a file system may choke on, e.g. `a/b`, `a:b`, or `NUL`.
    """

    DEFAULT = V0_13
    """Alias for the profile write functions use by default, currently :attr:`V0_13`."""

    _version: Version

    if TYPE_CHECKING:  # https://github.com/python/mypy/issues/16712

        def __init__(self, version: WriteCompat | WriteCompatStr, /) -> None: ...

    else:

        def __init__(self, version: str, /) -> None:
            self._version = Version(version)

    def __str__(self) -> str:
        return str(self._version)

    @staticmethod
    def _coerce(other: object) -> Version | None:
        if isinstance(other, Version):
            return other
        if isinstance(other, WriteCompat):
            return other._version
        return None

    def __eq__(self, other: object) -> bool:
        v = self._coerce(other)
        return NotImplemented if v is None else self._version == v

    def __lt__(self, other: object) -> bool:
        v = self._coerce(other)
        return NotImplemented if v is None else self._version < v

    def __hash__(self) -> int:
        return hash(self._version)
