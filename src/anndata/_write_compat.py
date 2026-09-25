from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Literal

from packaging.version import Version


class WriteCompat(Version, Enum):
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

    if TYPE_CHECKING:
        # the enum lookup also takes members, but type checkers see `Version.__init__(str)`
        def __init__(self, value: WriteCompat | WriteCompatStr, /) -> None: ...


type WriteCompatStr = Literal["0.13", "0.14"]
"""The values of :class:`WriteCompat`, which write functions also accept."""
