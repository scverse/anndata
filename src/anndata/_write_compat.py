from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Literal

from packaging.version import Version

from ._settings import settings


class WriteCompat(Version, Enum):
    """Which anndata version’s write behavior to target.

    Each member is the oldest anndata version whose I/O behavior we promise to be
    compatible with: encodings it cannot read are disallowed,
    and restrictions added after it are relaxed.
    See :doc:`/write-compat` for the whole table.
    """

    V0_12 = "0.12"
    """anndata 0.12 (2025-07-16).

    Relaxes the restrictions added in 0.13, i.e. allows writing
    3-dimensional :attr:`~anndata.AnnData.X`/:attr:`~anndata.AnnData.layers` and
    `/` in `h5ad` keys of `obs`, `var` and `uns`.
    """

    V0_13 = "0.13"
    """anndata 0.13 (2026-07-07).

    No relaxations.
    """

    DEFAULT = V0_13
    """Alias for the profile write functions use by default, currently :attr:`V0_13`."""

    if TYPE_CHECKING:
        # the enum lookup also takes members, but type checkers see `Version.__init__(str)`
        def __init__(self, value: WriteCompat | WriteCompatStr, /) -> None: ...


type WriteCompatStr = Literal["0.12", "0.13"]
"""The values of :class:`WriteCompat`, which write functions also accept."""


def forward_slash_disallowed(compat: WriteCompat) -> bool:
    """Whether `/` is disallowed in `h5ad` keys, honoring the deprecated setting."""
    # bypass the descriptor, as our own reads shouldn’t emit its deprecation warning
    if (explicit := settings.__dict__["disallow_forward_slash_in_h5ad"]) is not None:
        return explicit
    return compat >= WriteCompat.V0_13
