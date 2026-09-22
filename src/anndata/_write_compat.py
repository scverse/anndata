from __future__ import annotations

import operator
from dataclasses import dataclass, fields
from enum import Enum
from types import MethodType
from typing import TYPE_CHECKING

from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Self

    from .types import WriteCompatStr


def _min_version(version: str) -> Callable[[], bool]:
    """Make a flag’s default factory: enabled in profiles from `version` on.

    Called without arguments, it checks :attr:`WriteCompat.DEFAULT`.
    """

    def enabled(profile: str | None = None) -> bool:
        if profile is None:
            profile = WriteCompat.DEFAULT._value_
        return Version(profile) >= Version(version)

    return enabled


class _CompareFeatures:
    """Comparison method applying the `operator` of the same name to both configs’ features."""

    def __set_name__(self, owner: type, name: str) -> None:
        self.op: Callable[[frozenset[str], frozenset[str]], bool] = getattr(
            operator, name
        )

    def __get__(
        self, obj: WriteConfig | None, owner: type | None = None
    ) -> Callable[..., bool]:
        return self if obj is None else MethodType(self, obj)

    def __call__(self, obj: WriteConfig, other: object) -> bool:
        if not isinstance(other, WriteConfig):
            return NotImplemented
        return self.op(obj._features, other._features)


@dataclass(kw_only=True, eq=False)
class WriteConfig:
    """Fine-grained write behavior: one flag per on-disk feature.

    Pass this to a write function’s `compat` argument to pick features individually.
    Flags you leave out are taken from :attr:`WriteCompat.DEFAULT`.
    :class:`WriteCompat` profiles are preset instances of this.

    Configs compare by their sets of enabled features,
    e.g. `a < b` means `b` enables every feature `a` does, and more.
    """

    @property
    def _features(self) -> frozenset[str]:
        return frozenset(f.name for f in fields(self) if getattr(self, f.name))

    def __hash__(self) -> int:
        return hash(self._features)

    __eq__ = _CompareFeatures()
    __lt__ = _CompareFeatures()
    __le__ = _CompareFeatures()
    __gt__ = _CompareFeatures()
    __ge__ = _CompareFeatures()


class WriteCompat(WriteConfig, Enum):
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
    """

    DEFAULT = V0_13
    """Alias for the profile write functions use by default, currently :attr:`V0_13`."""

    _value_: str

    if TYPE_CHECKING:  # https://github.com/python/mypy/issues/16712

        def __init__(self, version: WriteCompat | WriteCompatStr, /) -> None: ...

    else:

        def __new__(cls, version: str) -> Self:
            self = object.__new__(cls)
            self._value_ = version
            return self

        def __init__(self, version: str) -> None:
            # explicit `version`, since `DEFAULT` doesn’t exist yet
            for f in fields(WriteConfig):
                setattr(self, f.name, f.default_factory(version))

    @classmethod
    def or_config(cls, obj: WriteCompatStr | WriteConfig) -> WriteConfig:
        return cls(obj) if isinstance(obj, str) else obj

    def __str__(self) -> str:
        return self._value_

    def __repr__(self) -> str:  # Enum’s would use the dataclass repr for `_value_`
        return f"<{type(self).__name__}.{self._name_}: {self._value_!r}>"
