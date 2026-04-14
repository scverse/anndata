from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

from scverse_misc import ExtensionNamespace  # noqa: F401

if TYPE_CHECKING:
    from enum import Enum
    from typing import Any, Literal

    from array_api.latest import ArrayNamespace


@runtime_checkable
class SupportsArrayApi(Protocol):
    device: str
    shape: tuple[int, ...]
    size: int

    def __array_namespace__(
        self,
        *,
        api_version: Literal["2021.12", "2022.12", "2023.12", "2024.12"] | None = None,
    ) -> ArrayNamespace: ...
    def to_device(self, device: str, /, *, stream: int | Any | None = ...) -> Any: ...
    def __dlpack__(
        self,
        *,
        stream: int | Any | None = None,
        max_version: tuple[int, int] | None = None,
        dl_device: tuple[Enum, int] | None = None,
        copy: bool | None = None,
    ) -> Any: ...
    def __dlpack_device__(self) -> tuple[int, int]: ...
