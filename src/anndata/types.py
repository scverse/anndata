from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from enum import Enum
    from types import ModuleType
    from typing import Any, Literal

    from ._core.anndata import AnnData


@runtime_checkable
class ExtensionNamespace(Protocol):
    """Protocol for extension namespaces.

    Enforces that the namespace initializer accepts a class with the proper `__init__` method.
    Protocol's can't enforce that the `__init__` accepts the correct types. See
    `_check_namespace_signature` for that. This is mainly useful for static type
    checking with mypy and IDEs.
    """

    def __init__(self, adata: AnnData) -> None:
        """
        Used to enforce the correct signature for extension namespaces.
        """


@runtime_checkable
class SupportsArrayApi(Protocol):
    device: str
    shape: tuple[int, ...]

    def __array_namespace__(
        self,
        *,
        api_version: Literal["2021.12", "2022.12", "2023.12", "2024.12"] | None = None,
    ) -> ModuleType: ...
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
