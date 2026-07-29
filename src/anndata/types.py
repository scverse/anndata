from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from typing import Any, Literal, Self

    from array_api.latest import ArrayNamespace


@runtime_checkable
class SupportsArrayApiBase(Protocol):
    """The part of the array API standard’s array object that we rely on.

    The attributes are read-only, since implementations such as
    :class:`numpy.ndarray` expose them as properties.

    https://data-apis.org/array-api/latest/API_specification/array_object.html
    """

    @property
    def device(self) -> str: ...
    @property
    def dtype(self) -> Any: ...
    @property
    def ndim(self) -> int: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def size(self) -> int: ...

    def __array_namespace__(
        self,
        *,
        api_version: Literal["2021.12", "2022.12", "2023.12", "2024.12"] | None = None,
    ) -> ArrayNamespace: ...
    # `device` and the index are implementation-defined objects
    def to_device(self, device: Any, /, *, stream: int | Any | None = ...) -> Any: ...
    def __getitem__(self, k: Any, /) -> Self: ...


@runtime_checkable
class SupportsArrayApi(SupportsArrayApiBase, Protocol):
    # https://data-apis.org/array-api/latest/design_topics/data_interchange.html
    def __dlpack__(
        self,
        *,
        stream: int | Any | None = None,
        max_version: tuple[int, int] | None = None,
        # the DLPack device type is an `IntEnum`, i.e. usable as a plain `int`
        dl_device: tuple[int, int] | None = None,
        copy: bool | None = None,
    ) -> Any: ...
    def __dlpack_device__(self) -> tuple[int, int]: ...


def __getattr__(key: str):
    match key:
        case "ExtensionNamespace":
            from scverse_misc import ExtensionNamespace

            from .utils import warn

            msg = (
                "Importing ExtensionNamespace from `types` is deprecated. "
                "Please use scverse_misc instead."
            )
            warn(msg, FutureWarning)
            return ExtensionNamespace
        case "SupportsArrayApi":
            return SupportsArrayApi
        case "SupportsArrayApiBase":
            return SupportsArrayApiBase
        case _:
            msg = f"types has no attribute {key!r}"
            raise AttributeError(msg)
