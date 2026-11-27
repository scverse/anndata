from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol, overload, runtime_checkable

if TYPE_CHECKING:
    from collections.abc import MutableMapping
    from types import CapsuleType
    from typing import Any, Self

    import pandas as pd
    from array_api import latest as array_api

    from ._core.xarray import Dataset2D
    from .acc import (
        AdAcc,
        AdRef,
        Array,
        DataFrameLike,
        DFAcc,
        GraphAcc,
        ItemAcc,
        LayerAcc,
        MultiAcc,
        RefAcc,
    )
    from .typing import AlignedAny, AlignedArray, Index


__all__ = [
    "AnnDataLike",
    "AnnDataLikeMetadataOnly",
    "SupportsArrayApi",
    "SupportsArrayApiBase",
    "WriteCompatStr",
]


type WriteCompatStr = Literal["0.13", "0.14"]
"""The values of :class:`~anndata.WriteCompat`, which write functions also accept."""


@runtime_checkable
class SupportsArrayApiBase(Protocol):
    """The part of the array API standard’s array object that we rely on.

    The attributes are read-only, since implementations such as
    :class:`numpy.ndarray` expose them as properties.

    https://data-apis.org/array-api/latest/API_specification/array_object.html
    """

    @property
    def device(self) -> object: ...
    @property
    def dtype(self) -> object: ...
    @property
    def ndim(self) -> int: ...
    # lazy arrays may have unknown dimensions
    @property
    def shape(self) -> tuple[int | None, ...]: ...
    @property
    def size(self) -> int | None: ...
    @property
    def T(self) -> Self: ...

    def __array_namespace__(
        self,
        *,
        api_version: Literal["2021.12", "2022.12", "2023.12", "2024.12"] | None = None,
    ) -> array_api.ArrayNamespace[array_api.Array[object, object], object, object]:
        """See :meth:`array_api.array.__array_namespace__`."""

    # `device` and the index are implementation-defined objects
    def to_device(
        self,
        device: Any,  # noqa: ANN401  # protocol parameter
        /,
        *,
        stream: int | Any | None = ...,  # noqa: ANN401  # protocol parameter
    ) -> object: ...
    def __getitem__(self, k: Any, /) -> Self:  # noqa: ANN401  # protocol parameter
        """See :meth:`array_api.array.__getitem__`."""


@runtime_checkable
class SupportsArrayApi(SupportsArrayApiBase, Protocol):
    # https://data-apis.org/array-api/latest/design_topics/data_interchange.html
    def __dlpack__(
        self,
        *,
        stream: int | Any | None = None,  # noqa: ANN401  # protocol parameter
        max_version: tuple[int, int] | None = None,
        # the DLPack device type is an `IntEnum`, i.e. usable as a plain `int`
        dl_device: tuple[int, int] | None = None,
        copy: bool | None = None,
    ) -> CapsuleType:
        """See :meth:`array_api.array.__dlpack__`."""

    def __dlpack_device__(self) -> tuple[int, int]:
        """See :meth:`array_api.array.__dlpack_device__`."""


@runtime_checkable
class AnnDataLikeMetadataOnly(Protocol):
    """The AnnData API excluding data matrices.

    This can be used for `isinstance()` checks in code that only accesses metadata.
    Indexing with a bare `str` is not portable:
    :class:`~anndata.AnnData` selects an observation, :class:`~mudata.MuData` a modality.
    """

    # read-only, so implementations may have setters or not
    @property
    def shape(self) -> tuple[int, int, *tuple[int, ...]]:
        """`(n_obs, n_vars, ...)`, implementations may have more dimensions."""

    @property
    def n_obs(self) -> int: ...
    @property
    def n_vars(self) -> int: ...
    @property
    def obs(self) -> pd.DataFrame | Dataset2D: ...
    @property
    def var(self) -> pd.DataFrame | Dataset2D: ...
    @property
    def obs_names(self) -> pd.Index: ...
    @property
    def var_names(self) -> pd.Index: ...
    @property
    def obsm(self) -> MutableMapping[str, AlignedAny]: ...
    @property
    def varm(self) -> MutableMapping[str, AlignedAny]: ...
    @property
    def obsp(self) -> MutableMapping[str, AlignedArray]: ...
    @property
    def varp(self) -> MutableMapping[str, AlignedArray]: ...
    @property
    def uns(self) -> MutableMapping[str, Any]: ...

    @overload
    def __getitem__(self, index: DFAcc) -> DataFrameLike: ...
    @overload
    def __getitem__(self, index: ItemAcc) -> object: ...
    @overload
    def __getitem__(self, index: MultiAcc | GraphAcc) -> AlignedArray: ...
    @overload
    def __getitem__(self, index: AdRef) -> Array: ...
    # separate from `Index` so `MuData`’s `str` overload (a modality) matches: mypy matches overloads in order
    @overload
    def __getitem__(self, index: str) -> AnnDataLikeMetadataOnly: ...
    # separate from Index so `MuData` can return `MuData` for `Index` and `AnnData` for `AdAcc`
    @overload
    def __getitem__(self, index: AdAcc) -> AnnDataLikeMetadataOnly: ...
    @overload
    def __getitem__(self, index: Index) -> AnnDataLikeMetadataOnly: ...
    def __setitem__(self, index: AdRef | RefAcc | AdAcc, value: object) -> None: ...
    def __len__(self) -> int: ...


@runtime_checkable
class AnnDataLike(AnnDataLikeMetadataOnly, Protocol):
    """The full AnnData API."""

    X: AlignedArray | None

    @property
    def layers(self) -> MutableMapping[str | None, AlignedArray]: ...

    # overloads replace the parent’s as a whole, so repeat them and add `LayerAcc`
    @overload
    def __getitem__(self, index: DFAcc) -> DataFrameLike: ...
    @overload
    def __getitem__(self, index: ItemAcc) -> object: ...
    @overload
    def __getitem__(self, index: LayerAcc | MultiAcc | GraphAcc) -> AlignedArray: ...
    @overload
    def __getitem__(self, index: AdRef) -> Array: ...
    # separate from `Index` so `MuData`’s `str` overload (a modality) matches: mypy matches overloads in order
    @overload
    def __getitem__(self, index: str) -> AnnDataLike: ...
    # separate from Index so `MuData` can return `MuData` for `Index` and `AnnData` for `AdAcc`
    @overload
    def __getitem__(self, index: AdAcc) -> AnnDataLike: ...
    @overload
    def __getitem__(self, index: Index) -> AnnDataLike: ...


def __getattr__(key: str) -> object:
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
        case _:
            msg = f"types has no attribute {key!r}"
            raise AttributeError(msg)
