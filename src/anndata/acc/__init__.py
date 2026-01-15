"""Accessor classes for AnnData interface."""

from __future__ import annotations

import abc
from collections.abc import Hashable
from dataclasses import KW_ONLY, dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING, ClassVar, cast, overload

import pandas as pd
import scipy.sparse as sp
from numpy.typing import NDArray

from anndata import AnnData

if TYPE_CHECKING:
    from collections.abc import Callable, Collection
    from typing import Any, Literal, Self, TypeIs

    from numpy.typing import NDArray

    from anndata import AnnData

    from . import hv
else:
    # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    type P = AdPath
    type I = Hashable

type Vector = pd.api.extensions.ExtensionArray | NDArray[Any]

# full slices: e.g. a[:, 5], a[18, :], or a[:, :]
type Sf = slice[None, None, None]
type Idx2D[Idx: int | str] = tuple[Idx | Sf, Sf] | tuple[Sf, Idx | Sf]
type Idx2DList[Idx: int | str] = tuple[list[Idx], Sf] | tuple[Sf, list[Idx]]
type AdPathFunc[I] = Callable[[AnnData, I], Vector]
type Axes = Collection[Literal["obs", "var"]]


__all__ = [
    "A",
    "AdAcc",
    "AdPath",
    "GraphAcc",
    "GraphVecAcc",
    "LayerAcc",
    "LayerVecAcc",
    "MetaVecAcc",
    "MultiAcc",
    "MultiVecAcc",
    "VecAcc",
    "hv",
]


class AdPath[I: Hashable]:
    """A path referencing an array in an AnnData object."""

    acc: VecAcc[Self, I]
    """The vector accessor containing information about this array."""

    idx: I
    """The index specifying the array."""

    __match_args__: ClassVar = ("acc", "idx")

    def __init__(self, acc: VecAcc[Self, I], idx: I) -> None:
        self.acc = acc
        self.idx = idx

    @cached_property
    def axes(self) -> Axes:
        """Which axes are spanned by the array?"""
        return self.acc.axes(self.idx)

    def __hash__(self) -> int:
        return hash((self.acc, self.idx))

    def __eq__(self, value: object) -> bool:
        match value:
            case AdPath():
                return self.acc == value.acc and self.idx == value.idx
            case str():
                # TODO: more elegant: access `A` through `acc`?
                a = (
                    A
                    if type(self).__eq__ is AdPath.__eq__
                    else AdAcc(path_class=type(self))
                )
                return self == a.resolve(value, strict=False)
        return False

    @cached_property
    def __repr(self) -> str:
        idx_repr = self.acc.idx_repr(self.idx).replace("slice(None, None, None)", ":")
        return f"{self.acc!r}{idx_repr}"

    def __str__(self) -> str:
        return self.__repr

    def __repr__(self) -> str:
        return self.__repr

    def __call__(self, adata: AnnData) -> Vector:
        """Retrieve referenced array from AnnData."""
        return self.acc.get(adata, self.idx)


@dataclass(frozen=True)
class VecAcc[P: AdPath[I], I](abc.ABC):  # type: ignore
    """Abstract base class for vector accessors."""

    _: KW_ONLY
    path_class: type[P]

    def process_idx(self, idx: Any, /) -> I:
        self.axes(idx)
        return idx

    def __getitem__(self, idx: Any, /) -> P:
        idx = self.process_idx(idx)
        return self.path_class(self, idx)  # type: ignore

    @abc.abstractmethod
    def axes(self, idx: I, /) -> Axes: ...
    @abc.abstractmethod
    def __repr__(self, /) -> str: ...
    @abc.abstractmethod
    def idx_repr(self, idx: I, /) -> str: ...
    @abc.abstractmethod
    def get(self, adata: AnnData, idx: I, /) -> Vector: ...


@dataclass(frozen=True)
class LayerAcc[P: AdPath]:
    """Accessor for layers."""

    _: KW_ONLY
    path_class: type[P]

    def __getitem__(self, k: str, /) -> LayerVecAcc[P]:
        if not isinstance(k, str):
            msg = f"Unsupported layer {k!r}"
            raise TypeError(msg)
        return LayerVecAcc(k, path_class=self.path_class)

    def __repr__(self) -> str:
        return "A.layers"


@dataclass(frozen=True)
class LayerVecAcc[P: AdPath[Idx2D[str]]](VecAcc[P, "Idx2D[str]"]):
    """Accessor for layer vectors."""

    k: str | None

    @overload
    def __getitem__(self, idx: Idx2D[str], /) -> P: ...
    @overload
    def __getitem__(self, idx: Idx2DList[str], /) -> list[P]: ...
    def __getitem__(self, idx: Idx2D[str] | Idx2DList[str], /) -> P | list[P]:
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def axes(self, idx: Idx2D[str], /) -> Axes:
        return _idx2axes(idx)

    def __repr__(self) -> str:
        return f"A.layers[{self.k!r}]"

    def idx_repr(self, idx: Idx2D[str]) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def get(self, adata: AnnData, idx: Idx2D[str], /) -> Vector:
        ver_or_mat = adata[idx].X if self.k is None else adata[idx].layers[self.k]
        if isinstance(ver_or_mat, sp.spmatrix | sp.sparray):
            ver_or_mat = ver_or_mat.toarray()
        # TODO: pandas
        axes = self.axes(idx)
        return ver_or_mat.flatten() if len(axes) == 1 else ver_or_mat


@dataclass(frozen=True)
class MetaVecAcc[P: AdPath[str | type[pd.Index]]](VecAcc[P, str | type[pd.Index]]):
    """Accessor for metadata (obs/var)."""

    ax: Literal["obs", "var"]

    @property
    def index(self) -> P:
        """Index accessor."""
        return self[pd.Index]

    def process_idx(self, k: object, /) -> str | type[pd.Index]:
        if k is not pd.Index and not isinstance(k, str):
            msg = f"Unsupported {self.ax} column {k!r}"
            raise TypeError(msg)
        return k

    @overload
    def __getitem__[K: str | type[pd.Index]](self, k: K, /) -> P: ...
    @overload
    def __getitem__[K: str | type[pd.Index]](self, k: list[K], /) -> list[P]: ...
    def __getitem__[K: str | type[pd.Index]](self, k: K | list[K], /) -> P | list[P]:
        if isinstance(k, list):
            return [self[i] for i in k]

        return super().__getitem__(k)

    def axes(self, k: str, /) -> Axes:
        return {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}"

    def idx_repr(self, k: str | type[pd.Index]) -> str:
        return ".index" if k is pd.Index else f"[{k!r}]"

    def get(self, adata: AnnData, k: str | type[pd.Index], /) -> Vector:
        df = cast("pd.DataFrame", getattr(adata, self.ax))
        return df.index.values if k is pd.Index else df[k].values


@dataclass(frozen=True)
class MultiAcc[P: AdPath]:
    """Accessor for multi-dimensional array containers (obsm/varm)."""

    ax: Literal["obs", "var"]
    _: KW_ONLY
    path_class: type[P]

    def __getitem__(self, k: str, /) -> MultiVecAcc[P]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax}m key {k!r}"
            raise TypeError(msg)
        return MultiVecAcc(self.ax, k, path_class=self.path_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}m"


@dataclass(frozen=True)
class MultiVecAcc[P: AdPath[int]](VecAcc[P, int]):
    """Accessor for arrays from multi-dimensional containers (obsm/varm)."""

    ax: Literal["obs", "var"]
    k: str

    @staticmethod
    def process_idx[T](i: T | tuple[Sf, T], /) -> T:
        if isinstance(i, tuple):
            if len(i) != 2 or i[0] != slice(None):
                msg = f"Unsupported slice {i!r}"
                raise ValueError(msg)
            i = cast("T", i[1])
        if not isinstance(i, list | int):
            msg = f"Unsupported index {i!r}"
            raise TypeError(msg)
        return i

    @overload
    def __getitem__(self, i: int | tuple[Sf, int], /) -> P: ...
    @overload
    def __getitem__(self, i: list[int] | tuple[Sf, list[int]], /) -> list[P]: ...
    def __getitem__(
        self, i: int | list[int] | tuple[Sf, int | list[int]], /
    ) -> P | list[P]:
        i = self.process_idx(i)
        if isinstance(i, list):
            return [self[j] for j in i]
        return super().__getitem__(i)

    def axes(self, i: int, /) -> Axes:
        return {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}m[{self.k!r}]"

    def idx_repr(self, i: int) -> str:
        return f"[:, {i!r}]"

    def get(self, adata: AnnData, i: int, /) -> Vector:
        return getattr(adata, f"{self.ax}m")[self.k][:, i]


@dataclass(frozen=True)
class GraphAcc[P: AdPath]:
    """Accessor for graph containers (obsp/varp)."""

    ax: Literal["obs", "var"]
    _: KW_ONLY
    path_class: type[P]

    def __getitem__(self, k: str, /) -> GraphVecAcc[P]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax}p key {k!r}"
            raise TypeError(msg)
        return GraphVecAcc(self.ax, k, path_class=self.path_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}p"


@dataclass(frozen=True)
class GraphVecAcc[P: AdPath[Idx2D[str]]](VecAcc[P, "Idx2D[str]"]):
    """Accessor for arrays from graph containers (obsp/varp)."""

    ax: Literal["obs", "var"]
    k: str

    def process_idx(self, idx: Idx2D[str], /) -> Idx2D[str]:
        if not all(isinstance(i, str | slice) for i in idx):
            msg = f"Unsupported index {idx!r}"
            raise TypeError(msg)
        if sum(isinstance(i, slice) for i in idx) not in {1, 2}:
            msg = (
                f"Unsupported index {idx!r}: "
                f"should be ['c1', :], ['c1', 'c2'], or [:, 'c2']"
            )
            raise TypeError(msg)
        return idx

    @overload
    def __getitem__(self, idx: Idx2D[str], /) -> P: ...
    @overload
    def __getitem__(self, idx: Idx2DList[str], /) -> list[P]: ...
    def __getitem__(self, idx: Idx2D[str] | Idx2DList[str], /) -> P | list[P]:
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def axes(self, idx: Idx2D[str], /) -> Axes:
        n_slices = sum(isinstance(i, slice) for i in idx)
        return (self.ax,) * n_slices if n_slices > 1 else {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}p[{self.k!r}]"

    def idx_repr(self, idx: Idx2D[str]) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def get(self, adata: AnnData, idx: Idx2D[str], /) -> Vector:
        df = cast("pd.DataFrame", getattr(adata, self.ax))
        iloc = tuple(df.index.get_loc(i) if isinstance(i, str) else i for i in idx)
        return getattr(adata, f"{self.ax}p")[self.k][iloc].toarray()


@dataclass(frozen=True)
class AdAcc[P: AdPath](LayerVecAcc[P]):
    r"""Accessor to create :class:`AdPath`\ s."""

    k: None = field(init=False, default=None)

    layers: LayerAcc[P] = field(init=False)
    obs: MetaVecAcc[P] = field(init=False)
    var: MetaVecAcc[P] = field(init=False)
    obsm: MultiAcc[P] = field(init=False)
    varm: MultiAcc[P] = field(init=False)
    obsp: GraphAcc[P] = field(init=False)
    varp: GraphAcc[P] = field(init=False)

    ATTRS: ClassVar = frozenset({
        "layers",
        "obs",
        "var",
        "obsm",
        "varm",
        "obsp",
        "varp",
    })

    def __post_init__(self) -> None:
        object.__setattr__(self, "layers", LayerAcc(path_class=self.path_class))
        object.__setattr__(self, "obs", MetaVecAcc("obs", path_class=self.path_class))
        object.__setattr__(self, "var", MetaVecAcc("var", path_class=self.path_class))
        object.__setattr__(self, "obsm", MultiAcc("obs", path_class=self.path_class))
        object.__setattr__(self, "varm", MultiAcc("var", path_class=self.path_class))
        object.__setattr__(self, "obsp", GraphAcc("obs", path_class=self.path_class))
        object.__setattr__(self, "varp", GraphAcc("var", path_class=self.path_class))

    @overload
    def resolve(self, spec: str, *, strict: Literal[True] = True) -> P: ...
    @overload
    def resolve(self, spec: str, *, strict: Literal[False]) -> P | None: ...
    def resolve(self, spec: str, *, strict: bool = True) -> P | None:
        """Create accessor from string.

        Examples
        --------
        >>> A.resolve("X[:,:]")
        A[:, :]
        >>> A.resolve("layers.y[c,:]")
        A.layers['y']['c', :]
        >>> A.resolve("layers.y[:,g]")
        A.layers['y'][:, 'g']
        >>> A.resolve("obs.a")
        A.obs['a']
        >>> A.resolve("var.b")
        A.var['b']
        >>> A.resolve("obsm.c.0")
        A.obsm['c'][:, 0]
        >>> A.resolve("varm.d.1")
        A.varm['d'][:, 1]
        >>> A.resolve("obsp.g[c1,:]")
        A.obsp['g']['c1', :]
        >>> A.resolve("obsp.g[:,c2]")
        A.obsp['g'][:, 'c2']
        """
        from ._parse import parse

        return parse(self, spec, strict=strict)

    def __repr__(self) -> str:
        return "A"


A: AdAcc[AdPath] = AdAcc(path_class=AdPath)
r"""A global accessor to create :class:`AdPath`\ s."""


def _is_idx2d_list[Idx: int | str](
    idx: Idx2D[Idx] | Idx2DList[Idx],
) -> TypeIs[Idx2DList[Idx]]:
    """Check if a 2D index contains a list in one of its dimensions."""
    return any(isinstance(i, list) for i in idx)


def _expand_idx2d_list[Idx: int | str](
    idx: Idx2D[Idx] | Idx2DList[Idx],
) -> list[Idx2D[Idx]]:
    """Expand a 2D index containing a list in one of its dimensions.

    Also validates that the 2D index contains at most one list.
    """
    match idx:
        case list(), list():
            msg = "2D index can have at most one list: …[:, [...]] or …[[...], :]"
            raise TypeError(msg)
        case list() as ixs, iy:
            return [(ix, iy) for ix in ixs]
        case ix, list() as iys:
            return [(ix, iy) for iy in iys]
        case _:  # pragma: no cover
            msg = "Should have checked _is_idx2d_list before."
            raise AssertionError(msg)


def _idx2axes(idx: Idx2D[str]) -> set[Literal["obs", "var"]]:
    """Get along which axes the referenced vector is and validate the index."""
    for ax_idx in idx:
        if isinstance(ax_idx, str):
            continue
        if isinstance(ax_idx, slice) and ax_idx == slice(None):
            continue
        msg = (
            f"Unsupported axis index {ax_idx!r} in index {idx!r} (not `:` or a string)"
        )
        raise ValueError(msg)
    match idx:
        case slice(), str():
            return {"obs"}
        case str(), slice():
            return {"var"}
        case slice(), slice():
            return {"obs", "var"}
        case _:  # pragma: no cover
            msg = f"Invalid index: {idx}"
            raise TypeError(msg)


def __getattr__(name: Literal["hv"]) -> Any:
    if name == "hv":
        from . import hv

        return hv

    raise AttributeError(name)
