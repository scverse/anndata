"""Accessor classes for AnnData interface."""

from __future__ import annotations

import abc
import re
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


class AdPath[I]:
    """A path referencing an array in an AnnData object."""

    acc: VecAcc[Self, I]
    idx: I

    def __init__(self, acc: VecAcc[Self, I], idx: I) -> None:
        self.acc = acc
        self.idx = idx

    @cached_property
    def axes(self) -> Axes:
        """Which axes are spanned by the array?"""
        return self.acc.axes(self.idx)

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
        return self.acc(adata, self.idx)


class VecAcc[P: AdPath[I], I](abc.ABC):  # type: ignore
    path_class: type[P]

    def process_idx(self, idx: Any, /) -> I:
        return idx

    def __getitem__(self, idx: Any, /) -> P:
        idx = self.process_idx(idx)
        return self.path_class(self, idx)  # type: ignore

    @abc.abstractmethod
    def __call__(self, adata: AnnData, idx: I, /) -> Vector: ...
    @abc.abstractmethod
    def axes(self, idx: I, /) -> Axes: ...
    @abc.abstractmethod
    def __repr__(self, /) -> str: ...
    @abc.abstractmethod
    def idx_repr(self, idx: I, /) -> str: ...


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
class LayerVecAcc[P: AdPath](VecAcc[P, "Idx2D[str]"]):
    """Accessor for layer vectors."""

    k: str | None
    _: KW_ONLY
    path_class: type[P]

    @overload
    def __getitem__(self, idx: Idx2D[str], /) -> P: ...
    @overload
    def __getitem__(self, idx: Idx2DList[str], /) -> list[P]: ...
    def __getitem__(self, idx: Idx2D[str] | Idx2DList[str], /) -> P | list[P]:
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def __call__(self, adata: AnnData, idx: Idx2D[str], /) -> Vector:
        ver_or_mat = adata[idx].X if self.k is None else adata[idx].layers[self.k]
        if isinstance(ver_or_mat, sp.spmatrix | sp.sparray):
            ver_or_mat = ver_or_mat.toarray()
        # TODO: pandas
        axes = self.axes(idx)
        return ver_or_mat.flatten() if len(axes) == 1 else ver_or_mat

    def axes(self, idx: Idx2D[str], /) -> Axes:
        return _idx2axes(idx)

    def __repr__(self) -> str:
        return f"A.layers[{self.k!r}]"

    def idx_repr(self, idx: Idx2D[str]) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"


@dataclass(frozen=True)
class MetaVecAcc[P: AdPath](VecAcc[P, str | type[pd.Index]]):
    """Accessor for metadata (obs/var)."""

    ax: Literal["obs", "var"]
    _: KW_ONLY
    path_class: type[P]

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

    def __call__(self, adata: AnnData, k: str | type[pd.Index], /) -> Vector:
        df = cast("pd.DataFrame", getattr(adata, self.ax))
        return df.index.values if k is pd.Index else df[k].values

    def axes(self, k: str, /) -> Axes:
        return {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}"

    def idx_repr(self, k: str | type[pd.Index]) -> str:
        return ".index" if k is pd.Index else f"[{k!r}]"


@dataclass(frozen=True)
class MultiAcc[P: AdPath]:
    """Accessor for multi-dimensional array containers (obsm/varm)."""

    ax: Literal["obsm", "varm"]
    _: KW_ONLY
    path_class: type[P]

    def __getitem__(self, k: str, /) -> MultiVecAcc[P]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax} key {k!r}"
            raise TypeError(msg)
        return MultiVecAcc(self.ax, k, path_class=self.path_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}"


@dataclass(frozen=True)
class MultiVecAcc[P: AdPath](VecAcc[P, int]):
    """Accessor for arrays from multi-dimensional containers (obsm/varm)."""

    ax: Literal["obsm", "varm"]
    k: str
    _: KW_ONLY
    path_class: type[P]

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

    def __call__(self, adata: AnnData, i: int, /) -> Vector:
        return getattr(adata, self.ax)[self.k][:, i]

    def axes(self, i: int, /) -> Axes:
        return {cast("Literal['obs', 'var']", self.ax[:-1])}

    def __repr__(self) -> str:
        return f"A.{self.ax}[{self.k!r}]"

    def idx_repr(self, i: int) -> str:
        return f"[:, {i!r}]"


@dataclass(frozen=True)
class GraphAcc[P: AdPath]:
    """Accessor for graph containers (obsp/varp)."""

    ax: Literal["obsp", "varp"]
    _: KW_ONLY
    path_class: type[P]

    def __getitem__(self, k: str, /) -> GraphVecAcc[P]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax} key {k!r}"
            raise TypeError(msg)
        return GraphVecAcc(self.ax, k, path_class=self.path_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}"


@dataclass(frozen=True)
class GraphVecAcc[P: AdPath](VecAcc[P, "Idx2D[str]"]):
    """Accessor for arrays from graph containers (obsp/varp)."""

    ax: Literal["obsp", "varp"]
    k: str
    _: KW_ONLY
    path_class: type[P]

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

    def __call__(self, adata: AnnData, idx: Idx2D[str], /) -> Vector:
        df = cast("pd.DataFrame", getattr(adata, self.ax[:-1]))
        iloc = tuple(df.index.get_loc(i) if isinstance(i, str) else i for i in idx)
        return getattr(adata, self.ax)[self.k][iloc].toarray()

    def axes(self, idx: Idx2D[str], /) -> Axes:
        n_slices = sum(isinstance(i, slice) for i in idx)
        ax = cast("Literal['obs', 'var']", self.ax[:-1])
        return (ax,) * n_slices if n_slices > 1 else {ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}[{self.k!r}]"

    def idx_repr(self, idx: Idx2D[str]) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"


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
        object.__setattr__(self, "obsm", MultiAcc("obsm", path_class=self.path_class))
        object.__setattr__(self, "varm", MultiAcc("varm", path_class=self.path_class))
        object.__setattr__(self, "obsp", GraphAcc("obsp", path_class=self.path_class))
        object.__setattr__(self, "varp", GraphAcc("varp", path_class=self.path_class))

    @overload
    def resolve(self, spec: str, *, strict: Literal[True] = True) -> P: ...
    @overload
    def resolve(self, spec: str, *, strict: Literal[False]) -> P | None: ...
    def resolve(self, spec: str, *, strict: bool = True) -> P | None:
        """Create accessor from string."""
        if not strict:
            try:
                self.resolve(spec)
            except ValueError:
                return None

        if "." not in spec:
            msg = f"Cannot parse accessor {spec!r}"
            raise ValueError(msg)
        acc, rest = spec.split(".", 1)
        match getattr(self, acc, None):
            # TODO: X
            case LayerAcc() as layers:
                return self._parse_path_layer(layers, rest)
            case MetaVecAcc() as meta:
                return meta[rest]
            case MultiAcc() as multi:
                return self._parse_path_multi(multi, rest)
            case GraphAcc():
                msg = "TODO"
                raise NotImplementedError(msg)
            case None:  # pragma: no cover
                msg = (
                    f"Unknown accessor {spec!r}. "
                    f"We support '{self.ATTRS}.*' and `AdPath` instances."
                )
                raise ValueError(msg)
        msg = f"Unhandled accessor {spec!r}. This is a bug!"  # pragma: no cover
        raise AssertionError(msg)  # pragma: no cover

    def __repr__(self) -> str:
        return "A"

    def _parse_path_layer(self, layers: LayerAcc, spec: str) -> P:
        if not (
            m := re.fullmatch(r"([^\[]+)\[([^,]+),\s?([^\]]+)\]", spec)
        ):  # pragma: no cover
            msg = f"Cannot parse layer accessor {spec!r}: should be `name[i,:]`/`name[:,j]`"
            raise ValueError(msg)
        layer, i, j = m.groups("")  # "" just for typing
        return layers[layer][_parse_idx_2d(i, j, str)]

    def _parse_path_multi(self, multi: MultiAcc, spec: str) -> P:
        if not (m := re.fullmatch(r"([^.]+)\.([\d_]+)", spec)):  # pragma: no cover
            msg = f"Cannot parse multi accessor {spec!r}: should be `name.i`"
            raise ValueError(msg)
        key, i = m.groups("")  # "" just for typing
        return multi[key][int(i)]


A = AdAcc(path_class=AdPath)


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


def _parse_idx_2d[Idx: int | str](i: str, j: str, cls: type[Idx]) -> Idx2D[Idx]:
    match i, j:
        case _, ":":
            return cls(i), slice(None)
        case ":", _:
            return slice(None), cls(j)
        case _:  # pragma: no cover
            msg = f"Unknown indices {i!r}, {j!r}"
            raise ValueError(msg)


def __getattr__(name: Literal["hv"]) -> Any:
    if name == "hv":
        from . import hv

        return hv

    raise AttributeError(name)
