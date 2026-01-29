"""Accessor classes for AnnData interface."""

from __future__ import annotations

import abc
from collections.abc import Hashable
from dataclasses import KW_ONLY, dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING, ClassVar, cast, overload

import scipy.sparse as sp
from numpy.typing import NDArray
from pandas.api.extensions import ExtensionArray

from anndata import AnnData

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Sequence
    from typing import Any, Literal, Self, TypeIs

    import pandas as pd
    from numpy.typing import NDArray

    from anndata import AnnData

    if TYPE_CHECKING:  # for sphinx
        from . import hv

    type Axes = Collection[Literal["obs", "var"]]


type Array = ExtensionArray | NDArray[Any]

type Idx2D = tuple[str | slice, slice] | tuple[slice, str | slice]
"""Index along the full length of one or both AnnData axes, resulting in a 1D or 2D array.

E.g. `a[:, 5]`, `a[18, :]`, or `a[:, :]`
"""

type Idx2DList = tuple[list[str], slice] | tuple[slice, list[str]]
type AdRefFunc[I] = Callable[[AnnData, I], Array]


__all__ = [
    "A",
    "AdAcc",
    "AdRef",
    "GraphAcc",
    "GraphMapAcc",
    "Idx2D",
    "LayerAcc",
    "LayerMapAcc",
    "MetaAcc",
    "MultiAcc",
    "MultiMapAcc",
    "RefAcc",
    "hv",
]


class AdRef[I: Hashable]:
    """A reference to a 1D or 2D array along the axes of an AnnData object.

    Examples
    --------
    Use :data:`A` in the same way you would an :class:`~anndata.AnnData` object if you wanted to access a 1D or 2D array along the whole axes.
    You can then introspect the reference:

    >>> from anndata.acc import A
    >>> A.obs.index.axes
    {'obs'}
    >>> A.obs["x"].idx
    'x'
    >>> type(A.var["y"].acc)
    <class 'anndata.acc.MetaAcc'>

    See :mod:`anndata.acc` and :class:`AdAcc` for more examples.

    """

    acc: RefAcc[Self, I]
    """The accessor containing information about this array.

    See :ref:`here <reference-accessors>` for all possible types this can assume.
    """

    idx: I
    """The index specifying the array."""

    __match_args__: ClassVar = ("acc", "idx")

    def __init__(self, acc: RefAcc[Self, I], idx: I) -> None:
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
            case AdRef():
                return self.acc == value.acc and self.idx == value.idx
            case str():
                # TODO: more elegant: access `A` through `acc`?
                a = (
                    A
                    if type(self).__eq__ is AdRef.__eq__
                    else AdAcc(ref_class=type(self))
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

    def __call__(self, adata: AnnData) -> Array:
        """Retrieve referenced array from AnnData."""
        return self.acc.get(adata, self.idx)


@dataclass(frozen=True)
class RefAcc[R: AdRef[I], I](abc.ABC):  # type: ignore
    """Abstract base class for reference accessors.

    See :ref:`here <reference-accessors>` for all existing subclasses.
    """

    _: KW_ONLY
    ref_class: type[R]

    def process_idx(self, idx: Any, /) -> I:
        self.axes(idx)
        return idx

    def __getitem__(self, idx: Any, /) -> R:
        idx = self.process_idx(idx)
        return self.ref_class(self, idx)  # type: ignore

    @abc.abstractmethod
    def axes(self, idx: I, /) -> Axes: ...
    @abc.abstractmethod
    def __repr__(self, /) -> str: ...
    @abc.abstractmethod
    def idx_repr(self, idx: I, /) -> str: ...
    @abc.abstractmethod
    def get(self, adata: AnnData, idx: I, /) -> Array: ...


@dataclass(frozen=True)
class LayerMapAcc[R: AdRef]:
    r"""Accessor for layers (`A.`\ :attr:`~AdAcc.layers`)."""

    _: KW_ONLY
    ref_class: type[R]

    def __getitem__(self, k: str, /) -> LayerAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported layer {k!r}"
            raise TypeError(msg)
        return LayerAcc(k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return "A.layers"


@dataclass(frozen=True)
class LayerAcc[R: AdRef[Idx2D]](RefAcc[R, Idx2D]):
    r"""Reference accessor for arrays in layers (`A.`\ :attr:`~AdAcc.layers`).

    Examples
    --------
    >>> from anndata.acc import A, LayerAcc
    >>> assert isinstance(A.layers["counts"], LayerAcc)
    >>> A.layers["counts"]["cell-1", :]
    A.layers['counts']['cell-1', :]

    """

    k: str | None
    """Key this accessor refers to, e.g. `A.layers['counts'].k == 'counts'`."""

    @overload
    def __getitem__(self, idx: Idx2D, /) -> R: ...
    @overload
    def __getitem__(self, idx: Idx2DList, /) -> list[R]: ...
    def __getitem__(self, idx: Idx2D | Idx2DList, /) -> R | list[R]:
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def axes(self, idx: Idx2D, /) -> Axes:
        return _idx2axes(idx)

    def __repr__(self) -> str:
        return f"A.layers[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def get(self, adata: AnnData, idx: Idx2D, /) -> Array:
        ver_or_mat = adata[idx].X if self.k is None else adata[idx].layers[self.k]
        if isinstance(ver_or_mat, sp.spmatrix | sp.sparray):
            ver_or_mat = ver_or_mat.toarray()
        # TODO: pandas
        axes = self.axes(idx)
        return ver_or_mat.flatten() if len(axes) == 1 else ver_or_mat


@dataclass(frozen=True)
class MetaAcc[R: AdRef[str | None]](RefAcc[R, str | None]):
    r"""Reference accessor for arrays from metadata containers (`A.`\ :attr:`~AdAcc.obs`/`A.`\ :attr:`~AdAcc.var`).

    Examples
    --------
    You can refer to columns or the index in the way you would access them
    in a :class:`~pandas.DataFrame`:

    >>> from anndata.acc import A, MetaAcc
    >>> assert isinstance(A.obs, MetaAcc)
    >>> A.obs["type"]
    A.obs['type']
    >>> A.var.index
    A.var.index

    """

    ax: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obs.ax == 'obs'`."""

    @property
    def index(self) -> R:
        """Index :class:`AdRef`, i.e. `A.obs.index` or `A.var.index`."""
        return self[None]

    def process_idx(self, k: object, /) -> str | None:
        if k is not None and not isinstance(k, str):
            msg = f"Unsupported {self.ax} column {k!r}"
            raise TypeError(msg)
        return k

    @overload
    def __getitem__[K: str | None](self, k: K, /) -> R: ...
    @overload
    def __getitem__[K: str | None](self, k: list[K], /) -> list[R]: ...
    def __getitem__[K: str | None](self, k: K | list[K], /) -> R | list[R]:
        if isinstance(k, list):
            return [self[i] for i in k]

        return super().__getitem__(k)

    def axes(self, k: str, /) -> Axes:
        return {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}"

    def idx_repr(self, k: str | None) -> str:
        return ".index" if k is None else f"[{k!r}]"

    def get(self, adata: AnnData, k: str | None, /) -> Array:
        df = cast("pd.DataFrame", getattr(adata, self.ax))
        return df.index.values if k is None else df[k].values


@dataclass(frozen=True)
class MultiMapAcc[R: AdRef]:
    r"""Accessor for multi-dimensional array containers (`A.`\ :attr:`~AdAcc.obsm`/`A.`\ :attr:`~AdAcc.varm`)."""

    ax: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsm.ax == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]

    def __getitem__(self, k: str, /) -> MultiAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax}m key {k!r}"
            raise TypeError(msg)
        return MultiAcc(self.ax, k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}m"


@dataclass(frozen=True)
class MultiAcc[R: AdRef[int]](RefAcc[R, int]):
    r"""Reference accessor for arrays from multi-dimensional containers (`A.`\ :attr:`~AdAcc.obsm`/`A.`\ :attr:`~AdAcc.varm`).

    Examples
    --------
    >>> from anndata.acc import A, MultiAcc
    >>> assert isinstance(A.obsm["pca"], MultiAcc)
    >>> A.obsm["pca"][:, 0]
    A.obsm['pca'][:, 0]

    """

    ax: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsm[k].ax == 'obs'`."""

    k: str
    """Key this accessor refers to, e.g. `A.varm['x'].k == 'x'`."""

    @staticmethod
    def process_idx[T](i: T | tuple[slice, T], /) -> T:
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
    def __getitem__(self, i: int | tuple[slice, int], /) -> R: ...
    @overload
    def __getitem__(self, i: list[int] | tuple[slice, list[int]], /) -> list[R]: ...
    def __getitem__(
        self, i: int | list[int] | tuple[slice, int | list[int]], /
    ) -> R | list[R]:
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

    def get(self, adata: AnnData, i: int, /) -> Array:
        return getattr(adata, f"{self.ax}m")[self.k][:, i]


@dataclass(frozen=True)
class GraphMapAcc[R: AdRef]:
    r"""Accessor for graph containers (`A.`\ :attr:`~AdAcc.obsp`/`A.`\ :attr:`~AdAcc.varp`)."""

    ax: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsp.ax == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]

    def __getitem__(self, k: str, /) -> GraphAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.ax}p key {k!r}"
            raise TypeError(msg)
        return GraphAcc(self.ax, k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return f"A.{self.ax}p"


@dataclass(frozen=True)
class GraphAcc[R: AdRef[Idx2D]](RefAcc[R, Idx2D]):
    r"""Reference accessor for arrays from graph containers (`A.`\ :attr:`~AdAcc.obsp`/`A.`\ :attr:`~AdAcc.varp`).

    Examples
    --------
    >>> from anndata.acc import A, GraphAcc
    >>> assert isinstance(A.obsp["distances"], GraphAcc)
    >>> A.obsp["distances"][:, "cell-1"]
    A.obsp['distances'][:, 'cell-1']

    """

    ax: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.varp[k].axis == "var"`."""

    k: str
    """Key this accessor refers to, e.g. `A.obsp['x'].k == 'x'`."""

    def process_idx(self, idx: Idx2D, /) -> Idx2D:
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
    def __getitem__(self, idx: Idx2D, /) -> R: ...
    @overload
    def __getitem__(self, idx: Idx2DList, /) -> list[R]: ...
    def __getitem__(self, idx: Idx2D | Idx2DList, /) -> R | list[R]:
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def axes(self, idx: Idx2D, /) -> Axes:
        n_slices = sum(isinstance(i, slice) for i in idx)
        return (self.ax,) * n_slices if n_slices > 1 else {self.ax}

    def __repr__(self) -> str:
        return f"A.{self.ax}p[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def get(self, adata: AnnData, idx: Idx2D, /) -> Array:
        df = cast("pd.DataFrame", getattr(adata, self.ax))
        iloc = tuple(df.index.get_loc(i) if isinstance(i, str) else i for i in idx)
        return getattr(adata, f"{self.ax}p")[self.k][iloc].toarray()


@dataclass(frozen=True)
class AdAcc[R: AdRef](LayerAcc[R]):
    r"""Accessor to create :class:`AdRef`\ s (:data:`A`).

    See examples below and in :mod:`anndata.acc`.
    """

    k: None = field(init=False, default=None)
    """: `A[:, :]` is equivalent to `A.layers[None][:, :]`."""  # weird “:” is load-bearing 🤷

    layers: LayerMapAcc[R] = field(init=False)
    """Access complete layers or 1D vectors across observations or variables.

    >>> A.layers["counts"][:, :]
    >>> A.layers["counts"]["cell-1", :]
    >>> A.layers["counts"][:, "gene-5"]
    """

    obs: MetaAcc[R] = field(init=False)
    """Access 1D arrays along observations.

    >>> A.obs["cell-type"]
    >>> A.obs.index
    """

    var: MetaAcc[R] = field(init=False)
    """Access 1D arrays along variables.

    >>> A.var["symbols"]
    >>> A.obs.index
    """

    obsm: MultiMapAcc[R] = field(init=False)
    """Access 1D vectors along observations.

    >>> A.obsm["pca"][:, 0].idx
    0
    """

    varm: MultiMapAcc[R] = field(init=False)
    """Access 1D vectors along variables.

    >>> A.varm["loadings"][:, 0].idx
    0
    """

    obsp: GraphMapAcc[R] = field(init=False)
    """Access 1D or 2D vectors along observations.

    >>> A.layers["x"][:, :].axes
    ('obs', 'obs')
    >>> A.layers["x"]["cell-1", :].axes
    {'obs'}
    >>> A.layers["x"][:, "cell-1"].axes
    {'obs'}
    """

    varp: GraphMapAcc[R] = field(init=False)
    """Access 1D or 2D vectors along variables.

    >>> A.layers["x"][:, :].axes
    ('var', 'var')
    >>> A.layers["x"]["gene-1", :].axes
    {'var'}
    >>> A.layers["x"][:, "gene-1"].axes
    {'var'}
    """

    ATTRS: ClassVar = frozenset({
        *("layers", "obs", "var"),
        *("obsm", "varm", "obsp", "varp"),
    })

    def __post_init__(self) -> None:
        object.__setattr__(self, "layers", LayerMapAcc(ref_class=self.ref_class))
        object.__setattr__(self, "obs", MetaAcc("obs", ref_class=self.ref_class))
        object.__setattr__(self, "var", MetaAcc("var", ref_class=self.ref_class))
        object.__setattr__(self, "obsm", MultiMapAcc("obs", ref_class=self.ref_class))
        object.__setattr__(self, "varm", MultiMapAcc("var", ref_class=self.ref_class))
        object.__setattr__(self, "obsp", GraphMapAcc("obs", ref_class=self.ref_class))
        object.__setattr__(self, "varp", GraphMapAcc("var", ref_class=self.ref_class))

    def to_json(self, ref: R) -> list[str | int | None]:
        """Serialize :class:`AdRef` to a JSON-compatible list."""
        from ._parse_json import to_json

        return to_json(ref)

    def from_json(self, data: Sequence[str | int | None]) -> R:
        """Create :class:`AdRef` from a JSON sequence.

        Raises
        ------
        ValueError
            If parsing fails.
        """
        from ._parse_json import parse_json

        try:
            return parse_json(self, data)
        except Exception as e:
            msg = f"Failed to parse {data!r}"
            raise ValueError(msg) from e

    @overload
    def resolve(self, spec: str, *, strict: Literal[True] = True) -> R: ...
    @overload
    def resolve(self, spec: str, *, strict: Literal[False]) -> R | None: ...
    def resolve(self, spec: str, *, strict: bool = True) -> R | None:
        """Create :class:`AdRef` from a simplified string.

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
        from ._parse_str import parse

        return parse(self, spec, strict=strict)

    def __repr__(self) -> str:
        return "A"


A: AdAcc[AdRef] = AdAcc(ref_class=AdRef)
r"""A global accessor to create :class:`AdRef`\ s."""


def _is_idx2d_list(idx: Idx2D | Idx2DList) -> TypeIs[Idx2DList]:
    """Check if a 2D index contains a list in one of its dimensions."""
    return any(isinstance(i, list) for i in idx)


def _expand_idx2d_list(idx: Idx2D | Idx2DList) -> list[Idx2D]:
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


def _idx2axes(idx: Idx2D) -> set[Literal["obs", "var"]]:
    """Get along which axes the referenced array is and validate the index."""
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
    __tracebackhide__ = True

    if name == "hv":
        import anndata.acc.hv

        return anndata.acc.hv

    raise AttributeError(name)


if not TYPE_CHECKING:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    R = AdRef[Hashable]
    I = Hashable
