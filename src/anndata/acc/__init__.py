"""Accessor classes for AnnData interface."""

from __future__ import annotations

import abc
from collections.abc import Hashable
from dataclasses import KW_ONLY, dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING, ClassVar, cast, overload

import pandas as pd
import scipy.sparse as sp

from anndata import AnnData

from .._core.views import ArrayView
from .._core.xarray import Dataset2D
from ..compat import DaskArray, has_xp

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Sequence
    from typing import Any, Literal, Self, TypeGuard

    from .._core.aligned_mapping import AxisArrays
    from ..compat import XVariable
    from ..typing import InMemoryArray

    if TYPE_CHECKING:  # for sphinx
        from . import hv


type Axes = Collection[Literal["obs", "var"]]

type Idx2D = tuple[str | slice, slice] | tuple[slice, str | slice]
"""Index along the full length of one or both AnnData dimensions, resulting in a 1D or 2D array.

E.g. `a[:, "c"]`, `a["b", :]`, or `a[:, :]`
"""

type Idx2DList = (
    tuple[list[str] | pd.Index[str], slice] | tuple[slice, list[str] | pd.Index[str]]
)
type IdxMultiList = list[int] | pd.Index[int] | tuple[slice, list[int] | pd.Index[int]]
type AdRefFunc[I] = Callable[[AnnData, I], InMemoryArray]


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
    r"""A reference to a 1D or 2D array along one or two dimensions of an AnnData object.

    Examples
    --------
    Use :data:`A` in the same way you would an :class:`~anndata.AnnData` object if you wanted to access a 1D or 2D array along the whole dimension.
    You can then introspect the reference:

    >>> from anndata.acc import A
    >>> ref = A.obs["x"]
    >>> ref.dims
    {'obs'}
    >>> ref.idx
    'x'
    >>> type(ref.acc)
    <class 'anndata.acc.MetaAcc'>

    See :mod:`anndata.acc` and :class:`AdAcc` for more examples.

    """

    acc: RefAcc[Self, I]
    """The accessor containing information about this array.

    See :ref:`here <reference accessors>` for all possible types this can assume.
    """

    idx: I
    """The index specifying the array."""

    __match_args__: ClassVar = ("acc", "idx")

    def __init__(self, acc: RefAcc[Self, I], idx: I) -> None:
        self.acc = acc
        self.idx = idx

    @cached_property
    def dims(self) -> Axes:
        """Which dimensions are spanned by the array?"""
        return self.acc.dims(self.idx)

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


class MapAcc[R: RefAcc](abc.ABC):
    r"""Accessor for mapping containers."""

    ref_acc_cls: type[R]

    @abc.abstractmethod
    def __getitem__(self, k: str, /) -> R:
        """Get a reference accessor for mapped array."""


@dataclass(frozen=True)
class RefAcc[R: AdRef[I], I](abc.ABC):  # type: ignore
    """Abstract base class for reference accessors.

    See :ref:`here <reference accessors>` for all existing subclasses.
    """

    _: KW_ONLY
    ref_class: type[R]

    def process_idx(self, idx: Any, /) -> I:
        self.dims(idx)
        return idx

    def __getitem__(self, idx: Any, /) -> R:
        idx = self.process_idx(idx)
        return self.ref_class(self, idx)  # type: ignore

    @abc.abstractmethod
    def dims(self, idx: I, /) -> Axes:
        """Get along which dimensions the referenced array is."""

    @abc.abstractmethod
    def __repr__(self, /) -> str:
        """Get a string representation of the accessor."""

    @abc.abstractmethod
    def idx_repr(self, idx: I, /) -> str:
        """Get a string representation of the index."""

    @abc.abstractmethod
    def isin(self, adata: AnnData, idx: I | None = None, /) -> bool:
        """Check if the referenced array is in the AnnData object."""

    @abc.abstractmethod
    def get(self, adata: AnnData, idx: I, /) -> InMemoryArray:
        """Get the referenced array from the AnnData object."""

    def _maybe_flatten(self, idx: I, a: InMemoryArray) -> InMemoryArray:
        if len(self.dims(idx)) != 1:
            return a
        if isinstance(a, DaskArray):
            a = a.map_blocks(lambda x: self._maybe_flatten(idx, x))
        if isinstance(a, sp.sparray | sp.spmatrix | ArrayView):
            a = a.toarray()
        if has_xp(a):
            return a.__array_namespace__().reshape(a, (a.size,))
        return a.ravel()


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

    def dims(self, idx: Idx2D, /) -> set[Literal["obs", "var"]]:
        for dim_idx in idx:
            if isinstance(dim_idx, str):
                continue
            if isinstance(dim_idx, slice) and dim_idx == slice(None):
                continue
            msg = f"Unsupported dimension index {dim_idx!r} in index {idx!r} (not `:` or a string)"
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

    def __repr__(self) -> str:
        return f"A.layers[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def isin(self, adata: AnnData, idx: Idx2D | None = None) -> bool:
        if adata.X is None if self.k is None else self.k not in adata.layers:
            return False
        for i, dim in zip(idx or (), ("obs", "var"), strict=False):
            if isinstance(i, str):  # at most one str
                return i in getattr(adata, dim).index
        return True  # idx is None or [:, :]

    def get(self, adata: AnnData, idx: Idx2D, /) -> InMemoryArray:
        arr = adata[idx].X if self.k is None else adata[idx].layers[self.k]
        return self._maybe_flatten(idx, arr)


@dataclass(frozen=True)
class LayerMapAcc[R: AdRef](MapAcc[LayerAcc]):
    r"""Accessor for layers (`A.`\ :attr:`~AdAcc.layers`)."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[LayerAcc] = LayerAcc

    def __getitem__(self, k: str, /) -> LayerAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported layer {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return "A.layers"


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

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obs.dim == 'obs'`."""

    @property
    def index(self) -> R:
        """Index :class:`AdRef`, i.e. `A.obs.index` or `A.var.index`."""
        return self[None]

    def process_idx(self, k: object, /) -> str | None:
        if k is not None and not isinstance(k, str):
            msg = f"Unsupported {self.dim} column {k!r}"
            raise TypeError(msg)
        return k

    @overload
    def __getitem__[K: str | None](self, k: K, /) -> R: ...
    @overload
    def __getitem__[K: str | None](self, k: list[K] | pd.Index[str], /) -> list[R]: ...
    def __getitem__[K: str | None](
        self, k: K | list[K] | pd.Index[str], /
    ) -> R | list[R]:
        # _is_t_list doesn’t allow None, so be more liberal when checking for lists
        if isinstance(k, list) or _is_t_list(k, str):
            return [self[i] for i in k]

        return super().__getitem__(k)

    def dims(self, k: str, /) -> Axes:
        return {self.dim}

    def __repr__(self) -> str:
        return f"A.{self.dim}"

    def idx_repr(self, k: str | None) -> str:
        return ".index" if k is None else f"[{k!r}]"

    def isin(self, adata: AnnData, idx: str | None = None) -> bool:
        if idx is None:
            return True  # obs and var index always exist
        attr: pd.DataFrame | Dataset2D = getattr(adata, self.dim)
        return idx in attr

    def get(
        self, adata: AnnData, k: str | None, /
    ) -> pd.api.extensions.ExtensionArray | XVariable:
        match getattr(adata, self.dim), k:
            case pd.DataFrame() as df, None:
                return df.index.array
            case Dataset2D() as ds, None:
                return ds.xr_index.variable  # TODO: return cached index instead?
            case pd.DataFrame() as df, k:
                return df[k].array
            case Dataset2D() as ds, k:
                return ds[k].variable
            case _:
                msg = f"Unsupported {self.dim} container"
                raise TypeError(msg)


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

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsm[k].dim == 'obs'`."""

    k: str
    """Key this accessor refers to, e.g. `A.varm['x'].k == 'x'`."""

    @staticmethod
    def process_idx(i: object, /) -> int | list[int] | pd.Index[int]:
        if isinstance(i, tuple):
            if len(i) != 2 or i[0] != slice(None):
                msg = f"Unsupported slice {i!r}"
                raise ValueError(msg)
            i = i[1]
        if not isinstance(i, int) and not _is_t_list(i, int):
            msg = f"Unsupported index {i!r}"
            raise TypeError(msg)
        return i

    @overload
    def __getitem__(self, i: int | tuple[slice, int], /) -> R: ...
    @overload
    def __getitem__(self, i: IdxMultiList, /) -> list[R]: ...
    def __getitem__(self, i: int | tuple[slice, int] | IdxMultiList, /) -> R | list[R]:
        i = self.process_idx(i)
        if isinstance(i, list | pd.Index):
            return [self[j] for j in i]
        return super().__getitem__(i)

    def dims(self, i: int, /) -> Axes:
        return {self.dim}

    def __repr__(self) -> str:
        return f"A.{self.dim}m[{self.k!r}]"

    def idx_repr(self, i: int) -> str:
        return f"[:, {i!r}]"

    def isin(self, adata: AnnData, idx: int | None = None) -> bool:
        m: AxisArrays = getattr(adata, f"{self.dim}m")
        if self.k not in m:
            return False
        return idx is None or idx in range(m[self.k].shape[1])

    def get(self, adata: AnnData, i: int, /) -> InMemoryArray:
        # TODO: remove slicing when dropping scipy <1.14
        arr = getattr(adata, f"{self.dim}m")[self.k][:, i : i + 1]
        return self._maybe_flatten(i, arr)


@dataclass(frozen=True)
class MultiMapAcc[R: AdRef](MapAcc[MultiAcc]):
    r"""Accessor for multi-dimensional array containers (`A.`\ :attr:`~AdAcc.obsm`/`A.`\ :attr:`~AdAcc.varm`)."""

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsm.dim == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[MultiAcc] = MultiAcc

    def __getitem__(self, k: str, /) -> MultiAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.dim}m key {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(self.dim, k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return f"A.{self.dim}m"


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

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.varp[k].dim == "var"`."""

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

    def dims(self, idx: Idx2D, /) -> Axes:
        n_slices = sum(isinstance(i, slice) for i in idx)
        return (self.dim,) * n_slices if n_slices > 1 else {self.dim}

    def __repr__(self) -> str:
        return f"A.{self.dim}p[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def isin(self, adata: AnnData, idx: Idx2D | None = None) -> bool:
        if self.k not in getattr(adata, f"{self.dim}p"):
            return False
        if idx is None:
            return True
        [i] = (i for i in idx if isinstance(i, str))
        return i in getattr(adata, self.dim).index

    def get(self, adata: AnnData, idx: Idx2D, /) -> InMemoryArray:
        df = cast("pd.DataFrame", getattr(adata, self.dim))
        # TODO: remove wrapping in [] when dropping scipy <1.14
        iloc = tuple([df.index.get_loc(i)] if isinstance(i, str) else i for i in idx)
        arr = getattr(adata, f"{self.dim}p")[self.k][iloc]
        return self._maybe_flatten(idx, arr)


@dataclass(frozen=True)
class GraphMapAcc[R: AdRef](MapAcc[GraphAcc]):
    r"""Accessor for graph containers (`A.`\ :attr:`~AdAcc.obsp`/`A.`\ :attr:`~AdAcc.varp`)."""

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsp.dim == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[GraphAcc] = GraphAcc

    def __getitem__(self, k: str, /) -> GraphAcc[R]:
        if not isinstance(k, str):
            msg = f"Unsupported {self.dim}p key {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(self.dim, k, ref_class=self.ref_class)

    def __repr__(self) -> str:
        return f"A.{self.dim}p"


@dataclass(frozen=True)
class AdAcc[R: AdRef](LayerAcc[R]):
    r"""Accessor to create :class:`AdRef`\ s (:data:`A`).

    See examples below and in :mod:`anndata.acc`.
    """

    k: None = field(init=False, default=None)
    """: `A[:, :]` is equivalent to `A.layers[None][:, :]`."""  # weird “:” is load-bearing 🤷

    ref_class: type[R] = AdRef

    layer_cls: type[LayerAcc] = LayerAcc
    """Class to use for `layers` accessors.

    Note that :class:`!AdAcc` inherits from :class:`!LayerAcc`,
    so if you want `A[:, :]` to inherit the behavior,
    you need to create and use an :class:`!AdAcc` subclass.
    """

    meta_cls: type[MetaAcc] = MetaAcc
    """Class to use for `obs`/`var` accessors."""

    multi_cls: type[MultiAcc] = MultiAcc
    """Class to use for `obsm`/`varm` accessors."""

    graph_cls: type[GraphAcc] = GraphAcc
    """Class to use for `obsp`/`varp` accessors."""

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

    >>> A.layers["x"][:, :].dims
    ('obs', 'obs')
    >>> A.layers["x"]["cell-1", :].dims
    {'obs'}
    >>> A.layers["x"][:, "cell-1"].dims
    {'obs'}
    """

    varp: GraphMapAcc[R] = field(init=False)
    """Access 1D or 2D vectors along variables.

    >>> A.layers["x"][:, :].dims
    ('var', 'var')
    >>> A.layers["x"]["gene-1", :].dims
    {'var'}
    >>> A.layers["x"][:, "gene-1"].dims
    {'var'}
    """

    ATTRS: ClassVar = frozenset({
        *("layers", "obs", "var"),
        *("obsm", "varm", "obsp", "varp"),
    })

    def __post_init__(self) -> None:
        layers = LayerMapAcc(ref_class=self.ref_class, ref_acc_cls=self.layer_cls)
        object.__setattr__(self, "layers", layers)
        for dim in ("obs", "var"):
            meta = self.meta_cls(dim, ref_class=self.ref_class)
            multi = MultiMapAcc(
                dim, ref_class=self.ref_class, ref_acc_cls=self.multi_cls
            )
            graphs = GraphMapAcc(
                dim, ref_class=self.ref_class, ref_acc_cls=self.graph_cls
            )
            object.__setattr__(self, dim, meta)
            object.__setattr__(self, f"{dim}m", multi)
            object.__setattr__(self, f"{dim}p", graphs)

    def to_json(self, ref: R) -> list[str | int | None]:
        """Serialize :class:`AdRef` to a JSON-compatible list.

        Schema: `acc-schema-v1.json <../acc-schema-v1.json>`_
        """
        from ._parse_json import to_json

        return to_json(ref)

    def from_json(self, data: Sequence[str | int | None]) -> R:
        """Create :class:`AdRef` from a JSON sequence.

        Schema: `acc-schema-v1.json <../acc-schema-v1.json>`_

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


A: AdAcc[AdRef] = AdAcc()
r"""A global accessor to create :class:`AdRef`\ s."""


_checks: dict[type[int | str], Callable[..., bool]] = {
    str: pd.api.types.is_string_dtype,
    int: pd.api.types.is_integer_dtype,
}


def _is_t_list[T: (int, str)](
    idx: object, cls: type[T], /
) -> TypeGuard[list[T] | pd.Index[T]]:
    if isinstance(idx, pd.Index) and _checks[cls](idx.dtype):
        return True
    return isinstance(idx, list | pd.Index) and all(isinstance(j, cls) for j in idx)


def _is_idx2d_list(idx: Idx2D | Idx2DList) -> TypeGuard[Idx2DList]:
    """Check if a 2D index contains a list in one of its dimensions."""
    return any(_is_t_list(i, str) for i in idx)


def _expand_idx2d_list(idx: Idx2DList) -> list[Idx2D]:
    """Expand a 2D index containing a string list or pandas index in one of its dimensions.

    Also validates that the 2D index contains at most one list.
    """
    match idx:
        case list() | pd.Index(), list() | pd.Index():
            msg = "2D index can have at most one list/index: …[:, [...]] or …[[...], :]"
            raise TypeError(msg)
        case list() | pd.Index() as ixs, iy:
            return [(ix, iy) for ix in ixs]
        case ix, list() | pd.Index() as iys:
            return [(ix, iy) for iy in iys]
        case _:  # pragma: no cover
            msg = "Should have checked _is_idx2d_list before."
            raise AssertionError(msg)


def __getattr__(name: Literal["hv"]) -> Any:
    __tracebackhide__ = True

    if name == "hv":
        import anndata.acc.hv

        return anndata.acc.hv

    raise AttributeError(name)


if not TYPE_CHECKING:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    R = AdRef[Hashable]
    I = Hashable
