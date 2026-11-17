"""Accessor classes for AnnData interface."""

from __future__ import annotations

import abc
import sys
from collections.abc import Hashable
from dataclasses import KW_ONLY, dataclass, field, replace
from functools import cached_property
from typing import (
    TYPE_CHECKING,
    ClassVar,
    Never,
    Protocol,
    overload,
    runtime_checkable,
)

if sys.version_info < (3, 15):
    from typing_extensions import sentinel

import numpy as np
import pandas as pd
from scverse_misc import Deprecation, future_abstractmethod

from .. import AnnData
from .._core.views import ArrayView
from .._core.xarray import Dataset2D
from ..compat import (
    CSArray,
    CSMatrix,
    CupySparseMatrix,
    DaskArray,
    SupportsArrayApiBase,
    XVariable,
)
from ..typing import AlignedArray, InMemoryArray, _is_in_memory

if TYPE_CHECKING:
    from typing import TypeAlias

    from mudata import MuData
else:
    MuData = type("MuData", (), {"__module__": "mudata"})

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, MutableMapping, Sequence
    from typing import Any, Literal, Self, TypeGuard

    from .._core.aligned_mapping import (
        AxisArrays,
        AxisArraysView,
        Layers,
        LayersView,
        PairwiseArrays,
        PairwiseArraysView,
    )


type Axes = Collection[Literal["obs", "var"]]

type Idx2D = tuple[str | slice, slice] | tuple[slice, str | slice]
"""Index along the full length of one or both AnnData dimensions, resulting in a 1D or 2D array.

E.g. `a[:, "c"]`, `a["b", :]`, or `a[:, :]`
"""

type Idx2DList = (
    tuple[list[str] | pd.Index[str], slice] | tuple[slice, list[str] | pd.Index[str]]
)
type IdxMultiList = list[int] | pd.Index[int] | tuple[slice, list[int] | pd.Index[int]]

Array: TypeAlias = InMemoryArray | pd.api.extensions.ExtensionArray | XVariable  # noqa: UP040  # private, so autodoc must expand it

type AccJson = Sequence[str | int | AccJson | None]
"""JSON form of an accessor, see :meth:`AdAcc.to_json`."""

type DataFrameLike = pd.DataFrame | Dataset2D
"""A 2D dataframe-like container (pandas- or xarray-backed)."""

NO_IDX = sentinel("NO_IDX")
"""Sentinel object needed for implementing :meth:`anndata.acc.RefAcc.get` when subclassing."""

__all__ = [
    "NO_IDX",
    "A",
    "AccJson",
    "AdAcc",
    "AdRef",
    "DFAcc",
    "GraphAcc",
    "GraphMapAcc",
    "Idx2D",
    "ItemAcc",
    "LayerAcc",
    "LayerMapAcc",
    "MetaAcc",
    "MultiAcc",
    "MultiMapAcc",
    "RefAcc",
    "ScopeAcc",
    "UnsMapAcc",
]


class AdRef[I: Hashable, D: MuData | AnnData]:
    r"""A reference to a 1D or 2D array along one or two dimensions of an AnnData object.

    Examples
    --------
    Use :data:`A` in the same way you would an :class:`~anndata.AnnData` object if you wanted to access a 1D or 2D array along the whole dimension.
    You can then introspect the reference:

    >>> from anndata.acc import A
    >>> ref = A.obs["x"]
    >>> type(ref)
    <class 'anndata.acc.AdRef'>
    >>> ref.dims
    {'obs'}
    >>> ref.idx
    'x'
    >>> ref.acc
    A.obs

    See :mod:`anndata.acc` and :class:`AdAcc` for more examples.

    """

    acc: RefAcc[Self, I, D, object, Array]
    r"""The accessor containing information about this array.

    See :term:`reference accessor`\ s for all possible types this can assume.
    """

    idx: I
    """The index into the complete array referenced by :attr:`acc` making it 1D or 2D.

    E.g. a key of a `obs` or a :type:`Idx2D` into a 2D array.
    """

    __match_args__: ClassVar = ("acc", "idx")

    def __init__(self, acc: RefAcc[Self, I, D, object, Array], idx: I) -> None:
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
                if (other := a.resolve(value, strict=False)) is not None:
                    return self == other
        return NotImplemented

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
    scope: ScopeAcc | None
    """Accessor for the element whose `.AD` cast the mapping is in, see :attr:`anndata.acc.RefAcc.scope`."""

    @abc.abstractmethod
    def __getitem__(self, k: str, /) -> R:
        """Get a reference accessor for mapped array.

        Parameters
        ----------
        k
            The key of the mapped array.
        """


@dataclass(frozen=True)
class RefAcc[R: AdRef, I: Hashable, D: MuData | AnnData, E = object, V = object](
    abc.ABC
):
    r"""Abstract base class for reference accessors.

    See :term:`reference accessor`\ s for all existing subclasses.
    """

    _: KW_ONLY
    ref_class: type[R]
    scope: ScopeAcc | None = None
    """Accessor for the element whose `.AD` cast this accessor is in, e.g. `A.uns['a']` for `A.uns['a'].AD.obsm['x'].DF`.

    If `None`, this accessor refers into the object passed to :meth:`~anndata.acc.RefAcc.get`/:meth:`~anndata.acc.RefAcc.set`/:meth:`~anndata.acc.RefAcc.isin` itself.
    """

    def container(self, data: D, /) -> MuData | AnnData:
        """Get the AnnData object this accessor indexes into.

        That is `data` itself, or the element :attr:`~anndata.acc.RefAcc.scope` refers to in `data`.
        Call this in :meth:`~anndata.acc.RefAcc.get`, :meth:`~anndata.acc.RefAcc.set`, and :meth:`~anndata.acc.RefAcc.isin` so your subclass works when nested.

        Examples
        --------
        >>> from anndata import AnnData
        >>> from anndata.acc import A
        >>> adata = AnnData(shape=(3, 2), uns={"n": AnnData(shape=(5, 4))})
        >>> A.obs.container(adata) is adata
        True
        >>> A.uns["n"].AD.obs.container(adata) is adata.uns["n"]
        True

        Raises
        ------
        TypeError
            If the element :attr:`~anndata.acc.RefAcc.scope` refers to isn’t an AnnData object.
        """
        return _resolve(self.scope, data)

    def process_idx(self, idx: Any, /) -> I:  # noqa: ANN401  # subclasses may override `idx` however they want
        self.dims(idx)
        return idx

    def __getitem__(  # only accessors with array values (`V`) create references
        self: RefAcc[R, I, D, object, Array],
        idx: Any,  # noqa: ANN401  # subclasses may override `idx` however they want
        /,
    ) -> R:
        """Get a reference to an array.

        Parameters
        ----------
        idx
            Index into the referenced element.
            Subclasses define which indices are valid.
        """
        idx = self.process_idx(idx)
        return self.ref_class(self, idx)

    @abc.abstractmethod
    def dims(self, idx: I, /) -> Axes:
        """Get along which dimensions the referenced array is."""

    @abc.abstractmethod
    def __repr__(self, /) -> str:
        """Get a display representation of the accessor.

        This is not a simplified representation supported by :meth:`AdAcc.resolve`, but instead the code used to create the accessor,
        e.g. `A.layers['counts']`.
        """

    @abc.abstractmethod
    def idx_repr(self, idx: I, /) -> str:
        """Get a string representation of the index."""

    @abc.abstractmethod
    def isin(self, data: D, idx: I | None = None, /) -> bool:
        """Check if the referenced array is in the AnnData object.

        Only called if :attr:`~anndata.acc.RefAcc.scope` (or :attr:`~anndata.acc.DFAcc.frame_acc`) is `None`, or in `data` with the type this accessor casts it to (`.AD`/`.DF`), so :meth:`~anndata.acc.RefAcc.container` won’t fail.
        To check an accessor including its parents, use `acc in data`.
        """

    @overload
    def get(self, data: D, /) -> E: ...
    @overload
    def get(self, data: D, idx: I, /) -> V: ...
    @abc.abstractmethod
    def get(self, data: D, idx: I | NO_IDX = NO_IDX, /) -> E | V:
        """Get the indexed array or element from the AnnData object at `idx`.

        When `idx` is omitted (i.e., `idx` is :class:`~anndata.acc.NO_IDX`), return the full array one level up instead.
        This has the same semantics as the `AdRef` path but one level up:
        `adata[A.obs]` returns the full :class:`~pandas.DataFrame`, `adata[A.obsm["pca"]]` the full :class:`numpy.ndarray`, and `adata[A.layers]` the full :class:`Mapping`.
        The first two both have defined `shape`-like properties (or :class:`awkward.Array`), but :attr:`~anndata.AnnData.layers` or similar do not.
        """

    @future_abstractmethod(Deprecation("0.14.0"))
    def set(self, data: D, idx: I | NO_IDX, value: object, /) -> None:
        """Set the indexed array or element in the AnnData object at `idx`.

        When `idx` is :class:`~anndata.acc.NO_IDX`, replace the full element one level up (see :meth:`~anndata.acc.RefAcc.get`),
        otherwise assign into the existing element in-place.
        """
        msg = f"{type(self).__name__} does not support setting values."
        raise NotImplementedError(msg)

    def _maybe_flatten(
        self, idx: I, a: AlignedArray | pd.api.extensions.ExtensionArray
    ) -> Array:
        if not (_is_in_memory(a) or isinstance(a, pd.api.extensions.ExtensionArray)):
            msg = f"Cannot get {self.idx_repr(idx)} from {type(a).__name__}, load it into memory first."
            raise TypeError(msg)
        if len(self.dims(idx)) != 1:
            return a
        if isinstance(a, DaskArray):
            import dask.array as da

            a = da.map_blocks(lambda x: self._maybe_flatten(idx, x), a)
        if isinstance(a, CSMatrix | CSArray | CupySparseMatrix | ArrayView):
            a = a.toarray()
        if isinstance(a, SupportsArrayApiBase):
            return a.__array_namespace__().reshape(a, (a.size,))
        return a.ravel()


@dataclass(frozen=True)
class LayerAcc[R: AdRef[Idx2D, AnnData]](
    RefAcc[R, Idx2D, AnnData, AlignedArray, Array]
):
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
        # TODO: add docstring once sphinx-autodoc-typehints can render `pd.Index[str]` in `Idx2DList`
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def dims(self, idx: Idx2D, /) -> set[Literal["obs", "var"]]:
        if not isinstance(idx, tuple):
            msg = f"Unsupported index {idx!r} is not a tuple. Did you mean `{self}[{idx!r}, :]`?"
            raise TypeError(msg)
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
        root = _root(self.scope)
        return f"{root}.X" if self.k is None else f"{root}.layers[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def container(self, data: AnnData, /) -> AnnData:
        if not isinstance(inner := super().container(data), AnnData):
            msg = f"{self!r} can only refer into AnnData objects, not {type(inner).__name__}."
            raise TypeError(msg)
        return inner

    def isin(self, adata: AnnData, idx: Idx2D | None = None) -> bool:
        adata = self.container(adata)
        if adata.X is None if self.k is None else self.k not in adata.layers:
            return False
        for i, dim in zip(idx or (), ("obs", "var"), strict=False):
            if isinstance(i, str):  # at most one str
                return i in getattr(adata, dim).index
        return True  # idx is None or [:, :]

    @overload
    def get(self, adata: AnnData, /) -> AlignedArray: ...
    @overload
    def get(self, adata: AnnData, idx: Idx2D, /) -> InMemoryArray: ...
    def get(self, adata: AnnData, idx: Idx2D | NO_IDX = NO_IDX, /) -> AlignedArray:
        adata = self.container(adata)
        if idx is NO_IDX:
            return adata.X if self.k is None else adata.layers[self.k]
        # To keep things as lazy as possible, we don't reuse the full-array branch here
        arr = adata[idx].X if self.k is None else adata[idx].layers[self.k]
        return self._maybe_flatten(idx, arr)

    def set(self, adata: AnnData, idx: Idx2D | NO_IDX, value: object, /) -> None:
        inner = self.container(adata)
        if idx is NO_IDX:
            if self.k is None:
                inner.X = value
            else:
                inner.layers[self.k] = value
            return
        full = _check_settable(self.get(adata))
        full[_positional_idx(inner, idx, ("obs", "var"))] = value


@dataclass(frozen=True)
class LayerMapAcc[R: AdRef](MapAcc[LayerAcc]):
    r"""Accessor for layers (`A.`\ :attr:`~AdAcc.layers`)."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[LayerAcc] = LayerAcc
    scope: ScopeAcc | None = None

    def __getitem__(self, k: str | None, /) -> LayerAcc[R]:
        """Get a reference accessor for a layer.

        Parameters
        ----------
        k
            The layer’s key, or `None` for :attr:`~anndata.AnnData.X` (see :attr:`AdAcc.X`).
        """
        if not isinstance(k, str | None):
            msg = f"Unsupported layer {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(k, ref_class=self.ref_class, scope=self.scope)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.layers"

    def get(self, adata: AnnData) -> Layers | LayersView:
        return _resolve(self.scope, adata).layers


@dataclass(frozen=True)
class DFAcc[R: AdRef](RefAcc[R, str | None, MuData | AnnData, DataFrameLike, Array]):
    r"""Reference accessor for columns of a data frame in an AnnData object.

    Created by casting another accessor to a data frame (`.DF`), which becomes its :attr:`frame_acc`.
    :class:`MetaAcc` (`A.obs`/`A.var`) is a subclass.

    Examples
    --------
    >>> from anndata.acc import A, DFAcc
    >>> assert isinstance(A.obsm["df"].DF, DFAcc)
    >>> A.obsm["df"].DF["col"]
    A.obsm['df'].DF['col']
    >>> A.uns["df"].DF.index
    A.uns['df'].DF.index

    """

    _: KW_ONLY
    frame_acc: RefAcc | None = None
    """Accessor for the data frame this `.DF` cast refers into, e.g. `A.obsm['df']` for `A.obsm['df'].DF`."""

    def __post_init__(self) -> None:
        if self.frame_acc is not None and self.scope != self.frame_acc.scope:
            msg = f"A `.DF` cast must be in the same scope as {self.frame_acc!r}."
            raise ValueError(msg)

    @property
    def index(self) -> R:
        """Index :class:`AdRef`, e.g. `A.obs.index` or `A.uns['df'].DF.index`."""
        return self[None]

    def frame(self, data: MuData | AnnData, /) -> DataFrameLike:
        """Get the data frame this accessor refers to.

        :meth:`~anndata.acc.RefAcc.get`, :meth:`~anndata.acc.RefAcc.set`, and :meth:`~anndata.acc.RefAcc.isin` use this,
        so subclasses only need to override this to refer to another frame
        (e.g. :class:`MetaAcc` returns `.obs`/`.var` of :meth:`~anndata.acc.RefAcc.container`).

        Raises
        ------
        TypeError
            If the element :attr:`frame_acc` refers to isn’t a data frame.
        """
        frame = None if self.frame_acc is None else self.frame_acc.get(data)
        if not _is_frame(frame):
            msg = f"{self.frame_acc!r} is {type(frame).__name__}, not a data frame, so `.DF` can’t index into it."
            raise TypeError(msg)
        return frame

    def process_idx(self, k: object, /) -> str | None:
        if k is not None and not isinstance(k, str):
            msg = f"Unsupported column {k!r} for {self!r}"
            raise TypeError(msg)
        return k

    @overload
    def __getitem__[K: str | None](self, k: K, /) -> R: ...
    @overload
    def __getitem__[K: str | None](self, k: list[K] | pd.Index[str], /) -> list[R]: ...
    def __getitem__[K: str | None](
        self, k: K | list[K] | pd.Index[str], /
    ) -> R | list[R]:
        # TODO: add docstring once sphinx-autodoc-typehints can render `pd.Index[str]`
        # _is_t_list doesn’t allow None, so be more liberal when checking for lists
        if isinstance(k, list) or _is_t_list(k, str):
            return [self[i] for i in k]

        return super().__getitem__(k)

    def dims(self, k: str | None, /) -> Axes:
        # only frames in `obsm`/`varm` are aligned to a dimension
        return {self.frame_acc.dim} if isinstance(self.frame_acc, MultiAcc) else set()

    def __repr__(self) -> str:
        return f"{self.frame_acc!r}.DF"

    def idx_repr(self, k: str | None) -> str:
        return ".index" if k is None else f"[{k!r}]"

    def isin(self, data: MuData | AnnData, idx: str | None = None) -> bool:
        return idx is None or idx in self.frame(data)

    @overload
    def get(self, data: MuData | AnnData, /) -> DataFrameLike: ...
    @overload
    def get(
        self, data: MuData | AnnData, k: str | None, /
    ) -> pd.api.extensions.ExtensionArray | XVariable: ...
    def get(
        self, data: MuData | AnnData, k: str | NO_IDX | None = NO_IDX, /
    ) -> DataFrameLike | pd.api.extensions.ExtensionArray | XVariable:
        full = self.frame(data)
        if k is NO_IDX:
            return full
        match full, k:
            case pd.DataFrame() as df, None:
                return df.index.array
            case Dataset2D() as ds, None:
                return ds.true_index.array
            case pd.DataFrame() as df, k:
                return df[k].array
            case Dataset2D() as ds, k:
                return ds[k].variable
            case _:
                msg = f"Unsupported {self!r} container {type(full).__name__}"
                raise TypeError(msg)

    def set(
        self, data: MuData | AnnData, k: str | NO_IDX | None, value: object, /
    ) -> None:
        if k is NO_IDX:
            if self.frame_acc is None:
                msg = f"{self!r} has no `frame_acc` to set the data frame in."
                raise TypeError(msg)
            self.frame_acc.set(data, NO_IDX, value)
            return
        frame: Any = self.frame(data)  # value is validated by the frame itself
        if k is None:
            frame.index = value
        else:
            frame[k] = value


@dataclass(frozen=True)
class MetaAcc[R: AdRef[str | None, MuData | AnnData]](DFAcc[R]):
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

    def frame(self, data: MuData | AnnData, /) -> DataFrameLike:
        return getattr(self.container(data), self.dim)

    def dims(self, k: str | None, /) -> Axes:
        return {self.dim}

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.{self.dim}"

    def set(
        self, data: MuData | AnnData, k: str | NO_IDX | None, value: object, /
    ) -> None:
        # replacing the frame or index goes through AnnData to keep things aligned
        if k is NO_IDX:
            setattr(self.container(data), self.dim, value)
        elif k is None:
            setattr(self.container(data), f"{self.dim}_names", value)
        else:
            super().set(data, k, value)


@dataclass(frozen=True)
class MultiAcc[R: AdRef[int, MuData | AnnData]](
    RefAcc[R, int, MuData | AnnData, AlignedArray, Array]
):
    r"""Reference accessor for arrays from multi-dimensional containers (`A.`\ :attr:`~AdAcc.obsm`/`A.`\ :attr:`~AdAcc.varm`).

    Examples
    --------
    >>> from anndata.acc import A, MultiAcc
    >>> assert isinstance(A.obsm["pca"], MultiAcc)
    >>> A.obsm["pca"][:, 0]
    A.obsm['pca'][:, 0]

    Note that something like this is invalid,
    as it references an array that isn’t aligned to an anndata dimension.

    ..  code:: python

        A.varm["PCs"]["cell-0", :]

    """

    dim: Literal["obs", "var"]
    """Dimension this accessor refers to, e.g. `A.obsm[k].dim == 'obs'`."""

    k: str
    """Key this accessor refers to, e.g. `A.varm['x'].k == 'x'`."""

    @staticmethod
    def process_idx(i: object, /) -> int:
        # bool is an int subclass, but boolean indexing isn’t supported
        if not isinstance(i := _drop_multi_slice(i), int) or isinstance(i, bool):
            msg = f"Unsupported index {i!r}"
            raise TypeError(msg)
        return i

    @overload
    def __getitem__(self, i: int | tuple[slice, int], /) -> R: ...
    @overload
    def __getitem__(self, i: IdxMultiList, /) -> list[R]: ...
    def __getitem__(self, i: int | tuple[slice, int] | IdxMultiList, /) -> R | list[R]:
        # TODO: add docstring once sphinx-autodoc-typehints can render `pd.Index[int]` in `IdxMultiList`
        if _is_t_list(idx := _drop_multi_slice(i), int):
            return [self[j] for j in idx]
        return super().__getitem__(idx)

    def dims(self, i: int, /) -> Axes:
        return {self.dim}

    @property
    def DF(self) -> DFAcc[R]:
        """Treat the array as a data frame, to refer to its columns by name, e.g. `A.obsm['x'].DF['c']`."""
        return DFAcc(ref_class=self.ref_class, scope=self.scope, frame_acc=self)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.{self.dim}m[{self.k!r}]"

    def idx_repr(self, i: int) -> str:
        return f"[:, {i!r}]"

    def isin(self, data: MuData | AnnData, idx: int | None = None) -> bool:
        m: AxisArrays = getattr(self.container(data), f"{self.dim}m")
        if (arr := m.get(self.k)) is None:
            return False
        return idx is None or -arr.shape[1] <= idx < arr.shape[1]

    @overload
    def get(self, data: MuData | AnnData, /) -> AlignedArray: ...
    @overload
    def get(self, data: MuData | AnnData, i: int, /) -> InMemoryArray: ...
    def get(self, data: MuData | AnnData, i: int | NO_IDX = NO_IDX, /) -> AlignedArray:
        full: AlignedArray = getattr(self.container(data), f"{self.dim}m")[self.k]
        if i is NO_IDX:
            return full
        n_cols: int = full.shape[1]
        if not -n_cols <= i < n_cols:
            msg = (
                f"Column index `{i}` is out of range for {self} with {n_cols} columns."
            )
            raise IndexError(msg)
        i += n_cols * (i < 0)
        if isinstance(full, pd.DataFrame):
            return full.iloc[:, i].array
        if isinstance(full, Dataset2D):
            return full[full.columns[i]].variable
        # TODO: remove slicing when dropping scipy <1.14
        return self._maybe_flatten(i, full[:, i : i + 1])

    def set(self, data: MuData | AnnData, i: int | NO_IDX, value: object, /) -> None:
        if i is NO_IDX:
            getattr(self.container(data), f"{self.dim}m")[self.k] = value
            return
        full = _check_settable(self.get(data))
        if not -full.shape[1] <= i < full.shape[1]:
            msg = f"Column index `{i}` is out of range for {self} with {full.shape[1]} columns."
            raise IndexError(msg)
        full[:, i] = value


@dataclass(frozen=True)
class MultiMapAcc[R: AdRef](MapAcc[MultiAcc]):
    r"""Accessor for multi-dimensional array containers (`A.`\ :attr:`~AdAcc.obsm`/`A.`\ :attr:`~AdAcc.varm`)."""

    dim: Literal["obs", "var"]
    """Dimension this accessor refers to, e.g. `A.obsm.dim == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[MultiAcc] = MultiAcc
    scope: ScopeAcc | None = None

    def __getitem__(self, k: str, /) -> MultiAcc[R]:
        """Get a reference accessor for an array, e.g. `A.obsm["pca"]`.

        Parameters
        ----------
        k
            The array’s key.
        """
        if not isinstance(k, str):
            msg = f"Unsupported {self.dim}m key {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(self.dim, k, ref_class=self.ref_class, scope=self.scope)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.{self.dim}m"

    def get(self, adata: AnnData) -> AxisArrays | AxisArraysView:
        return getattr(_resolve(self.scope, adata), f"{self.dim}m")


@dataclass(frozen=True)
class GraphAcc[R: AdRef[Idx2D, MuData | AnnData]](
    RefAcc[R, Idx2D, MuData | AnnData, AlignedArray, Array]
):
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
        # TODO: add docstring once sphinx-autodoc-typehints can render `pd.Index[str]` in `Idx2DList`
        if _is_idx2d_list(idx):
            return [self[i] for i in _expand_idx2d_list(idx)]
        return super().__getitem__(idx)

    def dims(self, idx: Idx2D, /) -> Axes:
        n_slices = sum(isinstance(i, slice) for i in idx)
        return (self.dim,) * n_slices if n_slices > 1 else {self.dim}

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.{self.dim}p[{self.k!r}]"

    def idx_repr(self, idx: Idx2D) -> str:
        return f"[{idx[0]!r}, {idx[1]!r}]"

    def isin(self, data: MuData | AnnData, idx: Idx2D | None = None) -> bool:
        data = self.container(data)
        if self.k not in getattr(data, f"{self.dim}p"):
            return False
        if idx is None:
            return True
        index: pd.Index = getattr(data, self.dim).index
        return all(i in index for i in idx if isinstance(i, str))

    @overload
    def get(self, data: MuData | AnnData, /) -> AlignedArray: ...
    @overload
    def get(self, data: MuData | AnnData, idx: Idx2D, /) -> InMemoryArray: ...
    def get(
        self, data: MuData | AnnData, idx: Idx2D | NO_IDX = NO_IDX, /
    ) -> AlignedArray:
        data = self.container(data)
        full: AlignedArray = getattr(data, f"{self.dim}p")[self.k]
        if idx is NO_IDX:
            return full
        index: pd.Index = getattr(data, self.dim).index
        # a scalar index makes scipy return a 1D COO array, which cannot be flattened
        match idx:
            case str() as i, slice() as j:
                sub = full[np.array([index.get_loc(i)]), j]
            case slice() as i, str() as j:
                sub = full[i, np.array([index.get_loc(j)])]
            case slice() as i, slice() as j:
                sub = full[i, j]
            case _:  # pragma: no cover
                msg = f"Invalid index: {idx}"
                raise TypeError(msg)
        return self._maybe_flatten(idx, sub)

    def set(
        self, data: MuData | AnnData, idx: Idx2D | NO_IDX, value: object, /
    ) -> None:
        inner = self.container(data)
        if idx is NO_IDX:
            getattr(inner, f"{self.dim}p")[self.k] = value
            return
        full = _check_settable(self.get(data))
        full[_positional_idx(inner, idx, (self.dim, self.dim))] = value


@dataclass(frozen=True)
class GraphMapAcc[R: AdRef](MapAcc[GraphAcc]):
    r"""Accessor for graph containers (`A.`\ :attr:`~AdAcc.obsp`/`A.`\ :attr:`~AdAcc.varp`)."""

    dim: Literal["obs", "var"]
    """Axis this accessor refers to, e.g. `A.obsp.dim == 'obs'`."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[GraphAcc] = GraphAcc
    scope: ScopeAcc | None = None

    def __getitem__(self, k: str, /) -> GraphAcc[R]:
        """Get a reference accessor for a graph, e.g. `A.obsp["distances"]`.

        Parameters
        ----------
        k
            The graph’s key.
        """
        if not isinstance(k, str):
            msg = f"Unsupported {self.dim}p key {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls(self.dim, k, ref_class=self.ref_class, scope=self.scope)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.{self.dim}p"

    def get(self, adata: AnnData) -> PairwiseArrays | PairwiseArraysView:
        return getattr(_resolve(self.scope, adata), f"{self.dim}p")


@dataclass(frozen=True)
class ScopeAcc[R: AdRef, D: MuData | AnnData, E = object](RefAcc[R, Never, D, E, E]):
    r"""Abstract base class for accessors to elements that contain an AnnData-like object.

    Accessors for arrays in that object (created via :attr:`AD`)
    have it as their :attr:`~RefAcc.scope`.
    """

    @property
    @abc.abstractmethod
    def AD(self) -> AdAcc[R]:
        """Treat the element as a nested AnnData object, e.g. `A.uns['a'].AD.obs['x']`.

        Its :func:`repr` is used as the prefix of nested accessors’ :func:`repr`\\ s.
        """

    @abc.abstractmethod
    def get(self, data: D, idx: NO_IDX = NO_IDX, /) -> E:
        """Get the referenced element."""


@dataclass(frozen=True)
class ItemAcc[R: AdRef](ScopeAcc[R, MuData | AnnData]):
    r"""Accessor for (possibly nested) items in `.uns` (`A.`\ :attr:`~AdAcc.uns`\ `[k]`).

    Indexing it with a key descends further into the item (treating it as a mapping).
    It can’t create :term:`reference`\ s itself, as items in `.uns` aren’t aligned to any dimension,
    but it can be cast into a type with arrays inside (`.AD`, `.DF`).

    Examples
    --------
    >>> from anndata.acc import A, ItemAcc
    >>> assert isinstance(A.uns["foo"], ItemAcc)
    >>> A.uns["foo"]["bar"]
    A.uns['foo']['bar']
    >>> A.uns["baz"].AD.obs["x"]
    A.uns['baz'].AD.obs['x']
    >>> A.uns["df"].DF["col"]
    A.uns['df'].DF['col']

    """

    path: tuple[str, ...]
    """Keys leading to the item, e.g. `A.uns['a']['b'].path == ('a', 'b')`."""

    def __getitem__(self, k: str, /) -> Self:  # type: ignore[override]  # descends instead of creating refs
        """Get an accessor for an item nested in this one, e.g. `A.uns['a']['b']`.

        Parameters
        ----------
        k
            The key of the nested item.
        """
        if not isinstance(k, str):
            msg = f"Unsupported uns key {k!r}"
            raise TypeError(msg)
        return replace(self, path=(*self.path, k))

    @property
    def AD(self) -> AdAcc[R]:
        """Treat the item as a nested AnnData object, e.g. `A.uns['a'].AD.obs['x']`."""
        return AdAcc(ref_class=self.ref_class, scope=self)

    @property
    def DF(self) -> DFAcc[R]:
        """Treat the item as a data frame, e.g. `A.uns['a'].DF['col']`."""
        return DFAcc(ref_class=self.ref_class, scope=self.scope, frame_acc=self)

    def dims(self, idx: Never, /) -> Axes:
        msg = f"{self!r} doesn’t contain arrays aligned to a dimension. Cast it, e.g. with `.AD` or `.DF`."
        raise TypeError(msg)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.uns" + "".join(f"[{k!r}]" for k in self.path)

    def idx_repr(self, idx: Never) -> str:
        return f"[{idx!r}]"

    def isin(self, data: MuData | AnnData, idx: None = None) -> bool:
        item = self.container(data).uns
        try:
            for k in self.path:
                item = item[k]
        except (KeyError, IndexError, TypeError):  # missing, or not a mapping
            return False
        return True

    def get(self, data: MuData | AnnData, idx: NO_IDX = NO_IDX, /) -> object:
        if idx is not NO_IDX:
            msg = f"{self!r} can only be used to get the whole item."
            raise TypeError(msg)
        item = self.container(data).uns
        for k in self.path:
            item = item[k]
        return item

    def set(self, data: MuData | AnnData, idx: NO_IDX, value: object, /) -> None:
        if idx is not NO_IDX:
            msg = f"{self!r} can only be used to set the whole item."
            raise TypeError(msg)
        container = self.container(data).uns
        for k in self.path[:-1]:
            container = container[k]
        container[self.path[-1]] = value


@dataclass(frozen=True)
class UnsMapAcc[R: AdRef](MapAcc[ItemAcc]):
    r"""Accessor for unstructured metadata (`A.`\ :attr:`~AdAcc.uns`)."""

    _: KW_ONLY
    ref_class: type[R]
    ref_acc_cls: type[ItemAcc] = ItemAcc
    scope: ScopeAcc | None = None

    def __getitem__(self, k: str, /) -> ItemAcc[R]:
        """Get an accessor for an item, e.g. `A.uns["foo"]`.

        Parameters
        ----------
        k
            The item’s key.
        """
        if not isinstance(k, str):
            msg = f"Unsupported uns key {k!r}"
            raise TypeError(msg)
        return self.ref_acc_cls((k,), ref_class=self.ref_class, scope=self.scope)

    def __repr__(self) -> str:
        return f"{_root(self.scope)}.uns"

    def get(self, adata: AnnData) -> MutableMapping[str, Any]:
        return _resolve(self.scope, adata).uns


# `R`’s bound can’t double as its default before PEP 696 (Python 3.13)
_DEFAULT_REF_CLASS: Any = AdRef


@dataclass(frozen=True)
class AdAcc[R: AdRef]:
    r"""Accessor to create :class:`AdRef`\ s (:data:`A`).

    See examples below and in :mod:`anndata.acc`.
    """

    ref_class: type[R] = _DEFAULT_REF_CLASS

    layer_cls: type[LayerAcc] = LayerAcc
    """Class to use for `layers` accessors."""

    meta_cls: type[MetaAcc] = MetaAcc
    """Class to use for `obs`/`var` accessors."""

    multi_cls: type[MultiAcc] = MultiAcc
    """Class to use for `obsm`/`varm` accessors."""

    graph_cls: type[GraphAcc] = GraphAcc
    """Class to use for `obsp`/`varp` accessors."""

    item_cls: type[ItemAcc] = ItemAcc
    """Class to use for `uns` item accessors."""

    scope: ScopeAcc | None = field(default=None, kw_only=True)
    """Accessor for the element this object is a cast of, e.g. `A.uns['a']` for `A.uns['a'].AD`.

    `None` for a top-level accessor like :data:`A`.
    """

    X: LayerAcc[R] = field(init=False)
    """Access data matrix (shortcut for `A.layers[None]`).

    >>> A.X[:, :]
    >>> A.X["cell-1", :]
    >>> A.X[:, "gene-5"]
    """

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

    uns: UnsMapAcc[R] = field(init=False)
    """Access (nested) items in unstructured metadata.

    >>> A.uns["neighbors"]["params"]
    A.uns['neighbors']['params']
    >>> A.uns["nested"].AD.obs["x"]
    A.uns['nested'].AD.obs['x']
    """

    ATTRS: ClassVar = frozenset({
        *("layers", "obs", "var"),
        *("obsm", "varm", "obsp", "varp"),
        "uns",
    })

    def __post_init__(self) -> None:
        r, s = self.ref_class, self.scope
        x = self.layer_cls(None, ref_class=r, scope=s)
        layers = LayerMapAcc(ref_class=r, ref_acc_cls=self.layer_cls, scope=s)
        uns = UnsMapAcc(ref_class=r, ref_acc_cls=self.item_cls, scope=s)
        object.__setattr__(self, "X", x)
        object.__setattr__(self, "layers", layers)
        object.__setattr__(self, "uns", uns)
        for dim in ("obs", "var"):
            meta = self.meta_cls(dim, ref_class=r, scope=s)
            multi = MultiMapAcc(dim, ref_class=r, ref_acc_cls=self.multi_cls, scope=s)
            graphs = GraphMapAcc(dim, ref_class=r, ref_acc_cls=self.graph_cls, scope=s)
            object.__setattr__(self, dim, meta)
            object.__setattr__(self, f"{dim}m", multi)
            object.__setattr__(self, f"{dim}p", graphs)

    def to_json(
        self, ref: R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]
    ) -> list[str | int | AccJson | None]:
        """Serialize an :class:`AdRef` or a whole container accessor to a JSON-compatible list.

        Schema: `acc-schema-v1.json <../acc-schema-v1.json>`_
        (`#/$defs/ref` matches vectors, `#/$defs/acc` whole containers)

        Examples
        --------
        >>> A.to_json(A.obsm["pca"][0])
        ['obsm', 'pca', 0]
        >>> A.to_json(A.obsm["pca"])
        ['obsm', 'pca']
        >>> A.to_json(A.X)
        ['layers', None]

        Items in `.uns` have their path in a list,
        and casts (`.AD`, `.DF`) add a trailing `[cast, payload]` pair:

        >>> A.to_json(A.uns["a"]["b"])
        ['uns', ['a', 'b']]
        >>> A.to_json(A.obsm["df"].DF["c"])
        ['obsm', 'df', ['DF', 'c']]
        >>> A.to_json(A.uns["n"].AD.obs["x"])
        ['uns', ['n'], ['AD', ['obs', 'x']]]
        """
        from ._parse_json import to_json

        return to_json(ref)

    @overload
    def from_json(self, data: AccJson, *, vec: Literal[True]) -> R: ...
    @overload
    def from_json(
        self, data: AccJson, *, vec: Literal[False]
    ) -> LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]: ...
    @overload
    def from_json(
        self, data: AccJson, *, vec: None = None
    ) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]: ...
    def from_json(
        self, data: AccJson, *, vec: bool | None = None
    ) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | DFAcc[R]:
        """Create an :class:`AdRef` or a whole container accessor from a JSON sequence.

        `vec` works like in :meth:`resolve`:
        if `True`, `data` must refer to a vector (:class:`AdRef`),
        if `False`, to a whole container (:class:`LayerAcc`/:class:`MultiAcc`/:class:`GraphAcc`/:class:`ItemAcc`/:class:`DFAcc`),
        and if unset, both are accepted.

        Schema: `acc-schema-v1.json <../acc-schema-v1.json>`_
        (`#/$defs/ref` matches vectors, `#/$defs/acc` whole containers)

        Examples
        --------
        >>> A.from_json(["obsm", "pca", 0])
        A.obsm['pca'][:, 0]
        >>> A.from_json(["obsm", "pca"])
        A.obsm['pca']

        Raises
        ------
        ValueError
            If parsing fails.
        """
        from ._parse_json import parse_json

        try:
            return parse_json(self, data, vec=vec)
        except Exception as e:
            msg = f"Failed to parse {data!r}: {e}"
            raise ValueError(msg) from e

    @overload
    def resolve(
        self, spec: str, *, strict: Literal[True] = True, vec: Literal[True]
    ) -> R: ...
    @overload
    def resolve(
        self, spec: str, *, strict: Literal[False], vec: Literal[True]
    ) -> R | None: ...
    @overload
    def resolve(
        self, spec: str, *, strict: Literal[True] = True, vec: Literal[False]
    ) -> LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R]: ...
    @overload
    def resolve(
        self, spec: str, *, strict: Literal[False], vec: Literal[False]
    ) -> LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | None: ...
    @overload
    def resolve(
        self, spec: str, *, strict: Literal[True] = True, vec: None = None
    ) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R]: ...
    @overload
    def resolve(
        self, spec: str, *, strict: Literal[False], vec: None = None
    ) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | None: ...
    def resolve(
        self, spec: str, *, strict: bool = True, vec: bool | None = None
    ) -> R | LayerAcc[R] | MultiAcc[R] | GraphAcc[R] | ItemAcc[R] | None:
        """Create :class:`AdRef` from a simplified string.

        If `vec` is `True`, `spec` must be an indexed form yielding an :class:`AdRef`, e.g. `"X[:,:]"`, `"obs.a"`, or `"obsm.c.0"` (the current default behavior).
        If `vec` is `False`, `spec` must refer to a whole container instead (`"X"`, `"layers.<k>"`, `"obsm.<k>"`, `"varm.<k>"`, `"obsp.<k>"`, `"varp.<k>"`, or `"uns.<k1>.….<kn>"`), and a :class:`LayerAcc`/:class:`MultiAcc`/:class:`GraphAcc`/:class:`ItemAcc` is returned instead of an :class:`AdRef`.
        If `vec` is unset, both forms are accepted.

        Examples
        --------
        Indexed, yielding an `AdRef`:

        >>> A.resolve("X[:,:]")
        A.X[:, :]
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

        Whole containers, yielding a `LayerAcc`/`MultiAcc`/`GraphAcc`:

        >>> A.resolve("X", vec=False)
        A.X
        >>> A.resolve("layers.y", vec=False)
        A.layers['y']
        >>> A.resolve("obsm.c", vec=False)
        A.obsm['c']
        >>> A.resolve("obsp.g", vec=False)
        A.obsp['g']

        Items in `.uns`, yielding an :class:`ItemAcc`.
        Like for `"obsm.<k>.<i>"`, the path is split at periods.

        >>> A.resolve("uns.neighbors", vec=False)
        A.uns['neighbors']
        >>> A.resolve("uns.a.b.0", vec=False)
        A.uns['a']['b']['0']
        """
        from ._parse_str import parse

        return parse(self, spec, strict=strict, vec=vec)

    def __repr__(self) -> str:
        return "A" if self.scope is None else f"{self.scope!r}.AD"


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


def _drop_multi_slice(i: object, /) -> object:
    """Turn a `[:, i]` index into a plain `i` one."""
    if not isinstance(i, tuple):
        return i
    if len(i) != 2 or i[0] != slice(None):
        msg = f"Unsupported slice {i!r}"
        raise ValueError(msg)
    return i[1]


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


@runtime_checkable
class _Settable(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...
    def __setitem__(self, idx: tuple[int | slice, ...], value: object, /) -> None: ...


def _check_settable(a: AlignedArray, /) -> _Settable:
    if not isinstance(a, _Settable):
        msg = f"Cannot assign into {type(a).__name__}."
        raise TypeError(msg)
    return a


def _root(scope: ScopeAcc | None, /) -> str:
    """Get the repr of the :class:`AdAcc` an accessor in `scope` belongs to."""
    # the scope’s `.AD` decides how it’s written, e.g. `A.uns['a'].AD`, or mudata’s `A.mod['rna']`
    return "A" if scope is None else repr(scope.AD)


def _resolve(scope: ScopeAcc | None, data: MuData | AnnData, /) -> MuData | AnnData:
    """Get the AnnData object `scope` refers to in `data`, validating the `.AD` cast."""
    if scope is None:
        return data
    if not _is_anndata_like(elem := scope.get(data)):
        msg = f"{scope!r} is {type(elem).__name__}, not an AnnData object, so `.AD` can’t index into it."
        raise TypeError(msg)
    return elem


def _cast_in(acc: RefAcc | MapAcc, data: MuData | AnnData, /) -> bool:
    """Check if the element `acc` is cast from is in `data` and has the right type."""
    match acc:
        case DFAcc(frame_acc=RefAcc() as frame):  # `.DF`
            return frame in data and _is_frame(frame.get(data))
        case _ if (scope := acc.scope) is not None:  # `.AD`
            return scope in data and _is_anndata_like(scope.get(data))
    return True


def _is_frame(x: object, /) -> TypeGuard[DataFrameLike]:
    return isinstance(x, pd.DataFrame | Dataset2D)


def _is_anndata_like(x: object, /) -> TypeGuard[MuData | AnnData]:
    if isinstance(x, AnnData):
        return True
    # `MuData` is a placeholder at runtime, since importing mudata here would be circular.
    # A MuData object can only exist if mudata was imported, so look it up there.
    return (mudata := sys.modules.get("mudata")) is not None and isinstance(
        x, mudata.MuData
    )


def _positional_idx(
    data: MuData | AnnData, idx: Idx2D, dims: tuple[str, str], /
) -> tuple[int | slice, ...]:
    """Turn the names in a 2D index into positions along `dims`."""
    return tuple(
        getattr(data, dim).index.get_loc(i) if isinstance(i, str) else i
        for i, dim in zip(idx, dims, strict=True)
    )


if not TYPE_CHECKING:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    R = AdRef[Hashable, MuData | AnnData]
    I = Hashable
    D = MuData | AnnData
