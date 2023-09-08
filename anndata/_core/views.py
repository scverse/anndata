from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
from collections.abc import Sequence, KeysView, Callable, Iterable
from functools import reduce, singledispatch, wraps
from typing import Any, Literal
import warnings

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype
from scipy import sparse

import anndata
from anndata._warnings import ImplicitModificationWarning
from .access import ElementRef
from ..compat import (
    ZappyArray,
    AwkArray,
    DaskArray,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
)


@contextmanager
def view_update(adata_view: anndata.AnnData, attr_name: str, keys: tuple[str, ...]):
    """Context manager for updating a view of an AnnData object.

    Contains logic for "actualizing" a view. Yields the object to be modified in-place.

    Parameters
    ----------
    adata_view
        A view of an AnnData
    attr_name
        Name of the attribute being updated
    keys
        Keys to the attribute being updated

    Yields
    ------

    `adata.attr[key1][key2][keyn]...`
    """
    new = adata_view.copy()
    attr = getattr(new, attr_name)
    container = reduce(lambda d, k: d[k], keys, attr)
    yield container
    adata_view._init_as_actual(new)


class _SetItemMixin:
    """\
    Class which (when values are being set) lets their parent AnnData view know,
    so it can make a copy of itself.
    This implements copy-on-modify semantics for views of AnnData objects.
    """

    _view_args: ElementRef | None

    def __setitem__(self, idx: Any, value: Any):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            warnings.warn(
                f"Trying to modify attribute `.{self._view_args.attrname}` of view, "
                "initializing view as actual.",
                ImplicitModificationWarning,
                stacklevel=2,
            )
            with view_update(*self._view_args) as container:
                container[idx] = value


class _ViewMixin(_SetItemMixin):
    def __init__(
        self,
        *args,
        view_args: tuple["anndata.AnnData", str, tuple[str, ...]] = None,
        **kwargs,
    ):
        if view_args is not None:
            view_args = ElementRef(*view_args)
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    # TODO: This makes `deepcopy(obj)` return `obj._view_args.parent._adata_ref`, fix it
    def __deepcopy__(self, memo):
        parent, attrname, keys = self._view_args
        return deepcopy(getattr(parent._adata_ref, attrname))


_UFuncMethod = Literal["__call__", "reduce", "reduceat", "accumulate", "outer", "inner"]


class ArrayView(_SetItemMixin, np.ndarray):
    def __new__(
        cls,
        input_array: Sequence[Any],
        view_args: tuple["anndata.AnnData", str, tuple[str, ...]] = None,
    ):
        arr = np.asanyarray(input_array).view(cls)

        if view_args is not None:
            view_args = ElementRef(*view_args)
        arr._view_args = view_args
        return arr

    def __array_finalize__(self, obj: np.ndarray | None):
        if obj is not None:
            self._view_args = getattr(obj, "_view_args", None)

    def __array_ufunc__(
        self: ArrayView,
        ufunc: Callable[..., Any],
        method: _UFuncMethod,
        *inputs,
        out: tuple[np.ndarray, ...] | None = None,
        **kwargs,
    ) -> np.ndarray:
        """Makes numpy ufuncs convert all instances of views to plain arrays.

        See https://numpy.org/devdocs/user/basics.subclassing.html#array-ufunc-for-ufuncs
        """

        def convert_all(arrs: Iterable[np.ndarray]) -> Iterable[np.ndarray]:
            return (
                arr.view(np.ndarray) if isinstance(arr, ArrayView) else arr
                for arr in arrs
            )

        if out is None:
            outputs = (None,) * ufunc.nout
        else:
            out = outputs = tuple(convert_all(out))

        results = super().__array_ufunc__(
            ufunc, method, *convert_all(inputs), out=out, **kwargs
        )
        if results is NotImplemented:
            return NotImplemented

        if ufunc.nout == 1:
            results = (results,)
        results = tuple(
            (np.asarray(result) if output is None else output)
            for result, output in zip(results, outputs)
        )
        return results[0] if len(results) == 1 else results

    def keys(self) -> KeysView[str]:
        # it’s a structured array
        return self.dtype.names

    def copy(self, order: str = "C") -> np.ndarray:
        # we want a conventional array
        return np.array(self)

    def toarray(self) -> np.ndarray:
        return self.copy()


# Extends DaskArray
# Calls parent __new__ constructor since
# even calling astype on a dask array
# needs a .compute() call to actually happen.
# So no construction by view casting like ArrayView
class DaskArrayView(_SetItemMixin, DaskArray):
    def __new__(
        cls,
        input_array: DaskArray,
        view_args: tuple["anndata.AnnData", str, tuple[str, ...]] = None,
    ):
        arr = super().__new__(
            cls,
            dask=input_array.dask,
            name=input_array.name,
            chunks=input_array.chunks,
            dtype=input_array.dtype,
            meta=input_array._meta,
            shape=input_array.shape,
        )
        if view_args is not None:
            view_args = ElementRef(*view_args)
        arr._view_args = view_args

        return arr

    def __array_finalize__(self, obj: DaskArray | None):
        if obj is not None:
            self._view_args = getattr(obj, "_view_args", None)

    def keys(self) -> KeysView[str]:
        # it’s a structured array
        return self.dtype.names


# Unlike array views, SparseCSRView and SparseCSCView
# do not propagate through subsetting
class SparseCSRView(_ViewMixin, sparse.csr_matrix):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csr_matrix:
        return sparse.csr_matrix(self).copy()


class SparseCSCView(_ViewMixin, sparse.csc_matrix):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csc_matrix:
        return sparse.csc_matrix(self).copy()


class CupySparseCSRView(_ViewMixin, CupyCSRMatrix):
    def copy(self) -> CupyCSRMatrix:
        return CupyCSRMatrix(self).copy()


class CupySparseCSCView(_ViewMixin, CupyCSCMatrix):
    def copy(self) -> CupyCSCMatrix:
        return CupyCSCMatrix(self).copy()


class CupyArrayView(_ViewMixin, CupyArray):
    def __new__(
        cls,
        input_array: Sequence[Any],
        view_args: tuple["anndata.AnnData", str, tuple[str, ...]] = None,
    ):
        import cupy as cp

        arr = cp.asarray(input_array).view(type=cls)

        if view_args is not None:
            view_args = ElementRef(*view_args)
        arr._view_args = view_args
        return arr

    def copy(self) -> CupyArray:
        import cupy as cp

        return cp.array(self).copy()


class DictView(_ViewMixin, dict):
    pass


class DataFrameView(_ViewMixin, pd.DataFrame):
    _metadata = ["_view_args"]

    @wraps(pd.DataFrame.drop)
    def drop(self, *args, inplace: bool = False, **kw):
        if not inplace:
            return self.copy().drop(*args, **kw)
        with view_update(*self._view_args) as df:
            df.drop(*args, inplace=True, **kw)


@singledispatch
def as_view(obj, view_args):
    raise NotImplementedError(f"No view type has been registered for {type(obj)}")


@as_view.register(np.ndarray)
def as_view_array(array, view_args):
    return ArrayView(array, view_args=view_args)


@as_view.register(DaskArray)
def as_view_dask_array(array, view_args):
    return DaskArrayView(array, view_args=view_args)


@as_view.register(pd.DataFrame)
def as_view_df(df, view_args):
    return DataFrameView(df, view_args=view_args)


@as_view.register(sparse.csr_matrix)
def as_view_csr(mtx, view_args):
    return SparseCSRView(mtx, view_args=view_args)


@as_view.register(sparse.csc_matrix)
def as_view_csc(mtx, view_args):
    return SparseCSCView(mtx, view_args=view_args)


@as_view.register(dict)
def as_view_dict(d, view_args):
    return DictView(d, view_args=view_args)


@as_view.register(ZappyArray)
def as_view_zappy(z, view_args):
    # Previous code says ZappyArray works as view,
    # but as far as I can tell they’re immutable.
    return z


@as_view.register(AwkArray)
def as_view_awkarray(array, view_args):
    # We don't need any specific view behavior for awkward arrays. A slice of an awkward array is always a
    # shallow copy of the original. This implies that setting a record field on a slice never modifies the original.
    # Other fields than records are entirely immutable anyway.
    # See also https://github.com/scverse/anndata/issues/1035#issuecomment-1687619270.
    return AwkArray(array)


@as_view.register(CupyArray)
def as_view_cupy(array, view_args):
    return CupyArrayView(array, view_args=view_args)


@as_view.register(CupyCSRMatrix)
def as_view_cupy_csr(mtx, view_args):
    return CupySparseCSRView(mtx, view_args=view_args)


@as_view.register(CupyCSCMatrix)
def as_view_cupy_csc(mtx, view_args):
    return CupySparseCSCView(mtx, view_args=view_args)


def _resolve_idxs(old, new, adata):
    t = tuple(_resolve_idx(old[i], new[i], adata.shape[i]) for i in (0, 1))
    return t


@singledispatch
def _resolve_idx(old, new, l):
    return old[new]


@_resolve_idx.register(np.ndarray)
def _resolve_idx_ndarray(old, new, l):
    if is_bool_dtype(old):
        old = np.where(old)[0]
    return old[new]


@_resolve_idx.register(np.integer)
@_resolve_idx.register(int)
def _resolve_idx_scalar(old, new, l):
    return np.array([old])[new]


@_resolve_idx.register(slice)
def _resolve_idx_slice(old, new, l):
    if isinstance(new, slice):
        return _resolve_idx_slice_slice(old, new, l)
    else:
        return np.arange(*old.indices(l))[new]


def _resolve_idx_slice_slice(old, new, l):
    r = range(*old.indices(l))[new]
    # Convert back to slice
    start, stop, step = r.start, r.stop, r.step
    if len(r) == 0:
        stop = start
    elif stop < 0:
        stop = None
    return slice(start, stop, step)
