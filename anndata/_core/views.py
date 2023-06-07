from contextlib import contextmanager
from copy import deepcopy
from enum import Enum
from functools import reduce, singledispatch, wraps
from typing import Any, KeysView, Optional, Sequence, Tuple
import warnings

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype
from scipy import sparse

import anndata
from anndata._warnings import ImplicitModificationWarning
from .access import ElementRef
from ..compat import ZappyArray, AwkArray, DaskArray


class _SetItemMixin:
    """\
    Class which (when values are being set) lets their parent AnnData view know,
    so it can make a copy of itself.
    This implements copy-on-modify semantics for views of AnnData objects.
    """

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
            with self._update() as container:
                container[idx] = value

    @contextmanager
    def _update(self):
        adata_view, attr_name, keys = self._view_args
        new = adata_view.copy()
        attr = getattr(new, attr_name)
        container = reduce(lambda d, k: d[k], keys, attr)
        yield container
        adata_view._init_as_actual(new)


class _ViewMixin(_SetItemMixin):
    def __init__(
        self,
        *args,
        view_args: Tuple["anndata.AnnData", str, Tuple[str, ...]] = None,
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


class ArrayView(_SetItemMixin, np.ndarray):
    def __new__(
        cls,
        input_array: Sequence[Any],
        view_args: Tuple["anndata.AnnData", str, Tuple[str, ...]] = None,
    ):
        arr = np.asanyarray(input_array).view(cls)

        if view_args is not None:
            view_args = ElementRef(*view_args)
        arr._view_args = view_args
        return arr

    def __array_finalize__(self, obj: Optional[np.ndarray]):
        if obj is not None:
            self._view_args = getattr(obj, "_view_args", None)

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
        view_args: Tuple["anndata.AnnData", str, Tuple[str, ...]] = None,
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

    def __array_finalize__(self, obj: Optional[DaskArray]):
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


class DictView(_ViewMixin, dict):
    pass


class DataFrameView(_ViewMixin, pd.DataFrame):
    _metadata = ["_view_args"]

    @wraps(pd.DataFrame.drop)
    def drop(self, *args, inplace: bool = False, **kw):
        if not inplace:
            return self.copy().drop(*args, **kw)
        with self._update() as df:
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


try:
    from ..compat import awkward as ak
    import weakref

    # Registry to store weak references from AwkwardArrayViews to their parent AnnData container
    _registry = weakref.WeakValueDictionary()
    _PARAM_NAME = "_view_args"

    class AwkwardArrayView(_ViewMixin, AwkArray):
        @property
        def _view_args(self):
            """Override _view_args to retrieve the values from awkward arrays parameters.

            Awkward arrays cannot be subclassed like other python objects. Instead subclasses need
            to be attached as "behavior". These "behaviors" cannot take any additional parameters (as we do
            for other data types to store `_view_args`). Therefore, we need to store `_view_args` using awkward's
            parameter mechanism. These parameters need to be json-serializable, which is why we can't store
            ElementRef directly, but need to replace the reference to the parent AnnDataView container with a weak
            reference.
            """
            parent_key, attrname, keys = self.layout.parameter(_PARAM_NAME)
            parent = _registry[parent_key]
            return ElementRef(parent, attrname, keys)

        def __copy__(self) -> AwkArray:
            """
            Turn the AwkwardArrayView into an actual AwkwardArray with no special behavior.

            Need to override __copy__ instead of `.copy()` as awkward arrays don't implement `.copy()`
            and are copied using python's standard copy mechanism in `aligned_mapping.py`.
            """
            array = self
            # makes a shallow copy and removes the reference to the original AnnData object
            array = ak.with_parameter(self, _PARAM_NAME, None)
            array = ak.with_parameter(array, "__array__", None)
            return array

    @as_view.register(AwkArray)
    def as_view_awkarray(array, view_args):
        parent, attrname, keys = view_args
        parent_key = f"target-{id(parent)}"
        _registry[parent_key] = parent
        # TODO: See https://github.com/scverse/anndata/pull/647#discussion_r963494798_ for more details and
        # possible strategies to stack behaviors.
        # A better solution might be based on xarray-style "attrs", once this is implemented
        # https://github.com/scikit-hep/awkward/issues/1391#issuecomment-1412297114
        if type(array).__name__ != "Array":
            raise NotImplementedError(
                "Cannot create a view of an awkward array with __array__ parameter. "
                "Please open an issue in the AnnData repo and describe your use-case."
            )
        array = ak.with_parameter(array, _PARAM_NAME, (parent_key, attrname, keys))
        array = ak.with_parameter(array, "__array__", "AwkwardArrayView")
        return array

    ak.behavior["AwkwardArrayView"] = AwkwardArrayView

except ImportError:

    class AwkwardArrayView:
        pass


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
