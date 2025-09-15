from __future__ import annotations

import warnings
from contextlib import contextmanager
from copy import deepcopy
from functools import reduce, singledispatch, wraps
from typing import TYPE_CHECKING, Literal

import jax.numpy as jnp
import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype
from scipy import sparse

from anndata._warnings import ImplicitModificationWarning

from .._settings import settings
from ..compat import (
    AwkArray,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    DaskArray,
    ZappyArray,
    has_xp,
)
from .access import ElementRef
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, KeysView, Sequence
    from typing import Any, ClassVar

    from numpy.typing import NDArray

    from anndata import AnnData

    from ..compat import Index1DNorm


# def _to_numpy_if_immutable(x):
#     print("Initial x:", type(x))
#     if isinstance(x, np.ndarray | pd.DataFrame | pd.Series | DaskArray):
#         print("x is already a supported mutable type:", type(x))
#         return x
#     try:
#         # checking for mutability
#         print("Trying np.asarray(x)...")
#         x_array = np.asarray(x)
#         print("np.asarray(x) succeeded:", type(x_array))
#         if x_array.size > 0:
#             try:
#                 orig = x_array[0]
#                 print("Checking mutability: trying in-place write")
#                 x_array[0] = orig  # test no-op write
#                 print("In-place mutation succeeded, x is mutable:", type(x))
#                 return x  # mutation worked, so keep original
#             except (ValueError, TypeError) as e:
#                 print("In-place mutation failed:", type(x), "|", repr(e))
#                 # pass

#     except ValueError:
#         print("Trying np.from_dlpack(x)...")
#         result = np.from_dlpack(x)
#         print("np.from_dlpack(x) succeeded:", type(result))
#         # pass
#     # if it is not mutable, we convert to numpy
#     try:
#         # trying convert via from_dlpack first
#         return np.from_dlpack(x)
#     except TypeError:
#         try:
#             # fallback to asarray if from_dlpack not possible
#             print("Trying fallback np.asarray(x)...")
#             result = np.asarray(x)
#             print("Fallback np.asarray(x) succeeded:", type(result))
#             return result
#         except ValueError:
#             # Not an array-API object (or lib not available) = return unchanged
#             print("Final fallback np.asarray(x) failed:", type(x), "|", repr(e))
#             print("Returning x as-is:", type(x))
#             return x


@contextmanager
def view_update(adata_view: AnnData, attr_name: str, keys: tuple[str, ...]):
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
    container = reduce(lambda d, k: d[k], keys, attr)  # attr[k[0]][k[1]][k[2]]...
    # Yield it (not yet converted)
    yield container

    # # After yield, check if immutable and convert to mutable before reinserting
    # parent[key] = _to_numpy_if_immutable(container)
    adata_view._init_as_actual(new)


def _replace_field(container, idx, value):
    import awkward as ak
    import numpy as np

    # JAX-style immutable array
    if hasattr(container, "at"):
        return container.at[idx].set(value)

    # Awkward Array: replace a field or index
    if isinstance(container, ak.Array):
        # Awkward field key (e.g. container["c"] = value)
        if isinstance(idx, str):
            return ak.with_field(container, value, idx)
        # Positional index (e.g. container[0] = value)
        else:
            container_np = np.asarray(container)
            container_np[idx] = value
            return ak.Array(container_np)

    # NumPy or anything else: try normal in-place replacement
    try:
        container_copy = container.copy()
        container_copy[idx] = value
        return container_copy
    except TypeError as err:
        msg = f"Unsupported container type: {type(container)}"
        raise TypeError(msg) from err


class _SetItemMixin:
    """\
    Class which (when values are being set) lets their parent AnnData view know,
    so it can make a copy of itself.
    This implements copy-on-modify semantics for views of AnnData objects.
    """

    _view_args: ElementRef | None

    # def __setitem__(self, idx: Any, value: Any):
    #     # from anndata._core.merge import _to_numpy_if_immutable

    #     if self._view_args is None:
    #         super().__setitem__(idx, value)
    #     else:
    #         warnings.warn(
    #             f"Trying to modify attribute `.{self._view_args.attrname}` of view, "
    #             "initializing view as actual.",
    #             ImplicitModificationWarning,
    #             stacklevel=2,
    #         )
    #         with view_update(*self._view_args) as container:
    #             arr = _to_numpy_if_immutable(container)
    #             arr[idx] = value
    #             # manually assign back into the parent dict
    #             parent = reduce(lambda d, k: d[k], self._view_args[:-1])
    #             parent[self._view_args[-1]] = arr
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
            with view_update(*self._view_args) as (container):
                container[idx] = value
            # potential conversion to numpy
            # with view_update(*self._view_args) as (parent, key, container):
            #     try:
            #         container[idx] = value  # works for NumPy / mutable arrays
            #     except (TypeError, ValueError):
            #         new_container = _replace_field(container, idx, value)
            #         parent[key] = new_container


class _ViewMixin(_SetItemMixin):
    def __init__(
        self,
        *args,
        view_args: tuple[AnnData, str, tuple[str, ...]] | None = None,
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
        view_args: tuple[AnnData, str, tuple[str, ...]] | None = None,
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
            for result, output in zip(results, outputs, strict=True)
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
        view_args: tuple[AnnData, str, tuple[str, ...]] | None = None,
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


# Unlike array views, SparseCSRMatrixView and SparseCSCMatrixView
# do not propagate through subsetting
class SparseCSRMatrixView(_ViewMixin, sparse.csr_matrix):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csr_matrix:
        return sparse.csr_matrix(self).copy()


class SparseCSCMatrixView(_ViewMixin, sparse.csc_matrix):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csc_matrix:
        return sparse.csc_matrix(self).copy()


class SparseCSRArrayView(_ViewMixin, sparse.csr_array):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csr_array:
        return sparse.csr_array(self).copy()


class SparseCSCArrayView(_ViewMixin, sparse.csc_array):
    # https://github.com/scverse/anndata/issues/656
    def copy(self) -> sparse.csc_array:
        return sparse.csc_array(self).copy()


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
        view_args: tuple[AnnData, str, tuple[str, ...]] | None = None,
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
    _metadata: ClassVar = ["_view_args"]

    @wraps(pd.DataFrame.drop)
    def drop(self, *args, inplace: bool = False, **kw):
        if not inplace:
            return self.copy().drop(*args, **kw)
        with view_update(*self._view_args) as df:
            df.drop(*args, inplace=True, **kw)

    def __setattr__(self, key: str, value: Any):
        if key == "index":
            warnings.warn(
                f"Trying to modify {key} of attribute `.{self._view_args.attrname}` of view, "
                "initializing view as actual.",
                ImplicitModificationWarning,
                stacklevel=2,
            )
            with view_update(*self._view_args) as container:
                setattr(container, key, value)
        else:
            super().__setattr__(key, value)


@singledispatch
def as_view(obj, view_args):
    if has_xp(obj):
        # TODO: Determine if we need some sort of specific view object for array-api
        return obj
    msg = f"No view type has been registered for {type(obj)}"
    raise NotImplementedError(msg)


@as_view.register(np.ndarray)
def as_view_array(array, view_args):
    return ArrayView(array, view_args=view_args)


@as_view.register(DaskArray)
def as_view_dask_array(array, view_args):
    return DaskArrayView(array, view_args=view_args)


@as_view.register(pd.DataFrame)
def as_view_df(df, view_args):
    if settings.remove_unused_categories:
        for col in df.columns:
            if isinstance(df[col].dtype, pd.CategoricalDtype):
                with pd.option_context("mode.chained_assignment", None):
                    df[col] = df[col].cat.remove_unused_categories()
    return DataFrameView(df, view_args=view_args)


@as_view.register(sparse.csr_matrix)
def as_view_csr_matrix(mtx, view_args):
    return SparseCSRMatrixView(mtx, view_args=view_args)


@as_view.register(sparse.csc_matrix)
def as_view_csc_matrix(mtx, view_args):
    return SparseCSCMatrixView(mtx, view_args=view_args)


@as_view.register(sparse.csr_array)
def as_view_csr_array(mtx, view_args):
    return SparseCSRArrayView(mtx, view_args=view_args)


@as_view.register(sparse.csc_array)
def as_view_csc_array(mtx, view_args):
    return SparseCSCArrayView(mtx, view_args=view_args)


@as_view.register(dict)
def as_view_dict(d, view_args):
    return DictView(d, view_args=view_args)


@as_view.register(ZappyArray)
def as_view_zappy(z, view_args):
    # Previous code says ZappyArray works as view,
    # but as far as I can tell they’re immutable.
    return z


@as_view.register(CupyArray)
def as_view_cupy(array, view_args):
    return CupyArrayView(array, view_args=view_args)


@as_view.register(CupyCSRMatrix)
def as_view_cupy_csr(mtx, view_args):
    return CupySparseCSRView(mtx, view_args=view_args)


@as_view.register(CupyCSCMatrix)
def as_view_cupy_csc(mtx, view_args):
    return CupySparseCSCView(mtx, view_args=view_args)


@as_view.register(Dataset2D)
def _(a: Dataset2D, view_args):
    return a


try:
    import weakref

    from ..compat import awkward as ak

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
            array = ak.with_parameter(array, "__list__", None)
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
            msg = (
                "Cannot create a view of an awkward array with __array__ parameter. "
                "Please open an issue in the AnnData repo and describe your use-case."
            )
            raise NotImplementedError(msg)
        array = ak.with_parameter(array, _PARAM_NAME, (parent_key, attrname, keys))
        array = ak.with_parameter(array, "__list__", "AwkwardArrayView")
        return array

    ak.behavior["AwkwardArrayView"] = AwkwardArrayView

except ImportError:

    class AwkwardArrayView:
        pass


def _resolve_idxs(
    old: tuple[Index1DNorm, Index1DNorm],
    new: tuple[Index1DNorm, Index1DNorm],
    adata: AnnData,
) -> tuple[Index1DNorm, Index1DNorm]:
    o, v = (_resolve_idx(old[i], new[i], adata.shape[i]) for i in (0, 1))
    return o, v


@singledispatch
def _resolve_idx(old: Index1DNorm, new: Index1DNorm, l: Literal[0, 1]) -> Index1DNorm:
    raise NotImplementedError


@_resolve_idx.register(np.ndarray)
def _resolve_idx_ndarray(
    old: NDArray[np.bool_] | NDArray[np.integer], new: Index1DNorm, l: Literal[0, 1]
) -> NDArray[np.bool_] | NDArray[np.integer]:
    if is_bool_dtype(old) and is_bool_dtype(new):
        mask_new = np.zeros_like(old)
        mask_new[np.flatnonzero(old)[new]] = True
        return mask_new
    if is_bool_dtype(old):
        old = np.where(old)[0]
    return old[new]


@_resolve_idx.register(slice)
def _resolve_idx_slice(
    old: slice, new: Index1DNorm, l: Literal[0, 1]
) -> slice | NDArray[np.integer]:
    if isinstance(new, slice):
        return _resolve_idx_slice_slice(old, new, l)
    else:
        return np.arange(*old.indices(l))[new]


@_resolve_idx.register(jnp.ndarray)
def _resolve_idx_jnp(
    old: jnp.ndarray, new: Index1DNorm, l: Literal[0, 1]
) -> jnp.ndarray:
    # Boolean mask + index
    if old.dtype == jnp.bool_:
        old = jnp.where(old)[0]

    if isinstance(new, slice):
        new = jnp.arange(*new.indices(old.shape[0]))

    return old[new]


def _resolve_idx_slice_slice(old: slice, new: slice, l: Literal[0, 1]) -> slice:
    r = range(*old.indices(l))[new]
    # Convert back to slice
    start, stop, step = r.start, r.stop, r.step
    if len(r) == 0:
        stop = start
    elif stop < 0:
        stop = None
    return slice(start, stop, step)
