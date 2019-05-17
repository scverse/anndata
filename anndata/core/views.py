from copy import deepcopy
from functools import reduce, singledispatch
from typing import Any, KeysView, Optional, Sequence, Tuple, NamedTuple

import numpy as np
import pandas as pd
from scipy import sparse

from ..logging import anndata_logger as logger


class ViewArgs(NamedTuple):
    parent: "AnnData"
    attrname: str
    keys: Tuple[str, ...] = ()


class _SetItemMixin:
    def __setitem__(self, idx: Any, value: Any):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            adata_view, attr_name, keys = self._view_args
            logger.warning(
                'Trying to set attribute `.{}` of view, making a copy.'.format(attr_name))
            new = adata_view.copy()
            attr = getattr(new, attr_name)
            container = reduce(lambda d, k: d[k], keys, attr)
            container[idx] = value
            adata_view._init_as_actual(new)


class _ViewMixin(_SetItemMixin):
    def __init__(
        self,
        *args,
        view_args: Tuple['AnnData', str, Tuple[str, ...]] = None,
        **kwargs
    ):
        if view_args is not None:
            view_args = ViewArgs(*view_args)
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __deepcopy__(self, memo):
        parent, attrname, keys = self._view_args
        return deepcopy(getattr(parent._adata_ref, attrname))


class ArrayView(_SetItemMixin, np.ndarray):
    def __new__(
        cls,
        input_array: Sequence[Any],
        view_args: Tuple['AnnData', str, Tuple[str, ...]] = None,
    ):
        arr = np.asarray(input_array).view(cls)
        if view_args is not None:
            view_args = ViewArgs(*view_args)
        arr._view_args = view_args
        return arr

    def __array_finalize__(self, obj: Optional[np.ndarray]):
        if obj is not None:
            self._view_args = getattr(obj, '_view_args', None)

    def keys(self) -> KeysView[str]:
        # it's a structured array
        return self.dtype.names

    def copy(self, order: str = 'C') -> np.ndarray:
        # we want a conventional array
        return np.array(self)

    def toarray(self) -> np.ndarray:
        return self.copy()


class SparseCSRView(_ViewMixin, sparse.csr_matrix):
    pass


class SparseCSCView(_ViewMixin, sparse.csc_matrix):
    pass


class DictView(_ViewMixin, dict):
    pass


class DataFrameView(_ViewMixin, pd.DataFrame):
    _metadata = ['_view_args']


@singledispatch
def asview(obj, view_args):
    raise NotImplementedError(
        "No view type has been registered for {}".format(type(obj))
    )

@asview.register(np.ndarray)
def asview_array(array, view_args):
    return ArrayView(array, view_args=view_args)

@asview.register(pd.DataFrame)
def asview_df(df, view_args):
    return DataFrameView(df, view_args=view_args)

@asview.register(sparse.csr_matrix)
def asview_csr(mtx, view_args):
    return SparseCSRView(mtx, view_args=view_args)

@asview.register(sparse.csc_matrix)
def asview_csc(mtx, view_args):
    return SparseCSCView(mtx, view_args=view_args)

@asview.register(dict)
def asview_dict(d, view_args):
    return DictView(d, view_args=view_args)
