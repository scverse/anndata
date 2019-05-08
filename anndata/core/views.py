from copy import deepcopy
from typing import Any, KeysView, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from ..logging import anndata_logger as logger


class _SetItemMixin:
    def __setitem__(self, idx: Any, value: Any):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            logger.warning(
                'Trying to set attribute `.{}` of view, making a copy.'.format(attr_name))
            new = adata_view.copy()
            getattr(new, attr_name)[idx] = value
            adata_view._init_as_actual(new)


class _ViewMixin(_SetItemMixin):
    def __init__(self, *args, view_args: Tuple['AnnData', str] = None, **kwargs):
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __deepcopy__(self, memo):
        parent, k = self._view_args
        return deepcopy(getattr(parent._adata_ref, k))


class ArrayView(_SetItemMixin, np.ndarray):
    def __new__(
        cls,
        input_array: Sequence[Any],
        view_args: Tuple['AnnData', str] = None,
    ):
        arr = np.asarray(input_array).view(cls)
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
