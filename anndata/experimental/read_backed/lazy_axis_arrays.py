from typing import Mapping, Union, List
from anndata._core import anndata, raw
from anndata._core.aligned_mapping import AxisArraysBase, AxisArraysView

import pandas as pd
import numpy as np

from ..._core import AxisArrays


class AxisArraysRemote(AxisArrays):
    @property
    def dim_names(self) -> pd.Index:
        return (self.parent.obs_names, self.parent.var_names)[self._axis].compute()

    @property
    def columns(self) -> List:
        return list(self.keys())


def to_df_1d_axis_arrays(axis_arrays):
    """Convert to pandas dataframe."""
    df = pd.DataFrame(index=axis_arrays.dim_names[()])
    for key in axis_arrays.keys():
        if "index" not in key:
            df[key] = axis_arrays[key][()]
    return df


class AxisArraysRemote1dMixin:
    def to_df(self) -> pd.DataFrame:
        return to_df_1d_axis_arrays(self)

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                return self._view(self.parent, (idx,))

        return IlocDispatch()

    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)

    def _repr_html_(self):
        return self.__repr__()

    def _repr_latex_(self):
        return self.__repr__()


class AxisArrays1dRemote(AxisArraysRemote1dMixin, AxisArraysRemote):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AxisArrays1dRemoteView(AxisArraysRemote1dMixin, AxisArraysView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


AxisArrays1dRemote._view_class = AxisArrays1dRemoteView
