from anndata._core.aligned_mapping import AxisArraysView

import pandas as pd
import numpy as np

from ..._core import AxisArrays


class AxisArraysRemote(AxisArrays):
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                return self._view(self.parent, (idx,))

        return IlocDispatch()

    @property
    def dim_names(self) -> pd.Index:
        return (self.parent.obs_names, self.parent.var_names)[self._axis].compute()


def to_df_1d_axis_arrays(axis_arrays, idx=None):
    """Convert to pandas dataframe."""
    df = pd.DataFrame(index=axis_arrays.dim_names[() if idx is None else idx])
    for key in axis_arrays.keys():
        if "index" not in key:
            df[key] = axis_arrays[key][() if idx is None else idx]
    return df


class AxisArrays1dRemote(AxisArraysRemote):
    def to_df(self) -> pd.DataFrame:
        return to_df_1d_axis_arrays(self)


class AxisArraysRemoteView(AxisArraysView):
    def to_df(self) -> pd.DataFrame:
        return to_df_1d_axis_arrays(self)
    
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)


AxisArrays1dRemote._view_class = AxisArraysRemoteView