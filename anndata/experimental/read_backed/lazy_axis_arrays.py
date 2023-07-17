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


def to_df_1d_axis_arrays(axis_arrays: AxisArrays, exclude=[]) -> pd.DataFrame:
    """Convert a axis array to dataframe in mememort

    Args:
        axis_arrays (AxisArrays): AxisArrays to be converted
        exclude (list, optional): Keys to exclude from being loaded into the DataFrame/Memory. Defaults to [].

    Returns:
        pd.DataFrame: Potential subset of `axis_arrays` in memory
    """
    df = pd.DataFrame(index=axis_arrays.dim_names[...])
    for key in axis_arrays.keys():
        full_key = axis_arrays.attrname + "/" + key
        if "index" not in key and all(
            [full_key != exclude_key for exclude_key in exclude]
        ):
            df[key] = axis_arrays[key].data # all xarray DataArrays?
    return df


class AxisArraysRemote1dMixin:
    def to_df(self, exclude=[]) -> pd.DataFrame:
        """Convert to a DataFrame

        Args:
             exclude (list, optional): Keys to exclude from being loaded into the DataFrame/Memory. Defaults to [].

        Returns:
            pd.DataFrame: Potential subset of `axis_arrays` in memory
        """
        return to_df_1d_axis_arrays(self, exclude)

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                if type(idx) == list:
                    return self._view(self.parent, np.array(idx))
                return self._view(self.parent, idx)

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

    @property
    def attrname(self) -> str:
        return self.dim


class AxisArrays1dRemote(AxisArraysRemote1dMixin, AxisArraysRemote):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AxisArrays1dRemoteView(AxisArraysRemote1dMixin, AxisArraysView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


AxisArrays1dRemote._view_class = AxisArrays1dRemoteView
