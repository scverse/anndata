"""
This module implements on disk dataframes.
"""

import h5py
import pandas as pd


class DataFrame:
    def __init__(self, group: h5py.Group):
        assert (
            group.attrs["encoding-type"] == "dataframe"
        ), "HDF5 group at path '{group.name}' is not encoded as a dataframe"

        self._group = group
        self._attrs = self._group.attrs

        self._index = self._group[self._attrs["_index"]].asstr()
        self._columns = self._attrs["column-order"]

        for column in self.columns:
            # read_elem_partial(group)  # , items=obs, indices=(obs_idx, slice(None)))
            setattr(self, column, self._group[column])

    @property
    def index(self):
        return pd.Index(self._index[:])

    @property
    def columns(self):
        return pd.Index(self._columns)

    def __getitem__(self, index):
        if isinstance(index, str) and index in self.columns:
            return getattr(self, index)

        elif isinstance(index, slice):
            raise NotImplementedError("Slicing is not yet supported.")

        else:
            raise TypeError(f"Invalid index '{index}' of type {type(index)}")
