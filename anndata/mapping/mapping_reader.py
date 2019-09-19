from functools import wraps
from pathlib import Path
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd
from scipy import sparse

from .. import mapping as adm
from ..compat import _from_fixed_length_strings, _clean_uns


class AnnDataReadError(OSError):
    """Error caused while trying to read in AnnData."""

    pass


def report_key_on_error(func):
    """
    A decorator for zarr or hdf5 element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_key_on_error
    ... def read_arr(group):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr")
    >>> z["X"] = [1, 2, 3]
    >>> read_arr(z["X"])
    """

    @wraps(func)
    def func_wrapper(elem, *args, **kwargs):
        try:
            return func(elem, *args, **kwargs)
        except Exception as e:
            if isinstance(e, AnnDataReadError):
                raise e
            else:
                raise AnnDataReadError(str(e))

    return func_wrapper


class MappingReader:
    def __init__(self, dataset_class):
        self.dataset_class = dataset_class

    def create_file(self, filename, mode):
        pass

    def close_file(self, f):
        pass

    def read_attribute(self, value):
        if isinstance(value, adm.SparseDataset):
            return self.read_sparse_dataset(value)
        elif isinstance(value, self.dataset_class):
            return self.read_dataset(value)
        elif isinstance(value, adm.Group):
            return self.read_group(value)
        elif value is None:
            return None

        raise NotImplementedError(str(type(value)) + ' ' + str(value))

    @report_key_on_error
    def read_dataframe_legacy(self, dataset) -> pd.DataFrame:
        """
        Read pre-anndata 0.7 dataframes.
        """
        df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
        df.set_index(df.columns[0], inplace=True)
        return df

    @report_key_on_error
    def read_dataframe(self, group) -> pd.DataFrame:
        if not isinstance(group, adm.Group):
            return self.read_dataframe_legacy(group)
        df = pd.DataFrame({k: self.read_series(group[k]) for k in group.keys()})
        df.set_index(group.attrs["_index"], inplace=True)
        if group.attrs["_index"] == "_index":
            df.index.name = None
        if "column-order" in group.attrs:  # TODO: Should this be optional?
            assert set(group.attrs["column-order"]) == set(df.columns)
            df = df[group.attrs["column-order"]]
        return df

    @report_key_on_error
    def read_series(self, dataset) -> Union[np.ndarray, pd.Categorical]:
        if "categories" in dataset.attrs:
            return pd.Categorical.from_codes(
                dataset[...], dataset.attrs["categories"], ordered=False
            )
        else:
            return dataset[...]

    @report_key_on_error
    def read_group(self, group: adm.Group) -> Union[dict, pd.DataFrame]:
        if group.attrs.get("encoding-type", "") == "dataframe":
            return self.read_dataframe(group)
        d = dict()
        for sub_key, sub_value in group.items():
            d[sub_key] = self.read_attribute(sub_value)
        return d

    @report_key_on_error
    def read_dataset(self, dataset: 'Dataset'):

        value = dataset[()]
        if not hasattr(value, "dtype"):
            return value
        elif isinstance(value.dtype, str):
            pass
        elif issubclass(value.dtype.type, np.string_):
            value = value.astype(str)
            if (
                len(value) == 1
            ):  # Backwards compat, old datasets have strings written as one element 1d arrays
                return value[0]
        elif len(value.dtype.descr) > 1:  # Compound dtype
            value = _from_fixed_length_strings(
                value
            )  # For backwards compat, now strings are written as variable length
        if value.shape == ():
            value = value[()]
        return value

    @report_key_on_error
    def read_sparse_dataset(self, value) -> sparse.spmatrix:
        return value.value

    def read(
        self,
        filename: Union[str, Path],
        backed: Optional[Union[str, bool]] = None,
    ) -> 'AnnData':
        """Read ``.h5ad`` and ``.zarr`` formatted hdf5 file.

        Parameters
        ----------
        filename
            File name of data file.
        backed : {``None``, ``'r'``, ``'r+'``}
            If ``'r'``, load :class:`~anndata.AnnData` in ``backed`` mode instead
            of fully loading it into memory (`memory` mode). If you want to modify
            backed attributes of the AnnData object, you need to choose ``'r+'``.
        chunk_size
            Used only when loading sparse dataset that is stored as dense.
            Loading iterates through chunks of the dataset of this row size
            until it reads the whole dataset.
            Higher size means higher memory consumption and higher loading speed.
        """
        mode = None
        d = {}
        if backed not in {None, False}:
            mode = backed
            if mode is True:
                mode = "r+"
            assert mode in {"r", "r+"}

        raw = {}

        f = self.create_file(filename, mode if mode is not None else 'r')
        import anndata.zarr

        if mode is not None and not isinstance(f, anndata.zarr.File):
            # FIXME passing the filename only works for h5 files
            d = {"filename": filename, "filemode": mode}
        attributes = ["obsm", "varm", "obsp", "varp", "uns", "layers"]
        if mode is None:
            attributes.append('X')

        df_attributes = ["obs", "var"]

        d.update({k: self.read_attribute(f[k]) for k in attributes if k in f})
        for k in df_attributes:
            if k in f:  # Backwards compat
                d[k] = self.read_dataframe(f[k])

        if "raw" in f:
            if "raw/var" in f:
                raw["var"] = self.read_attribute(f["raw/var"])
            if "raw/varm" in f:
                raw["varm"] = self.read_attribute(f["raw/varm"])
            if mode is None and "raw/X" in f:
                raw["X"] = self.read_attribute(f["raw/X"])
        else:  # Legacy case
            if "raw.var" in f:
                raw["var"] = self.read_dataframe(
                    f["raw.var"]
                )  # Backwards compat
            if "raw.varm" in f:
                raw["varm"] = self.read_attribute(f["raw.varm"])
            if mode is None and "raw.X" in f:
                raw["X"] = self.read_attribute(f["raw.X"])
        if mode is not None and isinstance(f, anndata.zarr.File):
            # FIXME, backing manager only supports h5 files
            d['X'] = f['X']
        if len(raw) > 0:
            d["raw"] = raw

        if "X" in d:
            d["dtype"] = d["X"].dtype

        _clean_uns(d)
        if mode is not None and not isinstance(f, anndata.zarr.File):
            # FIXME, file will be opened again by backing manager
            self.close_file(f)
        from anndata import AnnData

        return AnnData(**d)
