from __future__ import annotations

from pathlib import Path
from typing import (
    TYPE_CHECKING,
)

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import MutableMapping

import dask.array as da
import h5py
import xarray as xr
import zarr

from ..._core.anndata import AnnData
from ..._core.sparse_dataset import BaseCompressedSparseDataset, sparse_dataset
from ...compat import DaskArray
from ...utils import convert_to_dict
from .. import read_dispatched
from ._lazy_arrays import LazyCategoricalArray, LazyMaskedArray
from ._xarray import Dataset2D


def to_memory(adata, exclude=[]):
    # nullable and categoricals need special handling because xarray will convert them to numpy arrays first with dtype object
    def get_nullable_and_categorical_cols(ds):
        cols = []
        for c in ds:
            dtype = ds[c].dtype
            if (
                isinstance(dtype, pd.CategoricalDtype)
                or dtype == pd.arrays.BooleanArray
                or dtype == pd.arrays.IntegerArray
            ):
                cols += [c]
        return cols

    def to_df(ds, exclude_vars=[]):
        nullable_and_categorical_df_cols = get_nullable_and_categorical_cols(ds)
        drop_vars = [
            k for k in set(exclude_vars + nullable_and_categorical_df_cols) if k in ds
        ]
        df = ds.drop_vars(drop_vars).to_dataframe()
        for c in nullable_and_categorical_df_cols:
            if c not in exclude_vars:
                df[c] = ds[c].data[()]
        df.index.name = None  # matches old AnnData object
        if len(exclude_vars) == 0:
            df = df[list(ds.keys())]
        return df

    # handling for AxisArrays
    def backed_dict_to_memory(d, prefix):
        res = {}
        for k, v in d.items():
            full_key = prefix + "/" + k
            if any([full_key == exclude_key for exclude_key in exclude]):
                continue
            if isinstance(v, DaskArray):
                res[k] = v.compute()
            elif isinstance(v, BaseCompressedSparseDataset):
                res[k] = v.to_memory()
            elif isinstance(v, Dataset2D):
                res[k] = to_df(v)
            else:
                res[k] = v
        return res

    exclude_obs = [key.replace("obs/", "") for key in exclude if key.startswith("obs/")]
    obs = to_df(adata.obs, exclude_obs)
    exclude_var = [key.replace("var/", "") for key in exclude if key.startswith("var/")]
    var = to_df(adata.var, exclude_var)
    obsm = backed_dict_to_memory(convert_to_dict(adata.obsm), "obsm")
    varm = backed_dict_to_memory(convert_to_dict(adata.varm), "varm")
    varp = backed_dict_to_memory(convert_to_dict(adata.varp), "varp")
    obsp = backed_dict_to_memory(convert_to_dict(adata.obsp), "obsp")
    layers = backed_dict_to_memory(dict(adata.layers), "layers")
    uns = backed_dict_to_memory(convert_to_dict(adata.uns), "uns")
    X = None
    if "X" not in exclude:
        if isinstance(adata.X, BaseCompressedSparseDataset):
            X = adata.X.to_memory()
        elif isinstance(adata.X, DaskArray):
            X = adata.X.compute()
        else:
            X = adata.X
    return AnnData(
        X=X,
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        layers=layers,
        uns=uns,
    )


def read_backed(
    store: str | Path | MutableMapping | zarr.Group | h5py.Dataset,
) -> AnnData:
    """Lazily read in on-disk/in-cloud AnnData stores.  A new, but familiar, AnnData object will be returned.
    No array data should need to be read into memory, with exception of non-obs/var dataframes and Awkward Arrays.

    Args:
        store (Union[str, Path, MutableMapping, zarr.Group, h5py.Dataset]): A store-like object to be read in.  If `zarr`, it is best
        for it to be consolidated.

    Returns:
        AnnData: A lazily read-in AnnData object.
    """
    is_h5 = False
    if isinstance(store, Path) or isinstance(store, str):
        store = str(store)
        if store.endswith("h5ad"):
            is_h5 = True

    has_keys = True  # true if consolidated or h5ad
    if not is_h5:
        try:
            f = zarr.open_consolidated(store, mode="r")
        except KeyError:
            has_keys = False
            f = zarr.open(store, mode="r")
    else:
        f = h5py.File(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            cols = [
                "obs",
                "var",
                "obsm",
                "varm",
                "obsp",
                "varp",
                "layers",
                "X",
                "raw",
                "uns",
            ]
            iter_object = (
                elem.items() if has_keys else [(k, elem[k]) for k in cols if k in elem]
            )
            return AnnData(**{k: read_dispatched(v, callback) for k, v in iter_object})
        elif elem_name.startswith("/raw"):
            return None
        elif iospec.encoding_type in {"dataframe"}:
            iter_object = [(k, elem[k]) for k in elem.attrs["column-order"]] + [
                (elem.attrs["_index"], elem[elem.attrs["_index"]])
            ]
            d = {k: read_dispatched(v, callback) for k, v in iter_object}
            d_with_xr = {}
            index_label = f'{elem_name.replace("/", "")}_names'
            for k in d:
                v = d[k]
                if type(v) == DaskArray and k != elem.attrs["_index"]:
                    d_with_xr[k] = xr.DataArray(
                        v, coords=[d[elem.attrs["_index"]]], dims=[index_label], name=k
                    )
                elif (
                    type(v) == LazyCategoricalArray or type(v) == LazyMaskedArray
                ) and k != elem.attrs["_index"]:
                    d_with_xr[k] = xr.DataArray(
                        xr.core.indexing.LazilyIndexedArray(v),
                        coords=[d[elem.attrs["_index"]]],
                        dims=[index_label],
                        name=k,
                    )
                elif k == elem.attrs["_index"]:
                    d_with_xr[index_label] = xr.DataArray(
                        v, coords=[v], dims=[index_label], name=index_label
                    )
                else:
                    d_with_xr[k] = v
            return Dataset2D(d_with_xr)
        elif iospec.encoding_type == "categorical":
            drop_unused_cats = False  # always don't because the `AnnData` object cannot drop them for us, so getting tests to pass means we need to leave this.
            return LazyCategoricalArray(
                elem["codes"], elem["categories"], elem.attrs, drop_unused_cats
            )
        elif "nullable" in iospec.encoding_type:
            return LazyMaskedArray(
                elem["values"],
                elem["mask"] if "mask" in elem else None,
                iospec.encoding_type,
            )
        elif iospec.encoding_type in {"array", "string-array"}:
            if is_h5:
                if iospec.encoding_type == "string-array":
                    if (
                        "read_dataset" not in dir()
                    ):  # avoid circular dependency, not sure what caused this all of a sudden after merging https://github.com/scverse/anndata/pull/949/commits/dc9f12fcbca977841e967c8414b9f1032e069250
                        from ..._io.h5ad import read_dataset
                    elem = read_dataset(elem)
                if not hasattr(elem, "chunks") or elem.chunks is None:
                    return da.from_array(elem, chunks=(1000,) * len(elem.shape))
                return da.from_array(elem)
            return da.from_zarr(elem)
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            return sparse_dataset(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
