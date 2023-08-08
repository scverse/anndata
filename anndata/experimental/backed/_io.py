from pathlib import Path
from typing import (
    MutableMapping,
    Union,
)

from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.h5ad import read_dataset
from anndata.compat import DaskArray
import h5py
import zarr
import xarray as xr
import dask.array as da

from .. import read_dispatched
from ._lazy_arrays import LazyCategoricalArray, LazyMaskedArray
from ._xarray import Dataset2D
from ._anndata_backed import AnnDataBacked


def read_backed(
    store: Union[str, Path, MutableMapping, zarr.Group, h5py.Dataset]
) -> AnnDataBacked:
    """Lazily read in on-disk/in-cloud AnnData stores.  A new, but familiar, AnnData object will be returned.
    No array data should need to be read into memory, with exception of non-obs/var dataframes and Awkward Arrays.

    Args:
        store (Union[str, Path, MutableMapping, zarr.Group, h5py.Dataset]): A store-like object to be read in.  If `zarr`, it is best
        for it to be consolidated.

    Returns:
        AnnDataBacked: A lazily read-in AnnData object.
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
            cols = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "raw"]
            iter_object = (
                elem.items() if has_keys else [(k, elem[k]) for k in cols if k in elem]
            )
            return AnnDataBacked(
                **{k: read_dispatched(v, callback) for k, v in iter_object}, file=elem
            )
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
            drop_unused_cats = not (
                elem_name.startswith("/obsm") or elem_name.startswith("/varm")
            )
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
