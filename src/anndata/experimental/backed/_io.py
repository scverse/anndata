from __future__ import annotations

import warnings
from pathlib import Path
from typing import (
    TYPE_CHECKING,
)

import h5py

from anndata._io.specs.registry import read_elem_as_dask

from ..._core.anndata import AnnData
from ...compat import DaskArray
from .. import read_dispatched
from ._compat import xr
from ._lazy_arrays import CategoricalArray, MaskedArray
from ._xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from ...compat import ZarrGroup


def read_backed(
    store: str | Path | MutableMapping | ZarrGroup | h5py.Dataset,
) -> AnnData:
    """Lazily read in on-disk/in-cloud AnnData stores, including `obs` and `var`.
    No array data should need to be read into memory with the exceptio of Awkward Arrays and some older-encoding string arrays.

    Params
    ------
        store: A store-like object to be read in.  If :doc:`zarr:index`, it is best
        for it to be consolidated.

    Returns
    -------
        A lazily read-in AnnData object.
    """
    is_h5 = False
    if isinstance(store, Path) or isinstance(store, str):
        store = str(store)
        if store.endswith("h5ad"):
            is_h5 = True

    has_keys = True  # true if consolidated or h5ad
    if not is_h5:
        import zarr

        try:
            f = zarr.open_consolidated(store, mode="r")
        except KeyError:
            warnings.warn(
                "Did not read zarr as consolidated. Consider consolidating your metadata."
            )
            has_keys = False
            f = zarr.open(store, mode="r")
    else:
        f = h5py.File(store, mode="r")

    def callback(func, elem_name: str, elem, iospec, dataset_kwargs):
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
            index = d[
                elem.attrs["_index"]
            ]  # no sense in reading this in multiple times
            for k in d:
                v = d[k]
                if type(v) == DaskArray and k != elem.attrs["_index"]:
                    d_with_xr[k] = xr.DataArray(
                        v, coords=[index], dims=[index_label], name=k
                    )
                elif (
                    type(v) == CategoricalArray or type(v) == MaskedArray
                ) and k != elem.attrs["_index"]:
                    variable = xr.Variable(
                        data=xr.core.indexing.LazilyIndexedArray(v), dims=[index_label]
                    )
                    d_with_xr[k] = xr.DataArray(
                        variable,
                        coords=[index],
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
            return CategoricalArray(
                codes=elem["codes"],
                categories=elem["categories"],
                ordered=elem.attrs["ordered"],
                drop_unused_cats=drop_unused_cats,
            )
        elif "nullable" in iospec.encoding_type:
            return MaskedArray(
                values=elem["values"],
                mask=elem["mask"] if "mask" in elem else None,
                dtype_str=iospec.encoding_type,
            )
        elif iospec.encoding_type in {
            "csr_matrix",
            "csc_matrix",
            "array",
            "string-array",
        }:
            chunks = None
            if "csr_matrix" == iospec.encoding_type:
                chunks = (1, elem.attrs["shape"][1])
            elif iospec.encoding_type == "csc_matrix":
                chunks = (elem.attrs["shape"][0], 1)
            return read_elem_as_dask(elem, chunks=chunks)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata