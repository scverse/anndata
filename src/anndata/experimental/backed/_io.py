from __future__ import annotations

import warnings
from pathlib import Path
from typing import (
    TYPE_CHECKING,
)

import h5py

from anndata._io.specs.registry import read_elem_lazy

from ..._core.anndata import AnnData
from .. import read_dispatched

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from ...compat import ZarrGroup


def read_backed(
    store: str | Path | MutableMapping | ZarrGroup | h5py.Dataset,
) -> AnnData:
    """
    Lazily read in on-disk/in-cloud AnnData stores, including `obs` and `var`.
    No array data should need to be read into memory with the exception of :class:`awkward.Array` and some older-encoding string arrays.

    Parameters
    ----------
    store
        A store-like object to be read in.  If :class:`zarr.hierarchy.Group`, it is best for it to be consolidated.

    Returns
    -------
        A lazily read-in :class:`~anndata.AnnData` object.
    """
    is_h5 = isinstance(store, (Path, str)) and Path(store).suffix == ".h5ad"

    has_keys = True  # true if consolidated or h5ad
    if not is_h5:
        import zarr

        try:
            f = zarr.open_consolidated(store, mode="r")
        except KeyError:
            msg = "Did not read zarr as consolidated. Consider consolidating your metadata."
            warnings.warn(msg)
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
        elif (
            iospec.encoding_type
            in {
                "csr_matrix",
                "csc_matrix",
                "array",
                "string-array",
                "dataframe",
                "categorical",
            }
            or "nullable" in iospec.encoding_type
        ):
            return read_elem_lazy(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
