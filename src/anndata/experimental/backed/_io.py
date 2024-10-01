from __future__ import annotations

import warnings
from pathlib import Path
from typing import (
    TYPE_CHECKING,
)

import h5py

from anndata._io.specs.registry import read_elem_lazy

from ..._core.anndata import AnnData
from ..._settings import settings
from .. import read_dispatched

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from ...compat import ZarrGroup


def read_lazy(
    store: str | Path | MutableMapping | ZarrGroup | h5py.Dataset,
    load_annotation_index: bool = True,
) -> AnnData:
    """
    Lazily read in on-disk/in-cloud AnnData stores, including `obs` and `var`.
    No array data should need to be read into memory with the exception of :class:`awkward.Array` and some older-encoding string arrays.

    Parameters
    ----------
    store
        A store-like object to be read in.  If :class:`zarr.hierarchy.Group`, it is best for it to be consolidated.
    load_annotation_index
        Whether or not to use a range index for the `{obs,var}` :class:`xr.Dataset` so as not to load the index into memory.
        If `False`, the real `index` will be inserted as `{obs,var}_names` in the object but not be one of the `coords` thereby preventing read operations.

    Returns
    -------
        A lazily read-in :class:`~anndata.AnnData` object.
    """
    is_h5_store = isinstance(store, (h5py.Dataset, h5py.File))
    is_h5 = (
        isinstance(store, (Path, str)) and Path(store).suffix == ".h5ad"
    ) or is_h5_store

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
        if is_h5_store:
            f = store
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
            if "dataframe" == iospec.encoding_type and elem_name in {"/obs", "/var"}:
                return read_elem_lazy(elem, use_range_index=not load_annotation_index)
            return read_elem_lazy(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        return func(elem)

    with settings.override(check_uniqueness=load_annotation_index):
        adata = read_dispatched(f, callback=callback)

    return adata