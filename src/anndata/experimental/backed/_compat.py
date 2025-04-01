from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

if find_spec("xarray") or TYPE_CHECKING:
    import xarray
    from xarray import DataArray
    from xarray.backends import BackendArray
    from xarray.backends.zarr import ZarrArrayWrapper


else:

    class DataArray:
        def __repr__(self) -> str:
            return "mock DataArray"

    xarray = None

    class ZarrArrayWrapper:
        def __repr__(self) -> str:
            return "mock ZarrArrayWrapper"

    class BackendArray:
        def __repr__(self) -> str:
            return "mock BackendArray"


from ._xarray import Dataset, Dataset2D  # noqa: F401

if TYPE_CHECKING:
    from anndata import AnnData


def has_dataset_2d(adata: AnnData) -> bool:
    if any(isinstance(annot_df, Dataset2D) for annot_df in [adata.obs, adata.var]):
        return True
    for annot_m_key in ["varm", "obsm"]:
        annot_m = getattr(adata, annot_m_key)
        if any(isinstance(maybe_df, Dataset2D) for maybe_df in annot_m.values()):
            return True
    return False
