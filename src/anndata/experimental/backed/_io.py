from __future__ import annotations

import typing
import warnings
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING

import h5py

from anndata._io.specs.registry import read_elem_lazy
from anndata._types import AnnDataElem
from testing.anndata._doctest import doctest_needs

from ..._core.anndata import AnnData
from ..._core.xarray import requires_xarray
from ..._settings import settings
from ...compat import ZarrGroup, is_zarr_v2
from .. import read_dispatched

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from anndata._io.specs.registry import IOSpec
    from anndata._types import Read, StorageType


@doctest_needs("xarray")
@requires_xarray
def read_lazy(
    store: PathLike[str] | str | MutableMapping | ZarrGroup | h5py.Dataset,
    *,
    load_annotation_index: bool = True,
) -> AnnData:
    """
    Lazily read in on-disk/in-cloud AnnData stores, including `obs` and `var`.
    No array data should need to be read into memory with the exception of :class:`ak.Array`, scalars, and some older-encoding arrays.

    Parameters
    ----------
    store
        A store-like object to be read in.  If :class:`zarr.Group`, it is best for it to be consolidated.
    load_annotation_index
        Whether or not to use a range index for the `{obs,var}` :class:`xarray.Dataset` so as not to load the index into memory.
        If `False`, the real `index` will be inserted as `{obs,var}_names` in the object but not be one of the `coords` thereby preventing read operations.
        Access to `adata.obs.index` will also only give the dummy index, and not the "real" index that is file-backed.

    Returns
    -------
    A lazily read-in :class:`~anndata.AnnData` object.

    Examples
    --------

    Preparing example objects

    >>> import anndata as ad
    >>> from urllib.request import urlretrieve
    >>> import scanpy as sc
    >>> base_url = "https://datasets.cellxgene.cziscience.com"
    >>> def get_cellxgene_data(id_: str):
    ...     out_path = sc.settings.datasetdir / f"{id_}.h5ad"
    ...     if out_path.exists():
    ...         return out_path
    ...     file_url = f"{base_url}/{id_}.h5ad"
    ...     sc.settings.datasetdir.mkdir(parents=True, exist_ok=True)
    ...     urlretrieve(file_url, out_path)
    ...     return out_path
    >>> path_b_cells = get_cellxgene_data("a93eab58-3d82-4b61-8a2f-d7666dcdb7c4")
    >>> path_fetal = get_cellxgene_data("d170ff04-6da0-4156-a719-f8e1bbefbf53")
    >>> b_cells_adata = ad.experimental.read_lazy(path_b_cells)
    >>> fetal_adata = ad.experimental.read_lazy(path_fetal)
    >>> print(b_cells_adata)
    AnnData object with n_obs × n_vars = 146 × 33452
        obs: 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id', ...
    >>> print(fetal_adata)
    AnnData object with n_obs × n_vars = 344 × 15585
        obs: 'nCount_Spatial', 'nFeature_Spatial', 'Cluster', 'adult_pred_type'...

    This functionality is compatible with :func:`anndata.concat`

    >>> ad.concat([b_cells_adata, fetal_adata], join="outer")
    AnnData object with n_obs × n_vars = 490 × 33452
        obs: 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id'...
    """
    is_h5_store = isinstance(store, h5py.Dataset | h5py.File | h5py.Group)
    is_h5 = (
        isinstance(store, PathLike | str) and Path(store).suffix == ".h5ad"
    ) or is_h5_store

    has_keys = True  # true if consolidated or h5ad
    if not is_h5:
        import zarr

        if not isinstance(store, ZarrGroup):
            # v3 returns a ValueError for consolidated metadata not found
            err_cls = KeyError if is_zarr_v2() else ValueError
            try:
                f = zarr.open_consolidated(store, mode="r")
            except err_cls:
                msg = "Did not read zarr as consolidated. Consider consolidating your metadata."
                warnings.warn(msg, UserWarning, stacklevel=2)
                has_keys = False
                f = zarr.open_group(store, mode="r")
        else:
            f = store
    elif is_h5_store:
        f = store
    else:
        f = h5py.File(store, mode="r")

    def callback(func: Read, /, elem_name: str, elem: StorageType, *, iospec: IOSpec):
        if iospec.encoding_type in {"anndata", "raw"} or elem_name.endswith("/"):
            iter_object = (
                dict(elem).items()
                if has_keys
                else (
                    (k, v)
                    for k, v in (
                        (k, elem.get(k, None)) for k in typing.get_args(AnnDataElem)
                    )
                    if v
                    is not None  # need to do this instead of `k in elem` to prevent unnecessary metadata accesses
                )
            )
            return AnnData(**{k: read_dispatched(v, callback) for k, v in iter_object})
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
            if iospec.encoding_type == "dataframe" and (
                elem_name[:4] in {"/obs", "/var"}
                or elem_name[:8] in {"/raw/obs", "/raw/var"}
            ):
                return read_elem_lazy(elem, use_range_index=not load_annotation_index)
            return read_elem_lazy(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        elif iospec.encoding_type == "dict":
            return {
                k: read_dispatched(v, callback=callback) for k, v in dict(elem).items()
            }
        return func(elem)

    with settings.override(check_uniqueness=load_annotation_index):
        adata = read_dispatched(f, callback=callback)

    return adata
