from __future__ import annotations

from contextlib import nullcontext
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, cast

import h5py
import zarr

from anndata._core.file_backing import AnnDataFileManager
from anndata._io.specs.registry import read_elem, read_elem_lazy
from anndata._types import AnnDataElem
from testing.anndata._doctest import doctest_needs

from ..._core.anndata import AnnData
from ..._core.xarray import requires_xarray
from ..._settings import settings
from ...utils import get_literal_members, warn
from .. import read_dispatched

if TYPE_CHECKING:
    from collections.abc import MutableMapping
    from typing import Any

    from zarr.storage import StoreLike

    from anndata._io.specs.registry import IOSpec
    from anndata._types import Read, StorageType, _GroupStorageType


ANNDATA_ELEMS: tuple[AnnDataElem, ...] = get_literal_members(AnnDataElem)


def _as_group(elem: StorageType) -> _GroupStorageType:
    if not isinstance(elem, zarr.Group | h5py.Group):
        msg = f"Expected a group, got {type(elem)}"
        raise ValueError(msg)
    return elem


@doctest_needs("xarray")
@requires_xarray
def read_lazy(
    store: PathLike[str]
    | str
    | MutableMapping
    | StoreLike
    | zarr.Group
    | h5py.File
    | h5py.Group,
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
        If a path to an ``.h5ad`` file is provided, the open HDF5 file will be attached to the {class}`~anndata.AnnData` at the `file` attribute and it will be the user’s responsibility to close it when done with the returned object.
        For this reason, it is recommended to use an {class}`h5py.File` as the `store` argument when working with h5 files.
        It must remain open for at least as long as this returned object is in use.
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
    >>> import pooch
    >>> import scanpy as sc
    >>> base_url = "https://datasets.cellxgene.cziscience.com"
    >>> # To update hashes: pooch.retrieve(url, known_hash=None) prints the new hash
    >>> def get_cellxgene_data(id_: str, hash_: str):
    ...     return pooch.retrieve(
    ...         f"{base_url}/{id_}.h5ad",
    ...         known_hash=hash_,
    ...         fname=f"{id_}.h5ad",
    ...         path=sc.settings.datasetdir,
    ...     )
    >>> path_b_cells = get_cellxgene_data(
    ...     "a93eab58-3d82-4b61-8a2f-d7666dcdb7c4",
    ...     "sha256:dac90fe2aa8b78aee2c1fc963104592f8eff7b873ca21d01a51a5e416734651c",
    ... )
    >>> path_fetal = get_cellxgene_data(
    ...     "d170ff04-6da0-4156-a719-f8e1bbefbf53",
    ...     "sha256:d497eebca03533919877b6fc876e8c9d8ba063199ddc86dd9fbcb9d1d87a3622",
    ... )
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
    f: StorageType
    h5_file: h5py.File | None = None
    has_keys = True  # true if consolidated or h5ad
    if isinstance(store, h5py.Dataset | h5py.File | h5py.Group):
        f = store
    elif isinstance(store, PathLike | str) and Path(store).suffix == ".h5ad":
        f = h5_file = h5py.File(store, mode="r")
    elif isinstance(store, zarr.Group):
        f = store
    else:
        try:
            f = zarr.open_consolidated(store, mode="r")
        except ValueError:
            msg = "Did not read zarr as consolidated. Consider consolidating your metadata."
            warn(msg, UserWarning)
            has_keys = False
            f = zarr.open_group(store, mode="r")

    def callback(func: Read, /, elem_name: str, elem: StorageType, *, iospec: IOSpec):
        if iospec.encoding_type in {"anndata", "raw"} or elem_name.endswith("/"):
            group = _as_group(elem)
            iter_object = (
                dict(group).items()
                if has_keys
                else tuple(
                    (k, v)
                    for k, v in ((k, group.get(k, None)) for k in ANNDATA_ELEMS)
                    # need to do this instead of `k in elem` to prevent unnecessary metadata accesses
                    if v is not None
                )
            )
            kwargs: dict[str, Any] = {
                k: read_dispatched(v, callback) for k, v in iter_object
            }
            return AnnData(**kwargs)
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
            return read_elem(elem)
        elif iospec.encoding_type == "dict":
            return {
                k: read_dispatched(v, callback=callback)
                for k, v in dict(_as_group(elem)).items()
            }
        return func(elem)

    with (
        settings.override(check_uniqueness=False)
        if not load_annotation_index
        else nullcontext()
    ):
        adata = cast("AnnData", read_dispatched(f, callback=callback))
    if h5_file is not None:
        adata.file = AnnDataFileManager(adata, file_obj=h5_file)
    return adata
