from __future__ import annotations

import typing
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import h5py

from anndata._io.specs.registry import read_elem_lazy
from anndata._types import AnnDataElem

from ..._core.anndata import AnnData
from ..._settings import settings
from .. import read_dispatched

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from anndata._io.specs.registry import IOSpec
    from anndata._types import Read, StorageType

    from ...compat import ZarrGroup


def read_lazy(
    store: str | Path | MutableMapping | ZarrGroup | h5py.Dataset,
    *,
    load_annotation_index: bool = True,
) -> AnnData:
    """
    Lazily read in on-disk/in-cloud AnnData stores, including `obs` and `var`.
    No array data should need to be read into memory with the exception of :class:`ak.Array`, scalars, and some older-encoding arrays.

    Parameters
    ----------
    store
        A store-like object to be read in.  If :class:`zarr.hierarchy.Group`, it is best for it to be consolidated.
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
        obs: 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sample_uuid', 'sample_preservation_method', 'tissue_ontology_term_id', 'development_stage_ontology_term_id', 'suspension_uuid', 'suspension_type', 'library_uuid', 'assay_ontology_term_id', 'mapped_reference_annotation', 'is_primary_data', 'cell_type_ontology_term_id', 'author_cell_type', 'disease_ontology_term_id', 'sex_ontology_term_id', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'Phase', 'sample', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype'
        uns: 'default_embedding', 'schema_version', 'title'
        obsm: 'X_harmony', 'X_pca', 'X_umap'
    >>> print(fetal_adata)
    AnnData object with n_obs × n_vars = 344 × 15585
        obs: 'nCount_Spatial', 'nFeature_Spatial', 'Cluster', 'adult_pred_type', 'adult_pred_value', 'fetal_pred_type', 'fetal_pred_value', 'pDCs', 'Cell Cycle', 'Type 3 ILCs', 'DCs', 'Mast', 'Monocytes', 'Naive T-Cells', 'Venous (CP) 1', 'Venous (M) 2', 'Arterial (L)', 'Endothelium G2M-phase', 'Venous (CP) 2', 'Arterial (CP)', 'Arterial (M)', 'Endothelium S-phase', 'Proximal Progenitor', 'Proximal Mature Enterocytes', 'BEST4_OTOP2 Cells', 'Proximal TA', 'Proximal Early Enterocytes', 'Proximal Enterocytes', 'Proximal Stem Cells', 'EECs', 'Distal Enterocytes', 'Goblets', 'Distal TA', 'Distal Absorptive', 'Distal Stem Cells', 'Secretory Progenitors', 'Distal Mature Enterocytes', 'S1', 'S1 COL6A5+', 'S4 CCL21+', 'Proximal S2 (2)', 'S1 IFIT3+', 'Distal S2', 'Fibroblasts S-phase', 'Proximal S2 (1)', 'S3 Progenitor', 'Fibroblasts G2M-phase', 'S4 CXCL14+', 'Fibroblast Progenitor', 'S3 Transitional', 'Erythroid', 'S3 EBF+', 'S3 HAND1+', 'Pericytes G2M-phase', 'Pericyte Progenitors', 'Undifferentiated Pericytes', 'ICC PDGFRA+', 'MYOCD+ Muscularis', 'Muscularis S-phase', 'Muscularis G2M-phase', 'HOXP+ Proximal Muscularis', 'FOXF2+ Distal Muscularis', 'FOXF2- Muscularis', 'MORN5+ Distal Muscularis', 'Myofibroblast Progenitors', 'Myofibroblasts', 'Mesothelium SOX6+', 'Myofibroblasts S-phase', 'Myofibroblasts G2M-phase', 'Glial Progenitors', 'Excitory Motor Neuron', 'Interneuron', 'Differentiating Submucosal Glial', 'Inhibitory Motor Neuron Precursor', 'Neuroendocrine (1)', 'max', 'tissue_ontology_term_id', 'assay_ontology_term_id', 'disease_ontology_term_id', 'development_stage_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'cell_type_ontology_term_id', 'sex_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id', 'donor_id', 'suspension_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage'
        var: 'sct.detection_rate', 'sct.gmean', 'sct.variance', 'sct.residual_mean', 'sct.residual_variance', 'sct.variable', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype'
        uns: 'adult_pred_mat', 'fetal_pred_mat', 'schema_version', 'title'
        obsm: 'X_pca', 'X_spatial', 'X_umap'
        layers: 'counts', 'scale.data'

    This functionality is compatible with :func:`anndata.concat`

    >>> ad.concat([b_cells_adata, fetal_adata], join="outer")
    AnnData object with n_obs × n_vars = 490 × 33452
        obs: 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sample_uuid', 'sample_preservation_method', 'tissue_ontology_term_id', 'development_stage_ontology_term_id', 'suspension_uuid', 'suspension_type', 'library_uuid', 'assay_ontology_term_id', 'mapped_reference_annotation', 'is_primary_data', 'cell_type_ontology_term_id', 'author_cell_type', 'disease_ontology_term_id', 'sex_ontology_term_id', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'Phase', 'sample', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'nCount_Spatial', 'nFeature_Spatial', 'Cluster', 'adult_pred_type', 'adult_pred_value', 'fetal_pred_type', 'fetal_pred_value', 'pDCs', 'Cell Cycle', 'Type 3 ILCs', 'DCs', 'Mast', 'Monocytes', 'Naive T-Cells', 'Venous (CP) 1', 'Venous (M) 2', 'Arterial (L)', 'Endothelium G2M-phase', 'Venous (CP) 2', 'Arterial (CP)', 'Arterial (M)', 'Endothelium S-phase', 'Proximal Progenitor', 'Proximal Mature Enterocytes', 'BEST4_OTOP2 Cells', 'Proximal TA', 'Proximal Early Enterocytes', 'Proximal Enterocytes', 'Proximal Stem Cells', 'EECs', 'Distal Enterocytes', 'Goblets', 'Distal TA', 'Distal Absorptive', 'Distal Stem Cells', 'Secretory Progenitors', 'Distal Mature Enterocytes', 'S1', 'S1 COL6A5+', 'S4 CCL21+', 'Proximal S2 (2)', 'S1 IFIT3+', 'Distal S2', 'Fibroblasts S-phase', 'Proximal S2 (1)', 'S3 Progenitor', 'Fibroblasts G2M-phase', 'S4 CXCL14+', 'Fibroblast Progenitor', 'S3 Transitional', 'Erythroid', 'S3 EBF+', 'S3 HAND1+', 'Pericytes G2M-phase', 'Pericyte Progenitors', 'Undifferentiated Pericytes', 'ICC PDGFRA+', 'MYOCD+ Muscularis', 'Muscularis S-phase', 'Muscularis G2M-phase', 'HOXP+ Proximal Muscularis', 'FOXF2+ Distal Muscularis', 'FOXF2- Muscularis', 'MORN5+ Distal Muscularis', 'Myofibroblast Progenitors', 'Myofibroblasts', 'Mesothelium SOX6+', 'Myofibroblasts S-phase', 'Myofibroblasts G2M-phase', 'Glial Progenitors', 'Excitory Motor Neuron', 'Interneuron', 'Differentiating Submucosal Glial', 'Inhibitory Motor Neuron Precursor', 'Neuroendocrine (1)', 'max'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'sct.detection_rate', 'sct.gmean', 'sct.variance', 'sct.residual_mean', 'sct.residual_variance', 'sct.variable'
        obsm: 'X_harmony', 'X_pca', 'X_umap', 'X_spatial'
        layers: 'counts', 'scale.data'
    """
    try:
        import xarray  # noqa: F401
    except ImportError:
        msg = (
            "xarray is required to use the `read_lazy` function. Please install xarray."
        )
        raise ImportError(msg)
    is_h5_store = isinstance(store, h5py.Dataset | h5py.File | h5py.Group)
    is_h5 = (
        isinstance(store, Path | str) and Path(store).suffix == ".h5ad"
    ) or is_h5_store

    has_keys = True  # true if consolidated or h5ad
    if not is_h5:
        import zarr

        if not isinstance(store, zarr.hierarchy.Group):
            try:
                f = zarr.open_consolidated(store, mode="r")
            except KeyError:
                msg = "Did not read zarr as consolidated. Consider consolidating your metadata."
                warnings.warn(msg)
                has_keys = False
                f = zarr.open(store, mode="r")
        else:
            f = store
    else:
        if is_h5_store:
            f = store
        else:
            f = h5py.File(store, mode="r")

    def callback(func: Read, /, elem_name: str, elem: StorageType, *, iospec: IOSpec):
        if iospec.encoding_type in {"anndata", "raw"} or elem_name.endswith("/"):
            iter_object = (
                elem.items()
                if has_keys
                else [(k, elem[k]) for k in typing.get_args(AnnDataElem) if k in elem]
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
            if "dataframe" == iospec.encoding_type and elem_name in {"/obs", "/var"}:
                return read_elem_lazy(elem, use_range_index=not load_annotation_index)
            return read_elem_lazy(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        elif iospec.encoding_type == "dict":
            return {k: read_dispatched(v, callback=callback) for k, v in elem.items()}
        return func(elem)

    with settings.override(check_uniqueness=load_annotation_index):
        adata = read_dispatched(f, callback=callback)

    return adata
