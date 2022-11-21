from ctypes import Union
from typing import Any, Callable
import warnings
from pathlib import Path
from os import PathLike, fspath

import pandas as pd
import math
import numpy as np
from scipy.sparse import issparse
from anndata._core.anndata import StorageType

from anndata._io.read import GroupType

from .. import AnnData
from ..logging import get_logger
from anndata._warnings import WriteWarning
from anndata._io.specs.registry import _REGISTRY, IORegistry, IOSpec, write_elem

from ..utils import import_function

logger = get_logger(__name__)

def write_dispatched(
    group: GroupType,
    adata: AnnData,
    dispatch_element: Callable[
        [Callable[[GroupType, str, Any, dict], Any], GroupType, str, Any], Any
    ],
    add_element_dispatcher: Callable[
        [Callable[[IORegistry, StorageType, IOSpec, frozenset], Any]], Any
    ] = lambda x: x,
) -> AnnData:
    """
    Custom reading of `filename` via `dispatcher` - limited to zarr and h5ad.
    group
        Data group.
    adata
        AnnData object
    dispatch_element
        Function for writing data from AnnData slots onto disk.
    add_element_dispatcher:
        Function for adding your own custom writer to AnnData's core, including overriding current methods.
    """
    add_element_dispatcher(_REGISTRY.register_write)
    keys = ["X", "raw", "uns", "layers", "var", "obs", "varm", "obsm", "varp", "obsp"]
    dict_keys = set(["obsm", "varm", "obsp", "varp", "layers", "uns"])
    for key in keys:
        elem = getattr(adata, key)
        if key in dict_keys:
            elem = dict(elem)
        dispatch_element(write_elem, group, key, elem)


def write_csvs(
    dirname: PathLike, adata: AnnData, skip_data: bool = True, sep: str = ","
):
    """See :meth:`~anndata.AnnData.write_csvs`."""
    dirname = Path(dirname)
    if dirname.suffix == ".csv":
        dirname = dirname.with_suffix("")
    logger.info(f"writing .csv files to {dirname}")
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=True)
    dir_uns = dirname / "uns"
    if not dir_uns.is_dir():
        dir_uns.mkdir(parents=True, exist_ok=True)
    d = dict(
        obs=adata._obs,
        var=adata._var,
        obsm=adata._obsm.to_df(),
        varm=adata._varm.to_df(),
    )
    if not skip_data:
        d["X"] = pd.DataFrame(adata.X.toarray() if issparse(adata.X) else adata.X)
    d_write = {**d, **adata._uns}
    not_yet_raised_sparse_warning = True
    for key, value in d_write.items():
        if issparse(value):
            if not_yet_raised_sparse_warning:
                warnings.warn("Omitting to write sparse annotation.", WriteWarning)
                not_yet_raised_sparse_warning = False
            continue
        filename = dirname
        if key not in {"X", "var", "obs", "obsm", "varm"}:
            filename = dir_uns
        filename /= f"{key}.csv"
        df = value
        if not isinstance(value, pd.DataFrame):
            value = np.array(value)
            if np.ndim(value) == 0:
                value = value[None]
            try:
                df = pd.DataFrame(value)
            except Exception as e:
                warnings.warn(
                    f"Omitting to write {key!r} of type {type(e)}.",
                    WriteWarning,
                )
                continue
        df.to_csv(
            filename,
            sep=sep,
            header=key in {"obs", "var", "obsm", "varm"},
            index=key in {"obs", "var"},
        )


def write_loom(filename: PathLike, adata: AnnData, write_obsm_varm: bool = False):
    filename = Path(filename)
    row_attrs = {k: np.array(v) for k, v in adata.var.to_dict("list").items()}
    row_names = adata.var_names
    row_dim = row_names.name if row_names.name is not None else "var_names"
    row_attrs[row_dim] = row_names.values
    col_attrs = {k: np.array(v) for k, v in adata.obs.to_dict("list").items()}
    col_names = adata.obs_names
    col_dim = col_names.name if col_names.name is not None else "obs_names"
    col_attrs[col_dim] = col_names.values

    if adata.X is None:
        raise ValueError("loompy does not accept empty matrices as data")

    if write_obsm_varm:
        for key in adata.obsm.keys():
            col_attrs[key] = adata.obsm[key]
        for key in adata.varm.keys():
            row_attrs[key] = adata.varm[key]
    elif len(adata.obsm.keys()) > 0 or len(adata.varm.keys()) > 0:
        logger.warning(
            f"The loom file will lack these fields:\n"
            f"{adata.obsm.keys() | adata.varm.keys()}\n"
            f"Use write_obsm_varm=True to export multi-dimensional annotations"
        )

    layers = {"": adata.X.T}
    for key in adata.layers.keys():
        layers[key] = adata.layers[key].T

    from loompy import create

    if filename.exists():
        filename.unlink()
    create(fspath(filename), layers, row_attrs=row_attrs, col_attrs=col_attrs)


def _get_chunk_indices(za):
    # TODO: does zarr provide code for this?
    """\
    Return all the indices (coordinates) for the chunks in a zarr array,
    even empty ones.
    """
    return [
        (i, j)
        for i in range(int(math.ceil(float(za.shape[0]) / za.chunks[0])))
        for j in range(int(math.ceil(float(za.shape[1]) / za.chunks[1])))
    ]


def _write_in_zarr_chunks(za, key, value):
    if key != "X":
        za[:] = value  # don’t chunk metadata
    else:
        for ci in _get_chunk_indices(za):
            s0, e0 = za.chunks[0] * ci[0], za.chunks[0] * (ci[0] + 1)
            s1, e1 = za.chunks[1] * ci[1], za.chunks[1] * (ci[1] + 1)
            print(ci, s0, e1, s1, e1)
            if issparse(value):
                za[s0:e0, s1:e1] = value[s0:e0, s1:e1].todense()
            else:
                za[s0:e0, s1:e1] = value[s0:e0, s1:e1]
