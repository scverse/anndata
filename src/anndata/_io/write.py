from __future__ import annotations

import math
import warnings
from os import fspath
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from anndata._io.utils import no_write_dataset_2d

from .._warnings import WriteWarning
from ..compat import old_positionals
from ..logging import get_logger

if TYPE_CHECKING:
    from os import PathLike

    from .. import AnnData

logger = get_logger(__name__)


@no_write_dataset_2d
@old_positionals("skip_data", "sep")
def write_csvs(
    dirname: PathLike[str] | str,
    adata: AnnData,
    *,
    skip_data: bool = True,
    sep: str = ",",
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
        obsm=adata.obsm.to_df(),
        varm=adata.varm.to_df(),
    )
    if not skip_data:
        d["X"] = pd.DataFrame(adata.X.toarray() if issparse(adata.X) else adata.X)
    d_write = {**d, **adata._uns}
    not_yet_raised_sparse_warning = True
    for key, value in d_write.items():
        if issparse(value):
            if not_yet_raised_sparse_warning:
                msg = "Omitting to write sparse annotation."
                warnings.warn(msg, WriteWarning, stacklevel=2)
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
            except Exception as e:  # noqa: BLE001
                msg = f"Omitting to write {key!r} of type {type(e)}."
                warnings.warn(msg, WriteWarning, stacklevel=2)
                continue
        df.to_csv(
            filename,
            sep=sep,
            header=key in {"obs", "var", "obsm", "varm"},
            index=key in {"obs", "var"},
        )


@no_write_dataset_2d
@old_positionals("write_obsm_varm")
def write_loom(
    filename: PathLike[str] | str, adata: AnnData, *, write_obsm_varm: bool = False
) -> None:
    """See :meth:`~anndata.AnnData.write_loom`."""
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
        msg = "loompy does not accept empty matrices as data"
        raise ValueError(msg)

    if write_obsm_varm:
        col_attrs.update(adata.obsm)
        row_attrs.update(adata.varm)
    elif len(adata.obsm.keys()) > 0 or len(adata.varm.keys()) > 0:
        logger.warning(
            f"The loom file will lack these fields:\n"
            f"{adata.obsm.keys() | adata.varm.keys()}\n"
            f"Use write_obsm_varm=True to export multi-dimensional annotations"
        )

    layers = {"": adata.X.T}
    for key, layer in adata.layers.items():
        layers[key] = layer.T

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
        for i in range(math.ceil(float(za.shape[0]) / za.chunks[0]))
        for j in range(math.ceil(float(za.shape[1]) / za.chunks[1]))
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
