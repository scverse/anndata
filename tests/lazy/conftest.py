from __future__ import annotations

import typing
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData
from anndata._types import AnnDataElem
from anndata.experimental import read_lazy
from anndata.tests.helpers import (
    DEFAULT_COL_TYPES,
    DEFAULT_KEY_TYPES,
    AccessTrackingStore,
    AwkArray,
    as_dense_dask_array,
    gen_adata,
    gen_typed_df,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path
    from typing import Literal

ANNDATA_ELEMS = typing.get_args(AnnDataElem)


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
    scope="session",
)
def mtx_format(request):
    return request.param


@pytest.fixture(
    params=[True, False], ids=["vars_different", "vars_same"], scope="session"
)
def are_vars_different(request):
    return request.param


@pytest.fixture(params=["zarr", "h5ad"], scope="session")
def diskfmt(request) -> Literal["zarr", "h5ad"]:
    return request.param


@pytest.fixture(
    params=[True, False],
    scope="session",
    ids=["load-annotation-index", "dont-load-annotation-index"],
)
def load_annotation_index(request):
    return request.param


@pytest.fixture(params=["outer", "inner"], scope="session")
def join(request):
    return request.param


@pytest.fixture(
    params=[
        pytest.param(lambda x: x, id="full"),
        pytest.param(lambda x: x[0:10, :], id="subset"),
    ],
    scope="session",
)
def simple_subset_func(request):
    return request.param


@pytest.fixture(scope="session")
def adata_remote_orig_with_path(
    tmp_path_factory,
    diskfmt: str,
    mtx_format,
    worker_id: str = "serial",
) -> tuple[Path, AnnData]:
    """Create remote fixtures, one without a range index and the other with"""
    file_name = f"orig_{worker_id}.{diskfmt}"
    if diskfmt == "h5ad":
        orig_path = tmp_path_factory.mktemp("h5ad_file_dir") / file_name
    else:
        orig_path = tmp_path_factory.mktemp(file_name)
    orig = gen_adata(
        (100, 110),
        mtx_format,
        obs_dtypes=(*DEFAULT_COL_TYPES, pd.StringDtype),
        var_dtypes=(*DEFAULT_COL_TYPES, pd.StringDtype),
        obsm_types=(*DEFAULT_KEY_TYPES, AwkArray),
        varm_types=(*DEFAULT_KEY_TYPES, AwkArray),
    )
    orig.raw = orig.copy()
    with ad.settings.override(allow_write_nullable_strings=True):
        getattr(ad.io, f"write_{diskfmt}")(
            orig_path, orig, convert_strings_to_categoricals=False
        )
    return orig_path, orig


@pytest.fixture
def adata_remote(
    adata_remote_orig_with_path: tuple[Path, AnnData], *, load_annotation_index: bool
) -> AnnData:
    orig_path, _ = adata_remote_orig_with_path
    return read_lazy(orig_path, load_annotation_index=load_annotation_index)


@pytest.fixture
def adata_orig(adata_remote_orig_with_path: tuple[Path, AnnData]) -> AnnData:
    _, orig = adata_remote_orig_with_path
    return orig


@pytest.fixture(scope="session")
def adata_remote_with_store_tall_skinny_path(
    tmp_path_factory,
    mtx_format,
    worker_id: str = "serial",
) -> Path:
    orig_path = tmp_path_factory.mktemp(f"orig_{worker_id}.zarr")
    M = 100_000  # forces zarr to chunk `obs` columns multiple ways - that way 1 access to `int64` below is actually only one access
    N = 5
    obs_names = pd.Index(f"cell{i}" for i in range(M))
    var_names = pd.Index(f"gene{i}" for i in range(N))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    orig = AnnData(
        obs=obs,
        var=var,
        X=mtx_format(np.random.binomial(100, 0.005, (M, N)).astype(np.float32)),
    )
    orig.raw = orig.copy()
    orig.write_zarr(orig_path)
    return orig_path


@pytest.fixture(scope="session")
def adatas_paths_var_indices_for_concatenation(
    tmp_path_factory, *, are_vars_different: bool, worker_id: str = "serial"
) -> tuple[list[AnnData], list[Path], list[pd.Index]]:
    adatas = []
    var_indices = []
    paths = []
    M = 1000
    N = 50
    n_datasets = 3
    for dataset_index in range(n_datasets):
        orig_path = tmp_path_factory.mktemp(f"orig_{worker_id}_{dataset_index}.zarr")
        paths.append(orig_path)
        obs_names = pd.Index(f"cell_{dataset_index}_{i}" for i in range(M))
        var_names = pd.Index(
            f"gene_{i}{f'_{dataset_index}_ds' if are_vars_different and (i % 2) else ''}"
            for i in range(N)
        )
        var_indices.append(var_names)
        obs = gen_typed_df(M, obs_names)
        var = gen_typed_df(N, var_names)
        orig = AnnData(
            obs=obs,
            var=var,
            X=np.random.binomial(100, 0.005, (M, N)).astype(np.float32),
        )
        orig.write_zarr(orig_path)
        adatas.append(orig)
    return adatas, paths, var_indices


@pytest.fixture
def var_indices_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[pd.Index]:
    _, _, var_indices = adatas_paths_var_indices_for_concatenation
    return var_indices


@pytest.fixture
def adatas_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[AnnData]:
    adatas, _, _ = adatas_paths_var_indices_for_concatenation
    return adatas


@pytest.fixture
def stores_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[AccessTrackingStore]:
    _, paths, _ = adatas_paths_var_indices_for_concatenation
    return [AccessTrackingStore(path) for path in paths]


@pytest.fixture
def lazy_adatas_for_concat(
    stores_for_concat,
) -> list[AnnData]:
    return [read_lazy(store) for store in stores_for_concat]


@pytest.fixture
def adata_remote_with_store_tall_skinny(
    adata_remote_with_store_tall_skinny_path: Path,
) -> tuple[AnnData, AccessTrackingStore]:
    store = AccessTrackingStore(adata_remote_with_store_tall_skinny_path)
    remote = read_lazy(store)
    return remote, store


@pytest.fixture
def remote_store_tall_skinny(
    adata_remote_with_store_tall_skinny_path: Path,
) -> AccessTrackingStore:
    return AccessTrackingStore(adata_remote_with_store_tall_skinny_path)


@pytest.fixture
def adata_remote_tall_skinny(
    remote_store_tall_skinny: AccessTrackingStore,
) -> AnnData:
    remote = read_lazy(remote_store_tall_skinny)
    return remote


def get_key_trackers_for_columns_on_axis(
    adata: AnnData, axis: Literal["obs", "var"]
) -> Generator[str, None, None]:
    """Generate keys for tracking, using `codes` from categorical columns instead of the column name

    Parameters
    ----------
    adata
        Object to get keys from
    axis
        Axis to get keys from

    Yields
    ------
    Keys for tracking
    """
    for col in getattr(adata, axis).columns:
        yield f"{axis}/{col}" if "cat" not in col else f"{axis}/{col}/codes"


ANNDATA_ELEMS = typing.get_args(AnnDataElem)
