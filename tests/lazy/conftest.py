from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse
from zarr.storage import MemoryStore

import anndata as ad
from anndata import AnnData
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
    from pathlib import Path
    from typing import Literal


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
def adata_remote_orig_with_store(
    tmp_path_factory,
    diskfmt: str,
    mtx_format,
    worker_id: str = "serial",
) -> tuple[Path | MemoryStore, AnnData]:
    """Create remote fixtures, one without a range index and the other with"""
    if diskfmt == "h5ad":
        orig_store = tmp_path_factory.mktemp("h5ad_file_dir") / f"orig_{worker_id}.h5ad"
    else:
        orig_store = MemoryStore()
    orig = gen_adata(
        (100, 110),
        mtx_format,
        obs_dtypes=(*DEFAULT_COL_TYPES, pd.StringDtype()),
        var_dtypes=(*DEFAULT_COL_TYPES, pd.StringDtype()),
        obsm_types=(*DEFAULT_KEY_TYPES, AwkArray),
        varm_types=(*DEFAULT_KEY_TYPES, AwkArray),
    )
    orig.raw = orig.copy()
    with ad.settings.override(allow_write_nullable_strings=True):
        getattr(ad.io, f"write_{diskfmt}")(
            orig_store, orig, convert_strings_to_categoricals=False
        )
    return orig_store, orig


@pytest.fixture
def adata_remote(
    adata_remote_orig_with_store: tuple[Path | MemoryStore, AnnData],
    *,
    load_annotation_index: bool,
) -> AnnData:
    orig_store, _ = adata_remote_orig_with_store
    return read_lazy(orig_store, load_annotation_index=load_annotation_index)


@pytest.fixture
def adata_orig(
    adata_remote_orig_with_store: tuple[Path | MemoryStore, AnnData],
) -> AnnData:
    _, orig = adata_remote_orig_with_store
    return orig


@pytest.fixture(scope="session", params=[pytest.param(None, marks=pytest.mark.zarr_io)])
def adata_remote_tall_skinny_orig_store(mtx_format) -> MemoryStore:
    orig_store = MemoryStore()
    M = 1000
    N = 5
    # One named, one unnamed
    obs_names = pd.Index((f"cell{i}" for i in range(M)), name="obs_names")
    var_names = pd.Index(f"gene{i}" for i in range(N))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    orig = AnnData(
        obs=obs,
        var=var,
        X=mtx_format(np.random.binomial(100, 0.005, (M, N)).astype(np.float32)),
    )
    orig.raw = orig.copy()
    orig.write_zarr(orig_store)
    g = zarr.open_group(orig_store, mode="a", use_consolidated=False)
    ad.io.write_elem(
        g,
        "obs",
        obs,
        # No shards so we can track chunking exactly.
        dataset_kwargs=dict(chunks=(250,), shards=None),
    )
    # Catch the warning so we are alerted once it is no longer surfaced i.e., once consolidated metadata stabilizes.
    with pytest.warns(
        zarr.errors.ZarrUserWarning if hasattr(zarr, "errors") else UserWarning,
        match=r"Consolidated metadata",
    ):
        zarr.consolidate_metadata(g.store)
    return orig_store


@pytest.fixture(scope="session", params=[pytest.param(None, marks=pytest.mark.zarr_io)])
def adatas_stores_var_indices_for_concatenation(
    *, are_vars_different: bool
) -> tuple[list[AnnData], list[MemoryStore], list[pd.Index]]:
    adatas = []
    var_indices = []
    stores = []
    M = 1000
    N = 50
    n_datasets = 3
    for dataset_index in range(n_datasets):
        orig_store = MemoryStore()
        stores.append(orig_store)
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
        orig.write_zarr(orig_store)
        adatas.append(orig)
    return adatas, stores, var_indices


@pytest.fixture
def var_indices_for_concat(
    adatas_stores_var_indices_for_concatenation,
) -> list[pd.Index]:
    _, _, var_indices = adatas_stores_var_indices_for_concatenation
    return var_indices


@pytest.fixture
def adatas_for_concat(
    adatas_stores_var_indices_for_concatenation,
) -> list[AnnData]:
    adatas, _, _ = adatas_stores_var_indices_for_concatenation
    return adatas


@pytest.fixture
def stores_for_concat(
    adatas_stores_var_indices_for_concatenation,
) -> list[AccessTrackingStore]:
    _, stores, _ = adatas_stores_var_indices_for_concatenation
    return [AccessTrackingStore(store) for store in stores]


@pytest.fixture
def lazy_adatas_for_concat(
    stores_for_concat: list[AccessTrackingStore], *, load_annotation_index: bool
) -> list[AnnData]:
    return [
        read_lazy(store, load_annotation_index=load_annotation_index)
        for store in stores_for_concat
    ]


@pytest.fixture
def adata_remote_with_store_tall_skinny(
    adata_remote_tall_skinny_orig_store: MemoryStore,
) -> tuple[AnnData, AccessTrackingStore]:
    store = AccessTrackingStore(adata_remote_tall_skinny_orig_store)
    remote = read_lazy(store)
    return remote, store


@pytest.fixture
def remote_store_tall_skinny(
    adata_remote_tall_skinny_orig_store: MemoryStore,
) -> AccessTrackingStore:
    return AccessTrackingStore(adata_remote_tall_skinny_orig_store)


@pytest.fixture
def adata_remote_tall_skinny(
    remote_store_tall_skinny: AccessTrackingStore,
) -> AnnData:
    remote = read_lazy(remote_store_tall_skinny)
    return remote
