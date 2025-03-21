from __future__ import annotations

import atexit
from functools import partial
from typing import TYPE_CHECKING, cast

import dask
import joblib
import pytest
from dask.base import normalize_token, tokenize
from packaging.version import Version

if Version(dask.__version__) < Version("2024.8.0"):
    from dask.base import normalize_seq
else:
    from dask.tokenize import normalize_seq
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import subset_func  # noqa: F401

if TYPE_CHECKING:
    from types import EllipsisType

    from xdist.workermanage import WorkerController


_dask_cluster_addr: str


def pytest_configure(config: pytest.Config) -> None:
    # We use this hook because it is run on sequential sessions,
    # and both pytest-xdist’s controller and workers
    global _dask_cluster_addr  # noqa: PLW0603

    if "--collect-only" in config.args:
        return  # no need to do work

    _dask_cluster_addr = _get_cluster_address(config)


def pytest_configure_node(node: WorkerController) -> None:
    # send the cluster address to the workers
    node.workerinput["dask_cluster_addr"] = _dask_cluster_addr


def _get_cluster_address(config: pytest.Config) -> str:
    """Start the dask cluster or (in a pytest-xdist worker) get its address."""
    # If we’re on a worker, we can use the data sent in `pytest_configure_node` above
    if workerinput := cast("dict[str, str]", getattr(config, "workerinput", {})):
        return workerinput["dask_cluster_addr"]

    # if we’re on the controller or running sequentially, we start the cluster
    import dask.distributed as dd

    clust = dd.LocalCluster(n_workers=1, threads_per_worker=1)
    clust.__enter__()
    atexit.register(clust.close)
    return clust.scheduler_address


@pytest.fixture(scope="session")
def local_cluster_addr() -> str:
    """Get the dask cluster address"""
    return _dask_cluster_addr


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"


@pytest.fixture(
    params=[
        pytest.param((..., (slice(None), slice(None))), id="ellipsis"),
        pytest.param(((...,), (slice(None), slice(None))), id="ellipsis_tuple"),
        pytest.param(
            ((..., slice(0, 10)), (slice(None), slice(0, 10))), id="obs-ellipsis"
        ),
        pytest.param(
            ((slice(0, 10), ...), (slice(0, 10), slice(None))), id="var-ellipsis"
        ),
        pytest.param(
            ((slice(0, 10), slice(0, 10), ...), (slice(0, 10), slice(0, 10))),
            id="obs-var-ellipsis",
        ),
        pytest.param(
            ((..., slice(0, 10), slice(0, 10)), (slice(0, 10), slice(0, 10))),
            id="ellipsis-obs-var",
        ),
        pytest.param(
            ((slice(0, 10), ..., slice(0, 10)), (slice(0, 10), slice(0, 10))),
            id="obs-ellipsis-var",
        ),
    ]
)
def ellipsis_index_with_equivalent(
    request,
) -> tuple[tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]]:
    return request.param


@pytest.fixture
def ellipsis_index(
    ellipsis_index_with_equivalent: tuple[
        tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]
    ],
) -> tuple[EllipsisType | slice, ...] | EllipsisType:
    return ellipsis_index_with_equivalent[0]


@pytest.fixture
def equivalent_ellipsis_index(
    ellipsis_index_with_equivalent: tuple[
        tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]
    ],
) -> tuple[slice, slice]:
    return ellipsis_index_with_equivalent[1]


#####################
# Dask tokenization #
#####################
# TODO: Should we be exporting this?


# sparray classes don't have tokenize defined yet, see: https://github.com/dask/dask/issues/10375
def normalize_sparse_matrix(x, attrs):
    return (
        type(x).__name__,
        normalize_seq(normalize_token(getattr(x, key)) for key in attrs),
    )


for cls, attrs in [
    (sparse.dia_array, ("data", "offsets", "shape")),
    (sparse.bsr_array, ("data", "indices", "indptr", "blocksize", "shape")),
    (sparse.coo_array, ("data", "row", "col", "shape")),
    (sparse.csr_array, ("data", "indices", "indptr", "shape")),
    (sparse.csc_array, ("data", "indices", "indptr", "shape")),
    (sparse.lil_array, ("data", "rows", "shape")),
]:
    normalize_token.register(cls, partial(normalize_sparse_matrix, attrs=attrs))


@normalize_token.register(sparse.dok_array)
def normalize_dok_matrix(x):
    return type(x).__name__, normalize_token(sorted(x.items()))


@normalize_token.register(ad.AnnData)
def tokenize_anndata(adata: ad.AnnData):
    res = []
    if adata.X is not None:
        res.append(tokenize(adata.X))
    res.extend([tokenize(adata.obs), tokenize(adata.var)])
    for attr in ["obsm", "varm", "obsp", "varp", "layers"]:
        elem = getattr(adata, attr)
        res.append(tokenize(list(dict(elem).items())))
    res.append(joblib.hash(adata.uns))
    if adata.raw is not None:
        res.append(tokenize(adata.raw.to_adata()))
    return tuple(res)
