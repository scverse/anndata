from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import dask
import joblib
import pytest
from dask.base import normalize_token, tokenize
from packaging.version import Version

if Version(dask.__version__) < Version("2024.8.0"):
    from dask.base import normalize_seq
else:
    from dask.tokenize import normalize_seq

from filelock import FileLock
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import subset_func  # noqa: F401

if TYPE_CHECKING:
    from collections.abc import Generator
    from types import EllipsisType


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"


@pytest.fixture(
    params=[("h5ad", None), ("zarr", 2), ("zarr", 3)], ids=["h5ad", "zarr2", "zarr3"]
)
def diskfmt(request):
    if (fmt := request.param[0]) == "h5ad":
        yield fmt
    else:
        with ad.settings.override(zarr_write_format=request.param[1]):
            yield fmt


@pytest.fixture
def diskfmt2(diskfmt):
    if diskfmt == "h5ad":
        with ad.settings.override(zarr_write_format=2):
            yield "zarr"
    else:
        yield "h5ad"


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


@pytest.fixture(scope="session")
def local_cluster_addr(
    tmp_path_factory: pytest.TempPathFactory, worker_id: str
) -> Generator[str, None, None]:
    # Adapted from https://pytest-xdist.readthedocs.io/en/latest/how-to.html#making-session-scoped-fixtures-execute-only-once
    import dask.distributed as dd

    def make_cluster() -> dd.LocalCluster:
        return dd.LocalCluster(n_workers=1, threads_per_worker=1)

    if worker_id == "master":
        with make_cluster() as cluster:
            yield cluster.scheduler_address
            return

    # get the temp directory shared by all workers
    root_tmp_dir = tmp_path_factory.getbasetemp().parent

    fn = root_tmp_dir / "dask_scheduler_address.txt"
    lock = FileLock(str(fn) + ".lock")
    lock.acquire()  # can’t use context manager, because we need to release the lock before yielding
    address = fn.read_text() if fn.is_file() else None
    if address:
        lock.release()
        yield address
        return

    with make_cluster() as cluster:
        fn.write_text(cluster.scheduler_address)
        lock.release()
        yield cluster.scheduler_address


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
