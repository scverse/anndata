"""
For tests using dask
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from zarr.storage import MemoryStore

import anndata as ad
from anndata._core.anndata import AnnData
from anndata.compat import CSArray, CSMatrix, CupyArray, DaskArray
from anndata.experimental.merge import as_group
from anndata.tests.helpers import (
    BASE_MATRIX_PARAMS,
    DASK_CAN_SPARRAY,
    DASK_MATRIX_PARAMS,
    GEN_ADATA_DASK_ARGS,
    as_dense_cupy_dask_array,
    as_dense_dask_array,
    as_sparse_dask_array,
    as_sparse_dask_matrix,
    assert_equal,
    gen_adata,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal

    from numpy.typing import NDArray
    from zarr.storage import MemoryStore


pytest.importorskip("dask.array")


@pytest.fixture(
    params=[
        [(2000, 1000), (100, 100)],
        [(200, 100), (100, 100)],
        [(200, 100), (100, 100)],
        [(20, 10), (1, 1)],
        [(20, 10), (1, 1)],
    ]
)
def sizes(request):
    return request.param


@pytest.fixture
def adata(sizes):
    import dask.array as da
    import numpy as np

    (M, N), chunks = sizes
    X = da.random.random((M, N), chunks=chunks)
    obs = pd.DataFrame(
        {"batch": np.random.choice(["a", "b"], M)},
        index=[f"cell{i:03d}" for i in range(M)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(N)])

    return AnnData(X, obs=obs, var=var)


def test_dask_X_view():
    import dask.array as da

    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.ones((M, N))
    view = adata[:30]
    view.copy()


def test_dask_write(adata, diskfmt_store, diskfmt):
    import dask.array as da
    import numpy as np

    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = getattr(ad, f"read_{diskfmt}")

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, diskfmt_store)
    curr = read(diskfmt_store)

    with pytest.raises(AssertionError):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


@pytest.mark.xdist_group("dask")
@pytest.mark.dask_distributed
def test_dask_distributed_write(
    adata: AnnData,
    tmp_path: Path,
    diskfmt: Literal["h5ad", "zarr"],
    local_cluster_addr: str,
) -> None:
    import dask.array as da
    import dask.distributed as dd
    import numpy as np

    assert isinstance(adata.X, DaskArray)
    # A real path, not an in-memory store: the dask workers are separate
    # processes and have to see the same store as this one.
    pth = tmp_path / f"test_write.{diskfmt}"
    with as_group(pth, mode="w") as g, dd.Client(local_cluster_addr):
        M, N = adata.X.shape
        adata.obsm["a"] = da.random.random((M, 10))
        adata.obsm["b"] = da.random.random((M, 10))
        adata.varm["a"] = da.random.random((N, 10))
        orig = adata
        ad.io.write_elem(g, "/", orig)
        # TODO: See https://github.com/zarr-developers/zarr-python/issues/2716
        with as_group(pth, mode="r") as g:
            curr = ad.io.read_elem(g)

    assert isinstance(curr, AnnData)
    with pytest.raises(AssertionError):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])
    assert_equal(curr.X, orig.X)

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_to_memory_check_array_types(adata, diskfmt_store, diskfmt):
    import dask.array as da
    import numpy as np

    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = getattr(ad, f"read_{diskfmt}")

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, diskfmt_store)
    curr = read(diskfmt_store)

    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)

    mem = orig.to_memory()

    with pytest.raises(AssertionError):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])
    assert_equal(mem.obsm["a"], orig.obsm["a"])
    assert_equal(mem.varm["a"], orig.varm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(mem.X, np.ndarray)
    assert isinstance(mem.obsm["a"], np.ndarray)
    assert isinstance(mem.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_to_memory_copy_check_array_types(adata, diskfmt_store, diskfmt):
    import dask.array as da
    import numpy as np

    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = getattr(ad, f"read_{diskfmt}")

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, diskfmt_store)
    curr = read(diskfmt_store)

    mem = orig.to_memory(copy=True)

    with pytest.raises(AssertionError):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])
    assert_equal(mem.obsm["a"], orig.obsm["a"])
    assert_equal(mem.varm["a"], orig.varm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(mem.X, np.ndarray)
    assert isinstance(mem.obsm["a"], np.ndarray)
    assert isinstance(mem.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_copy_check_array_types(adata):
    import dask.array as da

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    curr = adata.copy()

    with pytest.raises(AssertionError):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])

    assert isinstance(curr.X, DaskArray)
    assert isinstance(curr.obsm["a"], DaskArray)
    assert isinstance(curr.varm["a"], DaskArray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_assign_X(adata):
    """Check if assignment works"""
    import dask.array as da
    import numpy as np

    from anndata.compat import DaskArray

    adata.X = da.ones(adata.X.shape)
    prev_type = type(adata.X)
    adata_copy = adata.copy()

    adata.X = -1 * da.ones(adata.X.shape)
    assert prev_type is DaskArray
    assert type(adata_copy.X) is DaskArray
    assert_equal(adata.X, -1 * np.ones(adata.X.shape))
    assert_equal(adata_copy.X, np.ones(adata.X.shape))


# Test if dask arrays turn into numpy arrays after to_memory is called
@pytest.mark.parametrize(
    ("array_func", "mem_type"),
    [
        pytest.param(as_dense_dask_array, np.ndarray, id="dense"),
        pytest.param(as_sparse_dask_matrix, CSMatrix, id="sparse_matrix"),
        pytest.param(
            as_sparse_dask_array,
            CSArray,
            marks=pytest.mark.skipif(
                not DASK_CAN_SPARRAY, reason="Dask does not support sparrays"
            ),
            id="sparse_array",
        ),
        pytest.param(
            as_dense_cupy_dask_array, CupyArray, id="cupy_dense", marks=pytest.mark.gpu
        ),
    ],
)
def test_dask_to_memory_unbacked(array_func, mem_type):
    orig = gen_adata((15, 10), X_type=array_func, **GEN_ADATA_DASK_ARGS)
    orig.uns = {"da": {"da": array_func(np.ones((4, 12)))}}

    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)

    curr = orig.to_memory()

    assert_equal(orig, curr)
    assert isinstance(curr.X, mem_type)
    assert isinstance(curr.obsm["da"], np.ndarray)
    assert isinstance(curr.varm["da"], np.ndarray)
    assert isinstance(curr.layers["da"], np.ndarray)
    assert isinstance(curr.uns["da"]["da"], mem_type)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)


@pytest.mark.parametrize("to_dask", [*BASE_MATRIX_PARAMS, *DASK_MATRIX_PARAMS])
def test_dask_to_disk_view(
    to_dask: Callable[[NDArray], DaskArray],
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
) -> None:
    random_state = np.random.default_rng()
    arr = random_state.binomial(100, 0.005, (20, 15)).astype("float32")

    # TODO: need to change type for cupy
    orig = ad.AnnData(to_dask(arr))
    orig = orig[orig.shape[0] // 2]
    getattr(orig, f"write_{diskfmt}")(diskfmt_store)
    roundtrip = getattr(ad.io, f"read_{diskfmt}")(diskfmt_store)
    assert_equal(roundtrip, orig)


# Test if dask arrays turn into numpy arrays after to_memory is called
def test_dask_to_memory_copy_unbacked():
    import numpy as np

    orig = gen_adata((15, 10), X_type=as_dense_dask_array, **GEN_ADATA_DASK_ARGS)
    orig.uns = {"da": {"da": as_dense_dask_array(np.ones(12))}}

    curr = orig.to_memory(copy=True)

    assert_equal(orig, curr)
    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["da"], np.ndarray)
    assert isinstance(curr.varm["da"], np.ndarray)
    assert isinstance(curr.layers["da"], np.ndarray)
    assert isinstance(curr.uns["da"]["da"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)


def test_to_memory_raw():
    import dask.array as da
    import numpy as np

    orig = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)
    orig.X = da.ones((20, 10))

    with_raw = orig[:, ::2].copy()
    with_raw.raw = orig.copy()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)

    curr = with_raw.to_memory()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)
    assert isinstance(curr.raw.X, np.ndarray)
    assert isinstance(curr.raw.varm["da"], np.ndarray)


def test_to_memory_copy_raw():
    import dask.array as da
    import numpy as np

    orig = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)
    orig.X = da.ones((20, 10))

    with_raw = orig[:, ::2].copy()
    with_raw.raw = orig.copy()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)

    curr = with_raw.to_memory(copy=True)

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)
    assert isinstance(curr.raw.X, np.ndarray)
    assert isinstance(curr.raw.varm["da"], np.ndarray)


SMALL_SHARD_WARNING = r"dask chunk size is less than 100MB"
NO_SUBCHUNK_WARNING = r"Could not find inner chunks"


@pytest.mark.parametrize(
    (
        "shape",
        "dask_chunks",
        "diskfmt",
        "expected_chunks",
        "expected_shards",
        "expected_warnings",
    ),
    [
        pytest.param(
            # two >100MB shards, the last one partial
            (12000, 3000),
            (9000, 3000),
            "zarr",
            (500, 500),
            (9000, 3000),
            [],
            id="zarr-subchunked",
        ),
        pytest.param(
            (100, 100),
            (50, 100),
            "zarr",
            (50, 100),
            (50, 100),
            [SMALL_SHARD_WARNING],
            id="zarr-shard-under-1mb",
        ),
        pytest.param(
            # both axes prime, so the only divisors are 1 and the axis itself
            (1994, 1009),
            (997, 1009),
            "zarr",
            (997, 1009),
            None,
            [SMALL_SHARD_WARNING, NO_SUBCHUNK_WARNING],
            id="zarr-no-good-divisor",
        ),
        pytest.param(
            (2000, 2000), (1000, 2000), "h5ad", (1000, 2000), None, [], id="h5ad"
        ),
    ],
)
def test_dask_write_chunks_and_shards(
    tmp_path: Path,
    shape: tuple[int, int],
    dask_chunks: tuple[int, int],
    diskfmt: Literal["h5ad", "zarr"],
    expected_chunks: tuple[int, int],
    expected_shards: tuple[int, int] | None,
    expected_warnings: list[str],
) -> None:
    import dask.array as da

    X = da.zeros(shape, chunks=dask_chunks, dtype=np.float32)
    store = tmp_path / "test.h5ad" if diskfmt == "h5ad" else MemoryStore()
    with as_group(store, mode="w") as g:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            ad.io.write_elem(g, "X", X)
        messages = [str(w.message) for w in record]
        assert len(messages) == len(expected_warnings), messages
        for pattern, message in zip(expected_warnings, messages, strict=True):
            assert re.search(pattern, message), message
        arr = g["X"]
        assert isinstance(arr, zarr.Array | h5py.Dataset)
        assert arr.chunks == expected_chunks
        assert getattr(arr, "shards", None) == expected_shards
