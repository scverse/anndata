from __future__ import annotations

import asyncio
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
import zarr.api.asynchronous as za
from scipy import sparse
from zarr.core.sync import sync
from zarr.storage import LocalStore, WrapperStore

import anndata as ad
from anndata._io.specs.registry import IORegistryError
from anndata.compat import awkward as ak
from anndata.compat import zarr_loop
from anndata.experimental import read_elem_async
from anndata.io import read_elem
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata

if TYPE_CHECKING:
    from collections.abc import Callable, Coroutine
    from pathlib import Path
    from typing import Any

pytestmark = pytest.mark.zarr_io


@pytest.mark.parametrize(
    "run",
    [
        pytest.param(asyncio.run, id="other_loop"),
        pytest.param(sync, id="zarr_loop"),
    ],
)
def test_many_stores_concurrently(
    tmp_path: Path, run: Callable[[Coroutine[Any, Any, list[Any]]], list[Any]]
) -> None:
    groups = []
    for i in range(8):
        adata = gen_adata((20 + i, 10), **GEN_ADATA_NO_XARRAY_ARGS)
        adata.write_zarr(tmp_path / f"{i}.zarr")
        groups.append(zarr.open_group(tmp_path / f"{i}.zarr", mode="r"))

    # Sync zarr calls like `g["obs"]` fail on zarr’s loop, so look up elements outside it
    elems = groups + [g["obs"] for g in groups]

    async def read_all():
        return await asyncio.gather(*map(read_elem_async, elems))

    results = run(read_all())
    for elem, from_async in zip(elems, results, strict=True):
        assert_equal(read_elem(elem), from_async)


def test_rejects_h5py(tmp_path: Path) -> None:
    with h5py.File(tmp_path / "x.h5", "w") as f:
        ad.io.write_elem(f, "a", np.arange(3))
        with pytest.raises(TypeError, match=r"reads zarr groups and arrays"):
            asyncio.run(read_elem_async(f["a"]))


def test_error_reports_key(tmp_path: Path) -> None:
    g = zarr.open_group(tmp_path / "x.zarr", mode="w")
    ad.io.write_elem(g, "df", ad.AnnData(np.ones((3, 2))).obs)
    g["df"].attrs["column-order"] = ["missing"]
    with pytest.raises(KeyError) as e:
        asyncio.run(read_elem_async(g["df"]))
    assert any("'df'" in note for note in e.value.__notes__)


def test_loop_affinity(tmp_path: Path) -> None:
    """Sync zarr objects are read on zarr’s loop, async ones on the running loop."""
    loops: set[asyncio.AbstractEventLoop] = set()

    class LoopRecordingStore(WrapperStore):
        async def get(self, *args, **kwargs):
            loops.add(asyncio.get_running_loop())
            return await super().get(*args, **kwargs)

    df = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
    ad.io.write_elem(zarr.open_group(tmp_path / "x.zarr", mode="w"), "df", df)
    store = LoopRecordingStore(LocalStore(tmp_path / "x.zarr", read_only=True))

    elem = zarr.open_group(store, mode="r")["df"]
    loops.clear()
    assert_equal(df, asyncio.run(read_elem_async(elem)))
    assert loops == {zarr_loop()}

    async def read_on_this_loop() -> asyncio.AbstractEventLoop:
        elem = await (await za.open_group(store, mode="r")).getitem("df")
        loops.clear()
        assert_equal(df, await read_elem_async(elem))
        return asyncio.get_running_loop()

    assert loops == {asyncio.run(read_on_this_loop())}


def test_anndata_ignores_unrelated_children(tmp_path: Path) -> None:
    """Like the sync reader, only look up AnnData’s own keys, without listing the group."""
    path = tmp_path / "x.zarr"
    ad.AnnData(np.ones((3, 2))).write_zarr(path)
    (path / "junk").mkdir()
    (path / "junk" / "zarr.json").write_text("not json")
    g = zarr.open_group(path, mode="r", use_consolidated=False)
    assert_equal(read_elem(g), asyncio.run(read_elem_async(g)))


def test_sparse_legacy_attrs(tmp_path: Path) -> None:
    """Like the sync reader, honor `h5sparse_shape` and `h5sparse_format`."""
    g = zarr.open_group(tmp_path / "x.zarr", mode="w")
    ad.io.write_elem(g, "m", sparse.random(5, 4, density=0.5, format="csc"))
    attrs = dict(g["m"].attrs)
    attrs["h5sparse_shape"] = attrs.pop("shape")
    attrs |= {"encoding-type": "csr_matrix", "h5sparse_format": "csc"}
    g["m"].attrs.clear()
    g["m"].attrs.update(attrs)
    from_async = asyncio.run(read_elem_async(g["m"]))
    assert isinstance(from_async, sparse.csc_matrix)
    assert_equal(read_elem(g["m"]), from_async)


def test_recarray_and_awkward(tmp_path: Path) -> None:
    adata = gen_adata((5, 4), **GEN_ADATA_NO_XARRAY_ARGS)
    adata.uns["rec"] = np.rec.fromarrays(
        [np.array(["a", "b"]), np.arange(2)], names=["name", "score"]
    )
    adata.obsm["awk"] = ak.Array([[1], [2, 3], [], [4], [5, 6]])
    adata.write_zarr(tmp_path / "x.zarr")
    g = zarr.open_group(tmp_path / "x.zarr", mode="r")
    assert_equal(read_elem(g), asyncio.run(read_elem_async(g)))


def test_no_async_reader(tmp_path: Path) -> None:
    """An element that only `read_elem` reads is an error, with the key in a note."""
    g = zarr.open_group(tmp_path / "x.zarr", mode="w")
    ad.io.write_elem(g, "uns", {"a": 1})
    uns = g.require_group("uns")
    uns.create_group("plain")
    from_sync: Any = read_elem(uns)
    assert from_sync["plain"] == {}
    with pytest.raises(NotImplementedError, match=r"No async read method") as e:
        asyncio.run(read_elem_async(uns))
    assert any("'plain'" in note for note in e.value.__notes__)

    a: Any = uns["a"]
    a.attrs["encoding-type"] = "unknown"
    with pytest.raises(IORegistryError, match=r"No read method registered"):
        asyncio.run(read_elem_async(a))


def test_error_cancels_other_reads(tmp_path: Path) -> None:
    """When one child fails, the reads of its siblings are cancelled."""
    ad.AnnData(np.ones((3, 2))).write_zarr(tmp_path / "x.zarr")
    obs = zarr.open_group(tmp_path / "x.zarr", mode="r+").require_group("obs")
    obs.attrs["column-order"] = ["missing"]
    cancelled: list[str] = []

    class Store(WrapperStore):
        async def get(self, key: str, *args, **kwargs):
            if key.startswith("X/"):
                # Stay pending until the error in `obs` cancels this read
                try:
                    await asyncio.sleep(60)
                except asyncio.CancelledError:
                    cancelled.append(key)
                    raise
            elif key.startswith("obs/missing/"):
                # Let the read of `X` start first
                await asyncio.sleep(0.2)
            return await super().get(key, *args, **kwargs)

    store = Store(LocalStore(tmp_path / "x.zarr", read_only=True))
    g = zarr.open_group(store, mode="r", use_consolidated=False)
    with pytest.raises(KeyError):
        asyncio.run(read_elem_async(g))
    assert cancelled
