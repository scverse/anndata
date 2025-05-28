from __future__ import annotations

import re
from typing import TYPE_CHECKING

import h5py
import zarr

import anndata as ad
from anndata._io.zarr import open_write_group
from anndata.compat import CSArray, CSMatrix, ZarrGroup, is_zarr_v2
from anndata.experimental import read_dispatched, write_dispatched
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path


def test_read_dispatched_w_regex(tmp_path: Path):
    def read_only_axis_dfs(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or re.match(
            r"^/((obs)|(var))?(/.*)?$", elem_name
        ):
            return func(elem)
        else:
            return None

    adata = gen_adata((1000, 100), **GEN_ADATA_NO_XARRAY_ARGS)
    z = open_write_group(tmp_path)

    ad.io.write_elem(z, "/", adata)
    # TODO: see https://github.com/zarr-developers/zarr-python/issues/2716
    if not is_zarr_v2() and isinstance(z, ZarrGroup):
        z = zarr.open(z.store)

    expected = ad.AnnData(obs=adata.obs, var=adata.var)
    actual = read_dispatched(z, read_only_axis_dfs)

    assert_equal(expected, actual)


def test_read_dispatched_dask(tmp_path: Path):
    import dask.array as da

    def read_as_dask_array(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in {
            "dataframe",
            "csr_matrix",
            "csc_matrix",
            "awkward-array",
        }:
            # Preventing recursing inside of these types
            return ad.io.read_elem(elem)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        else:
            return func(elem)

    adata = gen_adata((1000, 100), **GEN_ADATA_NO_XARRAY_ARGS)
    z = open_write_group(tmp_path)
    ad.io.write_elem(z, "/", adata)
    # TODO: see https://github.com/zarr-developers/zarr-python/issues/2716
    if not is_zarr_v2() and isinstance(z, ZarrGroup):
        z = zarr.open(z.store)

    dask_adata = read_dispatched(z, read_as_dask_array)

    assert isinstance(dask_adata.layers["array"], da.Array)
    assert isinstance(dask_adata.obsm["array"], da.Array)
    assert isinstance(dask_adata.uns["nested"]["nested_further"]["array"], da.Array)

    expected = ad.io.read_elem(z)
    actual = dask_adata.to_memory(copy=False)

    assert_equal(expected, actual)


def test_read_dispatched_null_case(tmp_path: Path):
    adata = gen_adata((100, 100), **GEN_ADATA_NO_XARRAY_ARGS)
    z = open_write_group(tmp_path)
    ad.io.write_elem(z, "/", adata)
    # TODO: see https://github.com/zarr-developers/zarr-python/issues/2716
    if not is_zarr_v2() and isinstance(z, ZarrGroup):
        z = zarr.open(z.store)
    expected = ad.io.read_elem(z)
    actual = read_dispatched(z, lambda _, __, x, **___: ad.io.read_elem(x))

    assert_equal(expected, actual)


def test_write_dispatched_chunks(tmp_path: Path):
    from itertools import chain, repeat

    def determine_chunks(elem_shape, specified_chunks):
        chunk_iterator = chain(specified_chunks, repeat(None))
        return tuple(
            e if c is None else c
            for e, c in zip(elem_shape, chunk_iterator, strict=False)
        )

    adata = gen_adata((1000, 100), **GEN_ADATA_NO_XARRAY_ARGS)

    def write_chunked(func, store, k, elem, dataset_kwargs, iospec):
        M, N = 13, 42

        def set_copy(d, **kwargs):
            d = dict(d)
            d.update(kwargs)
            return d

        # TODO: Should the passed path be absolute?
        path = "/" + store.path + "/" + k
        if hasattr(elem, "shape") and not isinstance(
            elem, CSMatrix | CSArray | ad.AnnData
        ):
            if re.match(r"^/((X)|(layers)).*", path):
                chunks = (M, N)
            elif path.startswith("/obsp"):
                chunks = (M, M)
            elif path.startswith("/obs"):
                chunks = (M,)
            elif path.startswith("/varp"):
                chunks = (N, N)
            elif path.startswith("/var"):
                chunks = (N,)
            else:
                chunks = dataset_kwargs.get("chunks", ())
            func(
                store,
                k,
                elem,
                dataset_kwargs=set_copy(
                    dataset_kwargs, chunks=determine_chunks(elem.shape, chunks)
                ),
            )
        else:
            func(store, k, elem, dataset_kwargs=dataset_kwargs)

    z = open_write_group(tmp_path)

    write_dispatched(z, "/", adata, callback=write_chunked)

    def check_chunking(k: str, v: ZarrGroup | zarr.Array):
        if (
            not isinstance(v, zarr.Array)
            or v.shape == ()
            or any(k.endswith(x) for x in ("data", "indices", "indptr"))
        ):
            return
        if re.match(r"obs[mp]?/\w+", k):
            assert v.chunks[0] == 13
        elif re.match(r"var[mp]?/\w+", k):
            assert v.chunks[0] == 42

    if is_zarr_v2():
        z.visititems(check_chunking)
    else:

        def visititems(
            z: ZarrGroup, visitor: Callable[[str, ZarrGroup | zarr.Array], None]
        ) -> None:
            for key in z:
                maybe_group = z[key]
                if isinstance(maybe_group, ZarrGroup):
                    visititems(maybe_group, visitor)
                else:
                    visitor(key, maybe_group)

        visititems(z, check_chunking)


def test_io_dispatched_keys(tmp_path: Path):
    h5ad_write_keys = []
    zarr_write_keys = []
    h5ad_read_keys = []
    zarr_read_keys = []

    h5ad_path = tmp_path / "test.h5ad"
    zarr_path = tmp_path / "test.zarr"

    def h5ad_writer(func, store, k, elem, dataset_kwargs, iospec):
        h5ad_write_keys.append(k if is_zarr_v2() else k.strip("/"))
        func(store, k, elem, dataset_kwargs=dataset_kwargs)

    def zarr_writer(func, store, k, elem, dataset_kwargs, iospec):
        zarr_write_keys.append(
            k if is_zarr_v2() else f"{store.name.strip('/')}/{k.strip('/')}".strip("/")
        )
        func(store, k, elem, dataset_kwargs=dataset_kwargs)

    def h5ad_reader(func, elem_name: str, elem, iospec):
        h5ad_read_keys.append(elem_name if is_zarr_v2() else elem_name.strip("/"))
        return func(elem)

    def zarr_reader(func, elem_name: str, elem, iospec):
        zarr_read_keys.append(elem_name if is_zarr_v2() else elem_name.strip("/"))
        return func(elem)

    adata = gen_adata((50, 100), **GEN_ADATA_NO_XARRAY_ARGS)

    with h5py.File(h5ad_path, "w") as f:
        write_dispatched(f, "/", adata, callback=h5ad_writer)
        _ = read_dispatched(f, h5ad_reader)

    f = open_write_group(zarr_path)
    write_dispatched(f, "/", adata, callback=zarr_writer)
    _ = read_dispatched(f, zarr_reader)

    assert sorted(h5ad_read_keys) == sorted(zarr_read_keys)
    assert sorted(h5ad_write_keys) == sorted(zarr_write_keys)
    for sub_sparse_key in ["data", "indices", "indptr"]:
        assert f"/X/{sub_sparse_key}" not in h5ad_read_keys
        assert f"/X/{sub_sparse_key}" not in h5ad_write_keys
