from __future__ import annotations

import re

import h5py
import zarr
from scipy import sparse

import anndata as ad
from anndata.compat import SpArray
from anndata.experimental import (
    read_dispatched,
    read_elem,
    write_dispatched,
    write_elem,
)
from anndata.tests.helpers import assert_equal, gen_adata


def test_read_dispatched_w_regex():
    def read_only_axis_dfs(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata":
            return func(elem)
        elif re.match(r"^/((obs)|(var))?(/.*)?$", elem_name):
            return func(elem)
        else:
            return None

    adata = gen_adata((1000, 100))
    z = zarr.group()

    write_elem(z, "/", adata)

    expected = ad.AnnData(obs=adata.obs, var=adata.var)
    actual = read_dispatched(z, read_only_axis_dfs)

    assert_equal(expected, actual)


def test_read_dispatched_dask():
    import dask.array as da

    def read_as_dask_array(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in {
            "dataframe",
            "csr_matrix",
            "csc_matrix",
            "awkward-array",
        }:
            # Preventing recursing inside of these types
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        else:
            return func(elem)

    adata = gen_adata((1000, 100))
    z = zarr.group()
    write_elem(z, "/", adata)

    dask_adata = read_dispatched(z, read_as_dask_array)

    assert isinstance(dask_adata.layers["array"], da.Array)
    assert isinstance(dask_adata.obsm["array"], da.Array)
    assert isinstance(dask_adata.uns["nested"]["nested_further"]["array"], da.Array)

    expected = read_elem(z)
    actual = dask_adata.to_memory(copy=False)

    assert_equal(expected, actual)


def test_read_dispatched_null_case():
    adata = gen_adata((100, 100))
    z = zarr.group()
    write_elem(z, "/", adata)

    expected = read_elem(z)

    def callback(read_func, elem_name, x, iospec):
        return read_elem(x)

    actual = read_dispatched(z, callback)

    assert_equal(expected, actual)


def test_write_dispatched_chunks():
    from itertools import chain, repeat

    def determine_chunks(elem_shape, specified_chunks):
        chunk_iterator = chain(specified_chunks, repeat(None))
        return tuple(e if c is None else c for e, c in zip(elem_shape, chunk_iterator))

    adata = gen_adata((1000, 100))

    def write_chunked(func, store, k, elem, dataset_kwargs, iospec):
        M, N = 13, 42

        def set_copy(d, **kwargs):
            d = dict(d)
            d.update(kwargs)
            return d

        # TODO: Should the passed path be absolute?
        path = "/" + store.path + "/" + k
        if hasattr(elem, "shape") and not isinstance(
            elem, (sparse.spmatrix, SpArray, ad.AnnData)
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

    z = zarr.group()

    write_dispatched(z, "/", adata, callback=write_chunked)

    def check_chunking(k, v):
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

    z.visititems(check_chunking)


def test_io_dispatched_keys(tmp_path):
    h5ad_write_keys = []
    zarr_write_keys = []
    h5ad_read_keys = []
    zarr_read_keys = []

    h5ad_path = tmp_path / "test.h5ad"
    zarr_path = tmp_path / "test.zarr"

    def h5ad_writer(func, store, k, elem, dataset_kwargs, iospec):
        h5ad_write_keys.append(k)
        func(store, k, elem, dataset_kwargs=dataset_kwargs)

    def zarr_writer(func, store, k, elem, dataset_kwargs, iospec):
        zarr_write_keys.append(k)
        func(store, k, elem, dataset_kwargs=dataset_kwargs)

    def h5ad_reader(func, elem_name: str, elem, iospec):
        h5ad_read_keys.append(elem_name)
        return func(elem)

    def zarr_reader(func, elem_name: str, elem, iospec):
        zarr_read_keys.append(elem_name)
        return func(elem)

    adata = gen_adata((50, 100))

    with h5py.File(h5ad_path, "w") as f:
        write_dispatched(f, "/", adata, callback=h5ad_writer)
        _ = read_dispatched(f, h5ad_reader)

    with zarr.open_group(zarr_path, "w") as f:
        write_dispatched(f, "/", adata, callback=zarr_writer)
        _ = read_dispatched(f, zarr_reader)

    assert h5ad_write_keys == zarr_write_keys
    assert h5ad_read_keys == zarr_read_keys

    assert sorted(h5ad_write_keys) == sorted(h5ad_read_keys)
