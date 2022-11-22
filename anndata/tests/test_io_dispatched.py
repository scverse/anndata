import re
from anndata._io.specs.registry import IOSpec
from anndata.compat import ZarrArray
from scipy import sparse
import zarr

import anndata as ad
from anndata import read_dispatched, write_dispatched
from anndata.tests.helpers import gen_adata, assert_equal


def test_read_dispatched_w_regex():
    def read_only_axis_dfs(func, group, elem_name, iospec):
        if elem_name in set(["var", "obs"]):
            return func(group[elem_name])
        else:
            return None

    adata = gen_adata((1000, 100))
    adata.write_zarr("test.zarr")
    expected = ad.AnnData(obs=adata.obs, var=adata.var)
    z = zarr.open("test.zarr")
    actual = read_dispatched(z, dispatch_element=read_only_axis_dfs)

    assert_equal(expected, actual)


def test_read_dispatched_dask():
    import dask.array as da

    def register_element_dispatchers(register):
        register(ZarrArray, IOSpec("array", "0.2.0"))(da.from_zarr)

    adata = gen_adata((1000, 100))
    store = zarr.DirectoryStore("test.zarr")
    z = zarr.group(store=store, overwrite=True)
    adata.write_zarr(store)

    dask_adata = read_dispatched(
        z, register_element_dispatchers=register_element_dispatchers
    )

    assert isinstance(dask_adata.layers["array"], da.Array)
    assert isinstance(dask_adata.obsm["array"], da.Array)
    assert isinstance(dask_adata.uns["nested"]["nested_further"]["array"], da.Array)

    actual = dask_adata.to_memory(copy=False)
    expected = adata

    assert_equal(expected, actual)


def test_read_dispatched_null_case():
    adata = gen_adata((100, 100))
    store = zarr.DirectoryStore("test.zarr")
    z = zarr.group(store=store, overwrite=True)
    adata.write_zarr(store)
    expected = adata
    actual = read_dispatched(z)

    assert_equal(expected, actual)


def test_write_dispatched_chunks():
    from itertools import repeat, chain

    def determine_chunks(elem_shape, specified_chunks):
        chunk_iterator = chain(specified_chunks, repeat(None))
        return tuple(e if c is None else c for e, c in zip(elem_shape, chunk_iterator))

    adata = gen_adata((1000, 100))

    dataset_kwargs = {}

    def write_chunked(func, group, key, elem):
        M, N = 13, 42

        def set_copy(d, **kwargs):
            d = dict(d)
            d.update(kwargs)
            return d

        if hasattr(elem, "shape") and not isinstance(
            elem, (sparse.spmatrix, ad.AnnData)
        ):
            if re.match(r"^/((X)|(layers)).*", key):
                chunks = (M, N)
            elif key.startswith("obsp"):
                chunks = (M, M)
            elif key.startswith("obs"):
                chunks = (M,)
            elif key.startswith("varp"):
                chunks = (N, N)
            elif key.startswith("var"):
                chunks = (N,)
            else:
                chunks = dataset_kwargs.get("chunks", ())
            func(
                group,
                key,
                elem,
                dataset_kwargs=set_copy(
                    dataset_kwargs, chunks=determine_chunks(elem.shape, chunks)
                ),
            )
        else:
            func(group, key, elem, dataset_kwargs=dataset_kwargs)

        store = zarr.DirectoryStore("test.zarr")
        z = zarr.group(store=store, overwrite=True)

        write_dispatched(z, adata, dispatch_element=write_chunked)

        def check_chunking(k, v):
            if not isinstance(v, zarr.Array) or v.shape == ():
                return
            if re.match(r"obs[mp]?/\w+", k):
                assert v.chunks[0] == 13
            elif re.match(r"var[mp]?/\w+", k):
                assert v.chunks[0] == 42

        z.visititems(check_chunking)
