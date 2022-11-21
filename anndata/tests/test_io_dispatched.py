import re

from scipy import sparse
import zarr

import anndata as ad
from anndata.experimental import (
    read_dispatched,
    write_dispatched,
    read_elem,
    write_elem,
)
from anndata.tests.helpers import gen_adata, assert_equal


def test_read_dispatched_w_regex():
    def read_only_axis_dfs(func, elem_name: str, elem, iospec):
        if elem_name == "/":
            # TODO: this case is only complicated because of AnnData.__init__ dtype arg
            d = {}
            for k, v in elem.items():
                v_read = read_dispatched(v, read_only_axis_dfs)
                if v_read is not None:
                    d[k] = v_read
            return ad.AnnData(**d)
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
        if iospec.encoding_type in ("dataframe", "csr_matrix", "csc_matrix"):
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
    actual = read_dispatched(z, lambda _, __, x, ___: read_elem(x))

    assert_equal(expected, actual)


def test_write_dispatched_chunks():
    from itertools import repeat, chain

    def determine_chunks(elem_shape, specified_chunks):
        chunk_iterator = chain(specified_chunks, repeat(None))
        return tuple(e if c is None else c for e, c in zip(elem_shape, chunk_iterator))

    adata = gen_adata((1000, 100))

    def write_chunked(func, store, k, elem, dataset_kwargs):
        M, N = 13, 42

        def set_copy(d, **kwargs):
            d = dict(d)
            d.update(kwargs)
            return d

        # TODO: Should the passed path be absolute?
        path = "/" + store.path + "/" + k
        if hasattr(elem, "shape") and not isinstance(
            elem, (sparse.spmatrix, ad.AnnData)
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
            if not isinstance(v, zarr.Array) or v.shape == ():
                return
            if re.match(r"obs[mp]?/\w+", k):
                assert v.chunks[0] == 13
            elif re.match(r"var[mp]?/\w+", k):
                assert v.chunks[0] == 42

        z.visititems(check_chunking)
