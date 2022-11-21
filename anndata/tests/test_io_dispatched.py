import re

import anndata as ad
from anndata.experimental import read_dispatched, read_elem, write_elem
from anndata.tests.helpers import gen_adata, assert_equal

import zarr


def test_read_dispatched_w_regex():
    def read_only_axis_dfs(func, elem_name: str, elem, iospec):
        if elem_name == "/":
            # TODO: this case is only needed due to dtype
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


def test_read_dispatched_null_case():
    adata = gen_adata((100, 100))
    z = zarr.group()
    write_elem(z, "/", adata)

    expected = read_elem(z)
    actual = read_dispatched(z, lambda _, __, x, ___: read_elem(x))

    assert_equal(expected, actual)
