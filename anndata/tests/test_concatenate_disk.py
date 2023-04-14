from typing import Mapping

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata._io.merge import concat_on_disk
from anndata import AnnData, concat
from anndata.tests.helpers import (
    assert_equal,
    gen_adata,
)


from anndata.utils import asarray


from anndata import read_h5ad, read_zarr


GEN_ADATA_OOC_CONCAT_ARGS = dict(
    obsm_types=(
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
    ),
    varm_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame),
    layers_types=(sparse.spmatrix, np.ndarray, pd.DataFrame),
)


@pytest.fixture(params=[0, 1])
def axis(request):
    return request.param


@pytest.fixture(
    params=["array", "sparse"],
)
def array_type(request):
    return request.param


@pytest.fixture(params=["inner"])
def join_type(request):
    return request.param


@pytest.fixture(params=["zarr", "h5ad"])
def file_format(request):
    return request.param


def _adatas_to_paths(adatas, tmp_path, file_format):
    """
    Gets list of adatas, writes them and returns their paths as zarr
    """
    paths = None

    def write_func(adata, path):
        if file_format == "h5ad":
            adata.write(path)
        else:
            adata.write_zarr(path)

    if isinstance(adatas, Mapping):
        paths = {}
        for k, v in adatas.items():
            p = tmp_path / (f"{k}." + file_format)
            write_func(v, p)
            paths[k] = p
    else:
        paths = []
        for i, a in enumerate(adatas):
            p = tmp_path / (f"{i}." + file_format)
            write_func(a, p)
            paths += [p]
    return paths


def assert_eq_concat_on_disk(adatas, tmp_path, file_format, *args, **kwargs):
    def read_func(path):
        if file_format == "h5ad":
            return read_h5ad(path)
        return read_zarr(path)

    # create one from the concat function
    res1 = concat(adatas, *args, **kwargs)
    # create one from the on disk concat function
    paths = _adatas_to_paths(adatas, tmp_path, file_format)
    out_name = tmp_path / ("out." + file_format)
    concat_on_disk(paths, out_name, *args, **kwargs)
    res2 = read_func(out_name)
    assert_equal(res1, res2, exact=False)


def get_array_type(array_type, axis):
    if array_type == "sparse":
        if axis == 0:
            return sparse.csr_matrix
        return sparse.csc_matrix
    if array_type == "array":
        return asarray
    else:
        raise NotImplementedError(f"array_type {array_type} not implemented")


def test_anndatas_without_reindex(axis, array_type, join_type, tmp_path, file_format):
    N = 50
    M = 50
    sparse_fmt = "csr"
    adatas = []
    for _ in range(5):
        if axis == 0:
            M = np.random.randint(1, 100)
        else:
            N = np.random.randint(1, 100)
            sparse_fmt = "csc"

        a = gen_adata(
            (M, N),
            X_type=get_array_type(array_type, axis),
            sparse_fmt=sparse_fmt,
            **GEN_ADATA_OOC_CONCAT_ARGS,
        )
        adatas.append(a)

    assert_eq_concat_on_disk(adatas, tmp_path, file_format, axis=axis, join=join_type)


def test_concat_ordered_categoricals_retained(tmp_path, file_format):
    a = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("aabcd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5)],
        ),
    )
    b = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("abcdd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5, 10)],
        ),
    )

    adatas = [a, b]
    assert_eq_concat_on_disk(adatas, tmp_path, file_format)


@pytest.fixture
def obsm_adatas():
    def gen_index(n):
        return [f"cell{i}" for i in range(n)]

    return [
        AnnData(
            X=sparse.csr_matrix((3, 5)),
            obs=pd.DataFrame(index=gen_index(3)),
            obsm={
                "dense": np.arange(6).reshape(3, 2),
                "sparse": sparse.csr_matrix(np.arange(6).reshape(3, 2)),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(3),
                        "b": list("abc"),
                        "c": pd.Categorical(list("aab")),
                    },
                    index=gen_index(3),
                ),
            },
        ),
        AnnData(
            X=sparse.csr_matrix((4, 10)),
            obs=pd.DataFrame(index=gen_index(4)),
            obsm=dict(
                dense=np.arange(12).reshape(4, 3),
                df=pd.DataFrame(dict(a=np.arange(3, 7)), index=gen_index(4)),
            ),
        ),
        AnnData(
            X=sparse.csr_matrix((2, 100)),
            obs=pd.DataFrame(index=gen_index(2)),
            obsm={
                "sparse": np.arange(8).reshape(2, 4),
                "dense": np.arange(4, 8).reshape(2, 2),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(7, 9),
                        "b": list("cd"),
                        "c": pd.Categorical(list("ab")),
                    },
                    index=gen_index(2),
                ),
            },
        ),
    ]


def test_concatenate_obsm_inner(obsm_adatas, tmp_path, file_format):
    assert_eq_concat_on_disk(obsm_adatas, tmp_path, file_format, join="inner")


@pytest.mark.parametrize("elem", ["sparse", "array", "df"])
def test_concat_outer_aligned_mapping(elem, tmp_path, file_format):
    a = gen_adata((5, 5), **GEN_ADATA_OOC_CONCAT_ARGS)
    b = gen_adata((3, 5), **GEN_ADATA_OOC_CONCAT_ARGS)

    del b.obsm[elem]
    if elem == "df":
        del a.obsm[elem]["bool"]

    adatas = concat({"a": a, "b": b}, join="outer", label="group")

    assert_eq_concat_on_disk(adatas, tmp_path, file_format)
