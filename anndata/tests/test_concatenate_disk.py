from collections.abc import Hashable
from copy import deepcopy
from itertools import chain, product
from functools import partial, singledispatch
from typing import Any, List, Callable, Mapping
import warnings

import numpy as np
from numpy import ma
import pandas as pd
from pandas.api.types import is_categorical_dtype
import pytest
from scipy import sparse
from boltons.iterutils import research, remap, default_exit

from anndata._io.merge import concat_on_disk
from anndata import AnnData, Raw, concat
from anndata._core.index import _subset
from anndata._core import merge
from anndata.tests import helpers
from anndata.tests.helpers import (
    assert_equal,
    gen_adata,
)
from anndata.utils import asarray
from anndata.compat import DaskArray, AwkArray


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
            adata.write_h5ad(path)
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
            p = tmp_path / (f"{i}."+file_format)
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
    out_name = tmp_path / ("out."+file_format)
    concat_on_disk(paths, out_name, *args, **kwargs)
    res2 = read_func(out_name)
    assert_equal(res1, res2)


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

        a = gen_adata((M, N), X_type=get_array_type(
            array_type, axis), sparse_fmt=sparse_fmt,
            **GEN_ADATA_OOC_CONCAT_ARGS)
        adatas.append(a)

    assert_eq_concat_on_disk(
        adatas, tmp_path, file_format, axis=axis, join=join_type)
