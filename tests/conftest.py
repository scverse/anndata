from __future__ import annotations

import warnings
from functools import partial

import joblib
import pytest
from dask.base import normalize_seq, normalize_token, tokenize
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import subset_func  # noqa: F401

# TODO: Should be done in pyproject.toml eventually
# See https://github.com/pytest-dev/pytest-cov/issues/437
warnings.filterwarnings("ignore", category=ad.OldFormatWarning)
warnings.filterwarnings("ignore", category=ad.ImplicitModificationWarning)


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"


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
        res.append(tokenize(list(elem.items())))
    res.append(joblib.hash(adata.uns))
    if adata.raw is not None:
        res.append(tokenize(adata.raw.to_adata()))
    return tuple(res)
