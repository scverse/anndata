"""Tests for the attribute .X"""
import numpy as np
from scipy import sparse

from anndata.utils import asarray

import pytest

from anndata.tests.helpers import gen_adata

UNLABELLED_ARRAY_TYPES = [
    pytest.param(sparse.csr_matrix, id="csr"),
    pytest.param(sparse.csc_matrix, id="csc"),
    pytest.param(asarray, id="ndarray"),
]
SINGULAR_SHAPES = [
    pytest.param(shape, id=str(shape)) for shape in [(1, 10), (10, 1), (1, 1)]
]


@pytest.mark.parametrize("shape", SINGULAR_SHAPES)
@pytest.mark.parametrize("orig_array_type", UNLABELLED_ARRAY_TYPES)
@pytest.mark.parametrize("new_array_type", UNLABELLED_ARRAY_TYPES)
def test_setter_singular_dim(shape, orig_array_type, new_array_type):
    adata = gen_adata(shape, X_type=orig_array_type)
    adata.X = new_array_type(np.ones(shape))
    np.testing.assert_equal(asarray(adata.X), 1)
