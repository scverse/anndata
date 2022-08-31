"""Tests related to awkward arrays"""
import pytest
import numpy as np
import numpy.testing as npt

from anndata.tests.helpers import gen_adata, gen_awkward
from anndata.compat import awkward as ak
from anndata import ImplicitModificationWarning


@pytest.mark.parametrize(
    "field,value,valid",
    [
        ["obsm", gen_awkward((10, 5)), True],
        ["obsm", gen_awkward((10, None)), True],
        ["obsm", gen_awkward((10, None, None)), True],
        ["obsm", gen_awkward((10, 5, None)), True],
        ["obsm", gen_awkward((8, 10)), False],
        ["obsm", gen_awkward((8, None)), False],
        ["varm", gen_awkward((20, 5)), True],
        ["varm", gen_awkward((20, None)), True],
        ["varm", gen_awkward((20, None, None)), True],
        ["varm", gen_awkward((20, 5, None)), True],
        ["varm", gen_awkward((8, 20)), False],
        ["varm", gen_awkward((8, None)), False],
        ["uns", gen_awkward((7,)), True],
        ["uns", gen_awkward((7, None)), True],
        ["uns", gen_awkward((7, None, None)), True],
    ],
)
def test_set_awkward(field, value, valid):
    """Check if we can set .X, .layers, .obsm, .varm and .uns with different types
    of awkward arrays and if error messages are properly raised when the dimensions do not align.
    """
    adata = gen_adata((10, 20), varm_types=(), obsm_types=(), layers_types=())

    def _assign():
        getattr(adata, field)["test"] = value

    if not valid:
        with pytest.raises(ValueError):
            _assign()
    else:
        _assign()


@pytest.mark.parametrize("key", ["obsm", "varm", "uns"])
def test_copy(key):
    """Check that modifying a copy does not modify the original"""
    adata = gen_adata((3, 3), varm_types=(), obsm_types=(), layers_types=())
    getattr(adata, key)["awk"] = ak.Array([{"a": [1], "b": [2], "c": [3]}] * 3)
    adata_copy = adata.copy()
    getattr(adata_copy, key)["awk"]["c"] = np.full((3, 1), 4)
    getattr(adata_copy, key)["awk"]["d"] = np.full((3, 1), 5)

    # values in copy were correctly set
    npt.assert_equal(getattr(adata_copy, key)["awk"]["c"], np.full((3, 1), 4))
    npt.assert_equal(getattr(adata_copy, key)["awk"]["d"], np.full((3, 1), 5))

    # values in original were not updated
    npt.assert_equal(getattr(adata, key)["awk"]["c"], np.full((3, 1), 3))
    with pytest.raises(IndexError):
        getattr(adata, key)["awk"]["d"]


@pytest.mark.parametrize("key", ["obsm", "varm"])
def test_view(key):
    """Check that modifying a view does not modify the original"""
    adata = gen_adata((3, 3), varm_types=(), obsm_types=(), layers_types=())
    getattr(adata, key)["awk"] = ak.Array([{"a": [1], "b": [2], "c": [3]}] * 3)
    adata_view = adata[:2, :2]

    with pytest.warns(ImplicitModificationWarning, match="initializing view as actual"):
        getattr(adata_view, key)["awk"]["c"] = np.full((2, 1), 4)
        getattr(adata_view, key)["awk"]["d"] = np.full((2, 1), 5)

    # values in view were correctly set
    npt.assert_equal(getattr(adata_view, key)["awk"]["c"], np.full((2, 1), 4))
    npt.assert_equal(getattr(adata_view, key)["awk"]["d"], np.full((2, 1), 5))

    # values in original were not updated
    npt.assert_equal(getattr(adata, key)["awk"]["c"], np.full((3, 1), 3))
    with pytest.raises(IndexError):
        getattr(adata, key)["awk"]["d"]
