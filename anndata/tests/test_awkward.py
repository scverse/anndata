"""Tests related to awkward arrays"""
import pytest
from anndata.tests.helpers import gen_adata, gen_awkward


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
