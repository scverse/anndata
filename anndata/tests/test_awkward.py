"""Tests related to awkward arrays"""
import pytest
import numpy as np
import numpy.testing as npt

from anndata.tests.helpers import assert_equal, gen_adata, gen_awkward
from anndata.compat import awkward as ak, awkward_version
from anndata import ImplicitModificationWarning
from anndata.utils import dim_len
from anndata import AnnData, read_h5ad


@pytest.mark.parametrize(
    "array,shape",
    [
        # numpy array
        [ak.Array(np.arange(2 * 3 * 4 * 5).reshape((2, 3, 4, 5))), (2, 3, 4, 5)],
        # record
        [ak.Array([{"a": 1, "b": 2}, {"a": 1, "b": 3}]), (2,)],
        # ListType, variable length
        [ak.Array([[1], [2, 3], [4, 5, 6]]), (3, None)],
        # ListType, happens to have the same length, but is not regular
        [ak.Array([[2], [3], [4]]), (3, None)],
        # RegularType + nested ListType
        [ak.to_regular(ak.Array([[[1, 2], [3]], [[2], [3, 4, 5]]]), 1), (2, 2, None)],
        # nested record
        [
            ak.to_regular(ak.Array([[{"a": 0}, {"b": 1}], [{"c": 2}, {"d": 3}]]), 1),
            (2, 2),
        ],
        # mixed types (variable length)
        [ak.Array([[1, 2], ["a"]]), (2, None)],
        # mixed types (but regular)
        [ak.to_regular(ak.Array([[1, 2], ["a", "b"]]), 1), (2, 2)],
        # zero-size edge cases
        [ak.Array(np.ones((0, 7))), (0, 7)],
        [ak.Array(np.ones((7, 0))), (7, 0)],
        # UnionType of two regular types with different dimensions
        [
            ak.concatenate([ak.Array(np.ones((2, 2))), ak.Array(np.ones((2, 3)))]),
            (4, None),
        ],
        # UnionType of two regular types with same dimension
        [
            ak.concatenate(
                [
                    ak.Array(np.ones((2, 2))),
                    ak.Array(np.array([["a", "a"], ["a", "a"]])),
                ]
            ),
            (4, 2),
        ],
    ],
)
def test_dim_len(array, shape):
    """Test that dim_len returns the right value for awkward arrays."""
    for axis, size in enumerate(shape):
        assert size == dim_len(array, axis)

    # Requesting the size for an axis higher than the array has dimensions should raise a TypeError
    with pytest.raises(TypeError):
        dim_len(array, len(shape))


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
@pytest.mark.skipif(
    awkward_version >= 2, reason="This test is only applies for awkward versions 1.x"
)
def test_no_awkward_v1(key):
    """Assigning an awkward v1 array to anndata should fail"""
    import awkward as akv1

    v1arr = akv1.Array([1, 2, 3, 4])

    adata = AnnData(np.ones((4, 4)))
    with pytest.raises(AttributeError):
        getattr(adata, key)["test"] = v1arr


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


@pytest.mark.parametrize(
    "array",
    [
        # numpy array
        ak.Array(np.arange(2 * 3 * 4 * 5).reshape((2, 3, 4, 5))),
        # record
        ak.Array([{"a": 1, "b": 2}, {"a": 1, "b": 3}]),
        # ListType, variable length
        ak.Array([[1], [2, 3], [4, 5, 6]]),
        # RegularType + nested ListType
        ak.to_regular(ak.Array([[[1, 2], [3]], [[2], [3, 4, 5]]]), 1),
        # nested record
        ak.to_regular(ak.Array([[{"a": 0}, {"b": 1}], [{"c": 2}, {"d": 3}]]), 1),
        # mixed types (variable length)
        ak.Array([[1, 2], ["a"]]),
        # zero-size edge cases
        ak.Array(np.ones((0, 7))),
        ak.Array(np.ones((7, 0))),
        # UnionType of two regular types with different dimensions
        ak.concatenate([ak.Array(np.ones((2, 2))), ak.Array(np.ones((2, 3)))]),
        # UnionType of two regular types with same dimension
        ak.concatenate(
            [
                ak.Array(np.ones((2, 2))),
                ak.Array(np.array([["a", "a"], ["a", "a"]])),
            ]
        ),
    ],
)
def test_awkward_io(tmp_path, array):
    adata = AnnData()
    adata.uns["awk"] = array
    adata_path = tmp_path / "adata.h5ad"
    adata.write_h5ad(adata_path)

    adata2 = read_h5ad(adata_path)

    assert_equal(adata.uns["awk"], adata2.uns["awk"])
