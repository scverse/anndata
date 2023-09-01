"""Tests related to awkward arrays"""
import pytest
import numpy as np
import numpy.testing as npt

from anndata.tests.helpers import assert_equal, gen_adata, gen_awkward
from anndata.compat import awkward as ak
from anndata import ImplicitModificationWarning
from anndata.utils import dim_len
from anndata import AnnData, read_h5ad
import anndata
import pandas as pd


@pytest.mark.parametrize(
    "array,shape",
    [
        # numpy array
        [ak.Array(np.arange(2 * 3 * 4 * 5).reshape((2, 3, 4, 5))), (2, 3, 4, 5)],
        # record
        [ak.Array([{"a": 1, "b": 2}, {"a": 1, "b": 3}]), (2, 2)],
        # ListType, variable length
        [ak.Array([[1], [2, 3], [4, 5, 6]]), (3, None)],
        # ListType, happens to have the same length, but is not regular
        [ak.Array([[2], [3], [4]]), (3, None)],
        # RegularType + nested ListType
        [ak.to_regular(ak.Array([[[1, 2], [3]], [[2], [3, 4, 5]]]), 1), (2, 2, None)],
        # nested record
        [
            ak.to_regular(ak.Array([[{"a": 0}, {"b": 1}], [{"c": 2}, {"d": 3}]]), 1),
            (2, 2, 4),
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
        # Array of string types
        [ak.Array(["a", "b", "c"]), (3,)],
        [ak.Array([["a", "b"], ["c", "d"], ["e", "f"]]), (3, None)],
        [ak.to_regular(ak.Array([["a", "b"], ["c", "d"], ["e", "f"]]), 1), (3, 2)],
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
    """Check if we can set obsm, .varm and .uns with different types
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


@pytest.mark.parametrize(
    "array,setter_slice,setter_value,expected",
    [
        # Non-records are immutable and setting on them results in a TypeError
        pytest.param(
            [[1], [2, 3], [4, 5, 6]],
            (slice(None), 2),
            42,
            TypeError,
            id="immutable_ragged_list",
        ),
        pytest.param(
            np.zeros((3, 3)),
            (slice(None), 1),
            42,
            TypeError,
            id="immutable_regular_type",
        ),
        pytest.param(
            [{"a": 1}, {"a": 2}, {"a": 3}],
            "a",
            [42, 43, 44],
            [{"a": 42}, {"a": 43}, {"a": 44}],
            id="updating_record",
        ),
        pytest.param(
            [{"a": 1}, {"a": 2}, {"a": 3}],
            "b",
            [42, 43, 44],
            [{"a": 1, "b": 42}, {"a": 2, "b": 43}, {"a": 3, "b": 44}],
            id="adding_record",
        ),
    ],
)
@pytest.mark.parametrize("key", ["obsm", "varm", "uns"])
def test_view(key, array, setter_slice, setter_value, expected):
    """Check that modifying a view does not modify the original.

    Parameters
    ----------
    key
        key in anndata, obsm, varm, or uns
    view_func
        a function that returns a view of an AnnData object
    array
        The array that is assigned to adata[key]["awk"] for testing.
    setter_slice
        The slice used for setting a value on the awkward array with `arr[slice] = ...`
    setter_value
        The value assigned to the array with `arr[slice] = setter_value`
    expected
        The expected array after setting the value. Can be an exception if setting the value is supposed
        to result in an error.
    """
    adata = gen_adata((3, 3), varm_types=(), obsm_types=(), layers_types=())
    getattr(adata, key)["awk"] = ak.Array(array)
    adata_view = adata[:]
    awk_view = getattr(adata_view, key)["awk"]

    if isinstance(expected, type):
        with pytest.raises(expected):
            awk_view[setter_slice] = setter_value
    else:
        awk_view[setter_slice] = setter_value
        # values in view are correctly set
        assert ak.all(awk_view[setter_slice] == setter_value)
        # values in original were not modified
        assert ak.to_list(getattr(adata, key)["awk"]) == array


def test_view_of_awkward_array_with_custom_behavior():
    """Ensure that a custom behavior persists when creating a view."""
    from uuid import uuid4

    BEHAVIOUR_ID = str(uuid4())

    class ReversibleArray(ak.Array):
        def reversed(self):
            return self[..., ::-1]

    ak.behavior[BEHAVIOUR_ID] = ReversibleArray
    ak.behavior["*", BEHAVIOUR_ID] = ReversibleArray

    adata = gen_adata((3, 3), varm_types=(), obsm_types=(), layers_types=())
    adata.obsm["awk_string"] = ak.with_name(
        ak.Array([{"a": "AAA"}, {"a": "BBB"}, {"a": "CCC"}]), BEHAVIOUR_ID
    )
    ak_view = adata[1:]

    assert ak.to_list(ak_view.obsm["awk_string"].reversed()["a"]) == ["CCC", "BBB"]


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
        # categorical array
        ak.to_categorical(ak.Array([["a", "b", "c"], ["a", "b"]])),
        ak.to_categorical(ak.Array([[1, 1, 2], [3, 3]])),
        # tyical record type with AIRR data consisting of different dtypes
        ak.Array(
            [
                [
                    {
                        "v_call": "TRV1",
                        "junction_aa": "ADDEEKK",
                        "productive": True,
                        "locus": None,
                        "consensus_count": 3,
                    },
                    {
                        "v_call": "TRV2",
                        "productive": False,
                        "locus": "TRA",
                        "consensus_count": 4,
                    },
                ],
                [
                    {
                        "v_call": None,
                        "junction_aa": "ADDEKK",
                        "productive": None,
                        "locus": "IGK",
                        "consensus_count": 3,
                    }
                ],
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

    assert_equal(adata.uns["awk"], adata2.uns["awk"], exact=True)


# @pytest.mark.parametrize("join", ["outer", "inner"])
@pytest.mark.parametrize(
    "arrays,join,expected",
    [
        pytest.param(
            [ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]), None],
            "inner",
            None,
            id="awk:recordoflists_null-inner",
        ),
        pytest.param(
            [ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]), None],
            "outer",
            ak.Array(
                [{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}, None, None, None]
            ),
            # maybe should return: ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}, {}, {}, {}]),
            id="awk:recordoflists_null-outer",
        ),
        pytest.param(
            [ak.Array([[{"a": 1}, {"a": 2}], []]), None],
            "outer",
            ak.Array([[{"a": 1}, {"a": 2}], [], None, None, None]),
            # maybe should return: ak.Array([[{"a": 1}, {"a": 2}], [], [], []]),
            id="awk:listofrecords_null-outer",
        ),
        pytest.param(
            [None, ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}])],
            "inner",
            None,
            id="null_awk-inner",
        ),
        pytest.param(
            [None, ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}])],
            "outer",
            ak.Array(
                [None, None, None, {"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]
            ),
            # maybe should return: ak.Array([{}, {}, {}, {"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
            id="null_awk:recordoflists-outer",
        ),
        pytest.param(
            [ak.Array([{"a": 1}, {"a": 2}]), ak.Array([{"a": 3}, {"a": 4}])],
            "inner",
            ak.Array([{"a": i} for i in range(1, 5)]),
            id="awk-simple-record",
        ),
        pytest.param(
            [
                ak.Array([{"a": 1, "b": 1}, {"a": 2, "b": 2}]),
                ak.Array([{"a": 3}, {"a": 4}]),
            ],
            "inner",
            ak.Array([{"a": i} for i in range(1, 5)]),
            id="awk-simple-record-inner",
        ),
        # TODO:
        # pytest.param(
        #     [
        #         ak.Array([{"a": 1, "b": 1}, {"a": 2, "b": 2}]),
        #         ak.Array([{"a": 3}, {"a": 4}]),
        #     ],
        #     "outer",
        #     ak.Array([{"a": 1, "b": 1}, {"a": 2, "b": 2}, {"a": 3}, {"a": 4},]),
        #     id="awk-simple-record-outer",
        # ),
        pytest.param(
            [
                None,
                ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
                pd.DataFrame(),
            ],
            "outer",
            NotImplementedError,  # TODO: ak.Array([{}, {}, {}, {"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
            id="null_awk_empty-pd",
        ),
        pytest.param(
            [
                ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
                pd.DataFrame(),
            ],
            "outer",
            NotImplementedError,  # TODO: ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
            id="awk_empty-pd",
        ),
        pytest.param(
            [
                ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
                pd.DataFrame().assign(a=[3, 4], b=[5, 6]),
            ],
            "outer",  # TODO: Should try inner too if implemented
            NotImplementedError,
        ),
        pytest.param(
            [
                ak.Array([{"a": [1, 2], "b": [1, 2]}, {"a": [3], "b": [4]}]),
                np.ones((3, 2)),
            ],
            "outer",
            NotImplementedError,
        ),
    ],
)
@pytest.mark.parametrize("key", ["obsm", "varm"])
def test_concat_mixed_types(key, arrays, expected, join):
    """Test that concatenation of AwkwardArrays with arbitrary types, but zero length dimension
    or missing values works."""
    axis = 0 if key == "obsm" else 1

    to_concat = []
    cell_id, gene_id = 0, 0
    for a in arrays:
        shape = np.array([3, 3])  # default shape (in case of missing array)
        if a is not None:
            length = dim_len(a, 0)
            shape[axis] = length

        tmp_adata = gen_adata(
            tuple(shape), varm_types=(), obsm_types=(), layers_types=()
        )
        prev_cell_id, prev_gene_id = cell_id, gene_id
        cell_id, gene_id = cell_id + shape[0], gene_id + shape[1]
        tmp_adata.obs_names = pd.RangeIndex(prev_cell_id, cell_id).astype(str)
        tmp_adata.var_names = pd.RangeIndex(prev_gene_id, gene_id).astype(str)
        if a is not None:
            if isinstance(a, pd.DataFrame):
                a.set_index(
                    tmp_adata.obs_names if key == "obsm" else tmp_adata.var_names,
                    inplace=True,
                )
            getattr(tmp_adata, key)["test"] = a

        to_concat.append(tmp_adata)

    if isinstance(expected, type) and issubclass(expected, Exception):
        with pytest.raises(expected):
            anndata.concat(to_concat, axis=axis, join=join)
    else:
        print(to_concat)
        result_adata = anndata.concat(to_concat, axis=axis, join=join)
        result = getattr(result_adata, key).get("test", None)
        assert_equal(expected, result, exact=True)
