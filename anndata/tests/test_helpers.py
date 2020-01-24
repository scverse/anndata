from string import ascii_letters

import pandas as pd
import pytest
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import assert_equal, report_name, gen_adata

# Testing to see if all error types can have the key name appended.
# Currently fails for 22/118 since they have required arguments. Not sure what to do about that.
#
# @singledispatch
# def iswarning(x):
#     return iswarning(type(x))

# @iswarning.register(type)
# def _notwarning(x):
#     return False

# @iswarning.register(Warning)
# def _iswarning(x):
#     return True

# @pytest.mark.parametrize("exception", list(filter(lambda t: not iswarning(t), Exception.__subclasses__())))
# def test_report_name_types(exception):
#     def throw(e):
#         raise e()
#     tag = "".join(np.random.permutation(list(ascii_letters)))

#     with pytest.raises(exception) as err:
#         report_name(throw)(exception, _elem_name=tag)
#     assert tag in str(err.value)


@pytest.fixture(scope="function")
def reusable_adata():
    """Reusable anndata for when tests shouldn’t mutate it"""
    return gen_adata((10, 10))


# Does this work for every warning?
def test_report_name():
    def raise_error():
        raise Exception("an error occured!")

    letters = np.array(list(ascii_letters))
    tag = "".join(np.random.permutation(letters))
    with pytest.raises(Exception) as e1:
        raise_error()
    with pytest.raises(Exception) as e2:
        report_name(raise_error)(_elem_name=tag)
    assert str(e2.value).startswith(str(e1.value))
    assert tag in str(e2.value)


def test_assert_equal():
    # ndarrays
    assert_equal(np.ones((10, 10)), np.ones((10, 10)))
    assert_equal(  # Should this require an exact test?
        np.ones((10, 10), dtype="i8"), np.ones((10, 10), dtype="f8")
    )
    assert_equal(
        np.array(list(ascii_letters)), np.array(list(ascii_letters)), exact=True
    )
    with pytest.raises(AssertionError):
        assert_equal(np.array(list(ascii_letters)), np.array(list(ascii_letters))[::-1])

    adata = gen_adata((10, 10))
    adata.raw = adata.copy()
    assert_equal(adata, adata.copy(), exact=True)
    # TODO: I’m not sure this is good behaviour, I’ve disabled in for now.
    # assert_equal(
    #     adata,
    #     adata[
    #         np.random.permutation(adata.obs_names),
    #         np.random.permutation(adata.var_names),
    #     ].copy(),
    #     exact=False,
    # )
    adata2 = adata.copy()
    to_modify = list(adata2.layers.keys())[0]
    del adata2.layers[to_modify]
    with pytest.raises(AssertionError) as missing_layer_error:
        assert_equal(adata, adata2)
    assert "layers" in str(missing_layer_error.value)
    # `to_modify` will be in pytest info
    adata2 = adata.copy()
    adata2.layers[to_modify][0, 0] = adata2.layers[to_modify][0, 0] + 1
    with pytest.raises(AssertionError) as changed_layer_error:
        assert_equal(adata, adata2)
    assert "layers" in str(changed_layer_error.value)
    assert to_modify in str(changed_layer_error.value)

    assert_equal(adata.obs, adata.obs.copy(), exact=True)

    csr = sparse.random(100, 100, format="csr")
    csc = csr.tocsc()
    dense = csr.toarray()
    assert_equal(csr, csc)
    assert_equal(csc, dense)
    assert_equal(dense, csc)


def test_assert_equal_raw():
    base = gen_adata((10, 10))
    orig = base.copy()
    orig.raw = base.copy()
    mod = base.copy()
    mod.X[0, 0] = mod.X[0, 0] + 1
    to_compare = base.copy()
    to_compare.raw = mod.copy()
    with pytest.raises(AssertionError):
        assert_equal(orig, to_compare)

    mod = base.copy()
    mod.var["new_val"] = 1
    to_compare = base.copy()
    to_compare.raw = mod.copy()
    with pytest.raises(AssertionError):
        assert_equal(orig, to_compare)


# TODO: Should views be equal to actual?
# Should they not be if an exact comparison is made?
def test_assert_equal_aligned_mapping():
    adata1 = gen_adata((10, 10))
    adata2 = adata1.copy()

    for attr in ["obsm", "varm", "layers", "obsp", "varp"]:
        assert_equal(getattr(adata1, attr), getattr(adata2, attr))

    # Checking that subsetting other axis only changes some attrs
    obs_subset = adata2[:5, :]
    for attr in ["obsm", "layers", "obsp"]:
        with pytest.raises(AssertionError):
            assert_equal(getattr(adata1, attr), getattr(obs_subset, attr))
    for attr in ["varm", "varp"]:
        assert_equal(getattr(adata1, attr), getattr(obs_subset, attr))

    var_subset = adata2[:, 5:]
    for attr in ["varm", "layers", "varp"]:
        with pytest.raises(AssertionError):
            assert_equal(getattr(adata1, attr), getattr(var_subset, attr))
    for attr in ["obsm", "obsp"]:
        assert_equal(getattr(adata1, attr), getattr(var_subset, attr))


def test_assert_equal_aligned_mapping_empty():
    chars = np.array(list(ascii_letters))
    adata = ad.AnnData(
        X=np.zeros((10, 10)),
        obs=pd.DataFrame([], index=np.random.choice(chars[:20], 10, replace=False)),
        var=pd.DataFrame([], index=np.random.choice(chars[:20], 10, replace=False)),
    )
    diff_idx = ad.AnnData(
        X=np.zeros((10, 10)),
        obs=pd.DataFrame([], index=np.random.choice(chars[20:], 10, replace=False)),
        var=pd.DataFrame([], index=np.random.choice(chars[20:], 10, replace=False)),
    )
    same_idx = ad.AnnData(adata.X, obs=adata.obs.copy(), var=adata.var.copy())

    for attr in ["obsm", "varm", "layers", "obsp", "varp"]:
        with pytest.raises(AssertionError):
            assert_equal(getattr(adata, attr), getattr(diff_idx, attr))
        assert_equal(getattr(adata, attr), getattr(same_idx, attr))
