from __future__ import annotations

from contextlib import suppress

import numpy as np
import pytest

import anndata as ad
from anndata import AnnData

# getting all available backends
# start with numpy as default
backends: list[tuple[str, callable]] = [("numpy", np.array)]

with suppress(ImportError):
    # for now it is here
    import jax.numpy as jnp

    backends.append(("jax", jnp.array))

with suppress(ImportError):
    import cubed.array_api as cxp

    # future implementation
    backends.append(("cubed", cxp.asarray))


def assert_same_value_and_type(out_value, expected_array, expected_type):
    """Helper to check values and backend type."""
    np.testing.assert_array_equal(np.asarray(out_value), np.asarray(expected_array))
    assert type(out_value) is expected_type


@pytest.mark.parametrize(("name1", "asarray1"), backends)
@pytest.mark.parametrize(("name2", "asarray2"), backends)
def test_concat_uns_same_equal(name1, asarray1, name2, asarray2):
    # Equal arrays: values match, first type wins.
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray1([1, 2, 3])
    a2.uns["k"] = asarray2([1, 2, 3])

    out = ad.concat([a1, a2], uns_merge="same")

    assert "k" in out.uns
    assert_same_value_and_type(out.uns["k"], [1, 2, 3], type(a1.uns["k"]))


@pytest.mark.parametrize(("name1", "asarray1"), backends)
@pytest.mark.parametrize(("name2", "asarray2"), backends)
def test_concat_uns_unique_equal_kept(name1, asarray1, name2, asarray2):
    # Equal arrays across backends: kept, first type wins.
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray1([7, 8])
    a2.uns["k"] = asarray2([7, 8])

    out = ad.concat([a1, a2], uns_merge="unique")

    assert "k" in out.uns
    assert_same_value_and_type(out.uns["k"], [7, 8], type(a1.uns["k"]))


@pytest.mark.parametrize(("name1", "asarray1"), backends)
@pytest.mark.parametrize(("name2", "asarray2"), backends)
def test_concat_uns_unique_conflict_dropped(name1, asarray1, name2, asarray2):
    """Different arrays: dropped."""
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray1([1, 2, 3])
    a2.uns["k"] = asarray2([9, 9, 9])

    out = ad.concat([a1, a2], uns_merge="unique")

    assert "k" not in out.uns


@pytest.mark.parametrize(("name1", "asarray1"), backends)
@pytest.mark.parametrize(("name2", "asarray2"), backends)
def test_concat_uns_first_conflict_keeps_first(name1, asarray1, name2, asarray2):
    """Conflict arrays: keep the first's values and type."""
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray1([1, 2, 3])
    a2.uns["k"] = asarray2([9, 9, 9])

    out = ad.concat([a1, a2], uns_merge="first")

    assert "k" in out.uns
    assert_same_value_and_type(out.uns["k"], [1, 2, 3], type(a1.uns["k"]))


@pytest.mark.parametrize(("name", "asarray"), backends)
def test_concat_uns_only_singleton_kept(name, asarray):
    """Present in only one AnnData: kept."""
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray([1, 2, 3])  # only in a1

    out = ad.concat([a1, a2], uns_merge="only")

    assert "k" in out.uns
    assert_same_value_and_type(out.uns["k"], [1, 2, 3], type(a1.uns["k"]))


@pytest.mark.parametrize(("name1", "asarray1"), backends)
@pytest.mark.parametrize(("name2", "asarray2"), backends)
def test_concat_uns_only_conflict_drops(name1, asarray1, name2, asarray2):
    # Present in both: dropped even if equal.
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = asarray1([1, 2, 3])
    a2.uns["k"] = asarray2([1, 2, 3])

    out = ad.concat([a1, a2], uns_merge="only")

    assert "k" not in out.uns
