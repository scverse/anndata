from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
from anndata import AnnData

# making sure Jax and cubed are imported
jax = pytest.importorskip("jax")
cubed = pytest.importorskip("cubed")
import cubed.array_api as cxp  # noqa: E402
import jax.numpy as jnp  # noqa: E402


def assert_same_value_and_type(out_value, expected_array, expected_type):
    # Checking that the values are the same
    np.testing.assert_array_equal(np.asarray(out_value), np.asarray(expected_array))
    # Type matches exactly
    assert type(out_value) is expected_type


def test_concat_uns_jax_jax():
    # JAX + JAX  (values equal), type should stay JAX
    a1 = AnnData(np.ones((1, 1)))
    a2 = AnnData(np.ones((1, 1)))

    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = jnp.array([1, 2, 3])

    out = ad.concat([a1, a2], uns_merge="same")

    assert "k" in out.uns
    arr = out.uns["k"]
    assert getattr(arr, "shape", None) == (3,)
    # Values match
    np.testing.assert_array_equal(np.asarray(arr), np.array([1, 2, 3]))
    # Type must be the same as the input (JAX array)
    assert type(arr) is type(a1.uns["k"])


def test_concat_uns_cubed_cubed():
    # Cubed + Cubed (values equal), type should stay Cubed
    a1 = AnnData(np.ones((1, 1)))
    a2 = AnnData(np.ones((1, 1)))

    a1.uns["k"] = cxp.asarray([1, 2, 3])
    a2.uns["k"] = cxp.asarray([1, 2, 3])

    out = ad.concat([a1, a2], uns_merge="same")
    assert "k" in out.uns

    assert_same_value_and_type(out.uns["k"], cxp.asarray([1, 2, 3]), type(a1.uns["k"]))


def test_concat_uns_jax_numpy():
    # JAX + NumPy (values equal), keeps the first type (JAX)
    a1 = AnnData(np.ones((1, 1)))
    a2 = AnnData(np.ones((1, 1)))

    a1.uns["k"] = jnp.array([7, 8])
    a2.uns["k"] = np.array([7, 8])

    out = ad.concat([a1, a2], uns_merge="same")
    assert "k" in out.uns

    # Keep first's value and type (JAX)
    assert_same_value_and_type(out.uns["k"], jnp.array([7, 8]), type(a1.uns["k"]))


def test_concat_uns_jax_cubed():
    # JAX + Cubed (values equal), first object type (JAX) in result
    a1 = AnnData(np.ones((1, 1)))
    a2 = AnnData(np.ones((1, 1)))

    a1.uns["k"] = jnp.array([4, 5, 6])
    a2.uns["k"] = cxp.asarray([4, 5, 6])

    out = ad.concat([a1, a2], uns_merge="same")
    assert "k" in out.uns

    assert_same_value_and_type(out.uns["k"], jnp.array([4, 5, 6]), type(a1.uns["k"]))


def test_concat_uns_cubed_numpy():
    # Cubed + NumPy (values equal), first object type (Cubed) in result
    a1 = AnnData(np.ones((1, 1)))
    a2 = AnnData(np.ones((1, 1)))

    a1.uns["k"] = cxp.asarray([9, 10, 11])
    a2.uns["k"] = np.array([9, 10, 11])

    out = ad.concat([a1, a2], uns_merge="same")
    assert "k" in out.uns

    print(out.uns["k"])
    print(type(out.uns["k"]))

    assert_same_value_and_type(
        out.uns["k"], cxp.asarray([9, 10, 11]), type(a1.uns["k"])
    )


def test_uns_unique_jax_jax_equal_kept():
    # unique merge, equal JAX arrays, keeps values and type
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = jnp.array([1, 2, 3])

    # Concatenation with unique merge
    out = ad.concat([a1, a2], uns_merge="unique")
    assert "k" in out.uns
    arr = out.uns["k"]
    assert arr.shape == (3,)
    # check that the values are the same
    np.testing.assert_array_equal(np.asarray(arr), np.array([1, 2, 3]))
    assert type(arr) is type(a1.uns["k"])  # stays JAX


def test_uns_unique_jax_jax_conflict_dropped():
    # unique merge with JAX arrays that have differing values, so will be dropped
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = jnp.array([3, 2, 5])

    out = ad.concat([a1, a2], uns_merge="unique")
    # differing values, so will be dropped
    assert "k" not in out.uns


def test_uns_unique_jax_numpy_equal_kept():
    # unique merge with JAX and NumPy arrays that are equal, so will be kept and first type (JAX)
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([7, 8])
    a2.uns["k"] = np.array([7, 8])

    out = ad.concat([a1, a2], uns_merge="unique")
    assert "k" in out.uns
    assert_same_value_and_type(
        out.uns["k"], jnp.array([7, 8]), type(a1.uns["k"])
    )  # first wins, JAX type


def test_uns_first_jax_jax_conflict_keeps():
    # first merge with JAX arrays that have differing values, so will keep the first
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = jnp.array([3, 2, 5])

    out = ad.concat([a1, a2], uns_merge="first")
    assert "k" in out.uns
    assert_same_value_and_type(
        out.uns["k"], jnp.array([1, 2, 3]), type(a1.uns["k"])
    )  # first wins, JAX type


def test_uns_first_jax_numpy_conflict_keeps_first():
    # first merge with JAX and NumPy arrays that have differing values, so will keep the first
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = np.array([9, 9, 9])

    out = ad.concat([a1, a2], uns_merge="first")
    assert "k" in out.uns
    assert_same_value_and_type(
        out.uns["k"], jnp.array([1, 2, 3]), type(a1.uns["k"])
    )  # JAX stays


def test_uns_first_cubed_numpy_conflict_keeps_first():
    # first merge with Cubed and NumPy arrays that have differing values, so will keep the first
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = cxp.asarray([1, 2, 3])
    a2.uns["k"] = np.array([9, 9, 9])

    out = ad.concat([a1, a2], uns_merge="first")
    assert "k" in out.uns
    assert_same_value_and_type(
        out.uns["k"], cxp.asarray([1, 2, 3]), type(a1.uns["k"])
    )  # Cubed stays


def test_uns_only_keeps_when_present_in_only_one_input_jax():
    # only merge with JAX arrays, keeps the one present in only one input
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    a1.uns["k"] = jnp.array([1, 2, 3])
    # a2 has no 'k'

    out = ad.concat([a1, a2], uns_merge="only")
    assert "k" in out.uns
    assert_same_value_and_type(
        out.uns["k"], jnp.array([1, 2, 3]), type(a1.uns["k"])
    )  # kept from only input


def test_uns_only_drops_when_present_in_both_inputs_even_equal_jax_without_normalization():
    # only merge with present in both inputs, so drops it
    a1, a2 = AnnData(np.ones((1, 1))), AnnData(np.ones((1, 1)))
    # present in both inputs
    a1.uns["k"] = jnp.array([1, 2, 3])
    a2.uns["k"] = jnp.array([1, 2, 3])

    out = ad.concat([a1, a2], uns_merge="only")
    assert "k" not in out.uns
