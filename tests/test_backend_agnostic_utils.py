from __future__ import annotations

import jax
import numpy as np
import pytest
from benchmarks.benchmarks.utils import gen_adata, gen_indexer

from anndata import AnnData


@pytest.mark.parametrize("backend", ["numpy", "jax"])
def test_gen_adata_and_indexing(backend):
    # Generate AnnData using backend
    if backend == "numpy":
        pass  # default backend used by gen_adata
    elif backend == "jax":
        jnp = jax.numpy
        _ = jnp.ones((1,))
    else:
        raise ValueError(f"Unsupported backend: {backend}")

    adata: AnnData = gen_adata(100, 50, {"X-dense", "obs,var"}, seed=42)

    # Check structure
    assert adata.shape == (100, 50)
    assert "a" in adata.obs.columns
    assert "a" in adata.var.columns

    # Test each index kind
    for kind in ["slice", "intarray", "boolarray", "strarray"]:
        subset = gen_indexer(adata, "obs", kind, 0.3, seed=123)
        assert isinstance(subset, tuple)
        assert len(subset) == 2

        index = subset[0]
        if kind == "slice":
            assert isinstance(index, slice)
        elif kind == "intarray":
            assert hasattr(index, "shape")
            assert 0 < index.shape[0] <= 100
        elif kind == "boolarray":
            assert index.shape == (100,)
            assert index.dtype == bool
        elif kind == "strarray":
            assert isinstance(index, (list, np.ndarray))
            assert all(isinstance(i, str) for i in index)
