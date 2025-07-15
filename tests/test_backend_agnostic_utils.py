# from __future__ import annotations

# import jax
# import numpy as np
# import pytest
# from benchmarks.tests.helpers import gen_adata, gen_indexer

# from anndata import AnnData


# @pytest.mark.parametrize("backend", ["numpy", "jax"])
# def test_gen_adata_and_indexing(backend):
#     # Generate AnnData using backend
#     if backend == "numpy":
#         pass  # default backend used by gen_adata
#     elif backend == "jax":
#         jnp = jax.numpy
#         _ = jnp.ones((1,))
#     else:
#         raise ValueError(f"Unsupported backend: {backend}")

#     adata: AnnData = gen_adata(100, 50, {"X-dense", "obs,var"}, seed=42)

#     # Check structure
#     assert adata.shape == (100, 50)
#     assert "a" in adata.obs.columns
#     assert "a" in adata.var.columns

#     # Test each index kind
#     for kind in ["slice", "intarray", "boolarray", "strarray"]:
#         subset = gen_indexer(adata, "obs", kind, 0.3, seed=123)
#         assert isinstance(subset, tuple)
#         assert len(subset) == 2

#         index = subset[0]
#         if kind == "slice":
#             assert isinstance(index, slice)
#         elif kind == "intarray":
#             assert hasattr(index, "shape")
#             assert 0 < index.shape[0] <= 100
#         elif kind == "boolarray":
#             assert index.shape == (100,)
#             assert index.dtype == bool
#         elif kind == "strarray":
#             assert isinstance(index, (list, np.ndarray))
#             assert all(isinstance(i, str) for i in index)

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from anndata._core.anndata import AnnData
from anndata.tests.helpers import gen_adata

# Try to import JAX if available
# flagging it as a separate import to avoid issues if JAX is not installed
try:
    import jax
    import jax.numpy as jnp

    jax_available = True
except ImportError:
    jax_available = False


# testing gen_adata with NumPy backend and various X types (dense, CSR, CSC) to ensure correct shapes and AnnData validity
@pytest.mark.parametrize("X_type", [np.array, sparse.csr_matrix, sparse.csc_matrix])
def test_gen_adata_numpy_backends(X_type):
    adata = gen_adata(
        shape=(20, 30),
        X_type=X_type,
        X_dtype=np.float32,
        xp=np,
        rng=np.random.default_rng(0),
    )
    assert isinstance(adata, AnnData)
    assert adata.X.shape == (20, 30)
    assert adata.obs.shape[0] == 20
    assert adata.var.shape[0] == 30


# testing that gen_adata works with JAX backend when X is omitted
@pytest.mark.skipif(not jax_available, reason="JAX is not available")
def test_gen_adata_jax_backend_no_X():
    adata = gen_adata(
        shape=(20, 30),
        X_type=None,  # skip X
        xp=jax.numpy,
        rng=(jax.random, jax.random.PRNGKey(0)),
    )
    assert isinstance(adata, AnnData)
    assert adata.X is None
    assert adata.obs.shape[0] == 20
    assert adata.var.shape[0] == 30


# checking if function correctly returns an Anndata with X as None
def test_gen_adata_X_none():
    adata = gen_adata(shape=(10, 10), X_type=None)
    assert isinstance(adata, AnnData)
    assert adata.X is None


# testing that passing a sparse matrix format to JAX correctly raises a ValueError
def test_gen_adata_invalid_jax_sparse():
    if not jax_available:
        pytest.skip("JAX not available")
    with pytest.raises(ValueError, match="JAX does not support sparse matrices"):
        gen_adata(
            shape=(5, 5),
            X_type=sparse.csr_matrix,
            xp=jax.numpy,
            rng=(jax.random, jax.random.PRNGKey(0)),
        )
