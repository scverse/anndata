from __future__ import annotations

import gc
import sys
from string import ascii_lowercase
from time import sleep

import jax
import numpy as np
import pandas as pd
from array_api_compat import get_namespace as array_api_get_namespace
from memory_profiler import memory_usage
from scipy import sparse

from anndata import AnnData


def get_namespace(x=None):
    return array_api_get_namespace(x)


def get_rng(xp, seed=None):
    """Return a backend-specific random number generator."""
    # RNG isn't standardized in the Array API spec,
    # so backends like JAX, PyTorch, and NumPy each handle randomness differently.
    if xp.__name__.startswith("jax"):
        return jax.random.PRNGKey(seed or 0)
    elif xp.__name__.startswith("numpy"):
        return np.random.default_rng(seed)
    else:
        raise NotImplementedError(f"RNG not implemented for backend: {xp.__name__}")


def get_actualsize(input_obj):
    """Using Python Garbage Collector to calculate the size of all elements attached to an object"""

    memory_size = 0
    ids = set()
    objects = [input_obj]
    while objects:
        new = []
        for obj in objects:
            if id(obj) not in ids:
                ids.add(id(obj))
                memory_size += sys.getsizeof(obj)
                new.append(obj)
        objects = gc.get_referents(*new)
    return memory_size


def get_anndata_memsize(adata):
    recording = memory_usage(
        (sedate(adata.copy, naplength=0.005), (adata,)), interval=0.001
    )
    diff = recording[-1] - recording[0]
    return diff


def get_peak_mem(op, interval=0.001):
    recording = memory_usage(op, interval=interval)
    xp = get_namespace()
    return xp.max(recording) - xp.min(recording)


def sedate(func, naplength=0.05):
    """Make a function sleepy, so we can sample the start and end state."""

    def wrapped_function(*args, **kwargs):
        sleep(naplength)
        val = func(*args, **kwargs)
        sleep(naplength)
        return val

    return wrapped_function


# TODO: Factor out the time it takes to generate these


def gen_indexer(adata, dim, index_kind, ratio, seed=None):
    dimnames = ("obs", "var")
    index_kinds = {"slice", "intarray", "boolarray", "strarray"}

    if index_kind not in index_kinds:
        msg = f"Argument 'index_kind' must be one of {index_kinds}. Was {index_kind}."
        raise ValueError(msg)

    xp = get_namespace(adata.X)
    rng = get_rng(xp, seed)
    axis = dimnames.index(dim)
    subset = [slice(None), slice(None)]
    axis_size = adata.shape[axis]
    n = int(xp.round(axis_size * ratio))

    if index_kind == "slice":
        subset[axis] = slice(0, n))
    elif index_kind == "intarray":
        if xp.__name__.startswith("jax"):
            subset[axis] = jax.random.choice(rng, xp.arange(axis_size), shape=(n,), replace=False)
        elif xp.__name__.startswith("numpy"):
            subset[axis] = xp.asarray(rng.choice(axis_size, n, replace=False))

    elif index_kind == "boolarray":
        mask = xp.zeros(axis_size, dtype=bool)
        if xp.__name__.startswith("jax"):
            idx = jax.random.choice(rng, xp.arange(axis_size), shape=(n,), replace=False)
        elif xp.__name__.startswith("numpy"):
            idx = rng.choice(axis_size, n, replace=False)
        mask[idx] = True
        subset[axis] = mask

    elif index_kind == "strarray":
        subset[axis] = rng.choice(
            getattr(adata, dim).index, n, replace=False
        )
    else:
        raise ValueError()
    return tuple(subset)


def gen_adata(n_obs, n_var, attr_set, seed=None):
    xp = get_namespace()
    rng = get_rng(xp, seed)
    if "X-csr" in attr_set:
        X = sparse.random(n_obs, n_var, density=0.1, format="csr")
    elif "X-dense" in attr_set:
        X = sparse.random(n_obs, n_var, density=0.1, format="csr")
        X = xp.asarray(X.toarray())
    else:
        # TODO: There's probably a better way to do this
        X = sparse.random(n_obs, n_var, density=0, format="csr")
    adata = AnnData(X)
    if "obs,var" in attr_set:
        if xp.__name__.startswith("jax"):
            obs = {k: jax.random.randint(rng, (n_obs,), 0, 100) for k in ascii_lowercase}
            var = {k: jax.random.randint(rng, (n_var,), 0, 100) for k in ascii_lowercase}
        elif xp.__name__.startswith("numpy"):
            obs = {k: rng.integers(0, 100, size=n_obs) for k in ascii_lowercase}
            var = {k: rng.integers(0, 100, size=n_var) for k in ascii_lowercase}
        adata.obs = pd.DataFrame(obs, index=[f"cell{i}" for i in range(n_obs)])
        adata.var = pd.DataFrame(var, index=[f"gene{i}" for i in range(n_var)])
    return adata
