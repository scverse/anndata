from __future__ import annotations

from string import ascii_lowercase

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData


@pytest.fixture
def adata() -> AnnData:
    gen = np.random.default_rng()
    x = gen.random((100, 50), dtype=np.float32)
    layers = dict(a=sp.random(100, 50, format="csr"))
    obs = pd.DataFrame(
        dict(type=gen.integers(0, 3, size=100)),
        index="cell-" + pd.array(range(100)).astype(str),
    )
    var_grp = pd.Categorical(
        gen.integers(0, 6, size=50), categories=list(ascii_lowercase[:5])
    )
    var = pd.DataFrame(
        dict(grp=var_grp),
        index="gene-" + pd.array(range(50)).astype(str),
    )
    obsm = dict(umap=gen.random((100, 2)))
    varp = dict(cons=sp.csr_array(sp.random(50, 50)))
    return AnnData(x, obs, var, layers=layers, obsm=obsm, varm={}, obsp={}, varp=varp)
