import re
from string import ascii_letters

import numpy as np
import pandas as pd
import pytest

import anndata as ad

ADATA_ATTRS = ("obs", "var", "varm", "obsm", "layers", "obsp", "varp", "uns")


@pytest.fixture
def adata():
    return ad.AnnData(
        np.zeros((20, 10)),
        obs=pd.DataFrame(
            dict(obs_key=list(ascii_letters[:20])),
            index=[f"cell{i}" for i in range(20)],
        ),
        var=pd.DataFrame(
            dict(var_key=np.arange(10)), index=[f"gene{i}" for i in range(10)]
        ),
        varm=dict(varm_key=np.zeros((10, 20))),
        obsm=dict(obsm_key=np.zeros((20, 20))),
        layers=dict(layers_key=np.zeros((20, 10))),
        obsp=dict(obsp_key=np.zeros((20, 20))),
        varp=dict(varp_key=np.zeros((10, 10))),
        uns=dict(uns_key=dict(zip("abc", range(3)))),
    )


@pytest.fixture(params=ADATA_ATTRS)
def adata_attr(request):
    return request.param


def test_anndata_repr(adata):
    assert f"{adata.n_obs} × {adata.n_vars}" in repr(adata)

    for idxr in [
        (slice(10, 20), 12),
        (12, 10),
        (["cell1", "cell2"], slice(10, 15)),
    ]:
        v = adata[idxr]
        v_repr = repr(v)
        assert f"{v.n_obs} × {v.n_vars}" in v_repr
        assert "View of" in v_repr
        for attr in ADATA_ATTRS:
            assert re.search(
                rf"^\s+{attr}:[^$]*{attr}_key.*$", v_repr, flags=re.MULTILINE
            )


def test_removal(adata, adata_attr):
    attr = adata_attr
    assert re.search(rf"^\s+{attr}:.*$", repr(adata), flags=re.MULTILINE)
    delattr(adata, attr)
    assert re.search(rf"^\s+{attr}:.*$", repr(adata), flags=re.MULTILINE) is None
