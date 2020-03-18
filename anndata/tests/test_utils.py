from contextlib import suppress
from itertools import repeat

import pytest
import zarr
import h5py
import pandas as pd
from scipy import sparse

import anndata as ad
from anndata.compat import _clean_uns
from anndata._io.utils import report_read_key_on_error, AnnDataReadError
from anndata.tests.helpers import gen_typed_df


@pytest.mark.parametrize(
    "group_fn",
    [
        pytest.param(lambda _: zarr.group(), id="zarr"),
        pytest.param(lambda p: h5py.File(p / "test.h5", mode="a"), id="h5py"),
    ],
)
def test_key_error(tmp_path, group_fn):
    @report_read_key_on_error
    def read_attr(_):
        raise NotImplementedError()

    group = group_fn(tmp_path)
    with group if hasattr(group, "__enter__") else suppress():
        group["X"] = [1, 2, 3]
        group.create_group("group")
        with pytest.raises(AnnDataReadError) as e:
            read_attr(group["X"])
        assert "'/X'" in str(e.value)
        with pytest.raises(AnnDataReadError) as e:
            read_attr(group["group"])
        assert "'/group'" in str(e.value)


def test_clean_uns():
    d = dict(
        uns=dict(species_categories=["a", "b"]),
        obs=dict(species=pd.Series([0, 1, 0])),
        var=dict(species=pd.Series([0, 1, 0, 2])),
    )
    _clean_uns(d)
    assert "species_categories" not in d["uns"]
    assert isinstance(d["obs"]["species"], pd.Categorical)
    assert d["obs"]["species"].tolist() == ["a", "b", "a"]
    # var’s categories were overwritten by obs’s,
    # which we can detect here because var has too high codes
    assert isinstance(d["var"]["species"], pd.Series)


def test_unique_indices():
    m, n = (10, 20)
    obs_index = pd.Index(repeat("a", m), name="obs")
    var_index = pd.Index(repeat("b", n), name="var")

    adata = ad.AnnData(
        X=sparse.random(m, n, format="csr"),
        obs=gen_typed_df(m, index=obs_index),
        var=gen_typed_df(n, index=var_index),
        obsm={"df": gen_typed_df(m, index=obs_index)},
        varm={"df": gen_typed_df(n, index=var_index)},
    )

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    assert adata.obs_names.name == "obs"
    assert adata.var_names.name == "var"

    assert len(pd.unique(adata.obs_names)) == m
    assert len(pd.unique(adata.var_names)) == n

    assert adata.obsm["df"].index.name == "obs"
    assert adata.varm["df"].index.name == "var"

    assert len(pd.unique(adata.obsm["df"].index)) == m
    assert len(pd.unique(adata.varm["df"].index)) == n

    v = adata[:5, :5]

    assert v.obs_names.name == "obs"
    assert v.var_names.name == "var"

    assert v.obsm["df"].index.name == "obs"
    assert v.varm["df"].index.name == "var"
