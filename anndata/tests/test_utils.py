from contextlib import suppress

import pytest
import zarr
import h5py
import pandas as pd

from anndata.compat import _clean_uns
from anndata._io.utils import report_read_key_on_error, AnnDataReadError


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
