from contextlib import suppress

import pytest
import zarr
import h5py
import pandas as pd

import anndata as ad
from anndata._io.specs.registry import IORegistryError
from anndata.compat import _clean_uns
from anndata._io.utils import (
    report_read_key_on_error,
)
from anndata.experimental import read_elem, write_elem
from anndata.tests.helpers import pytest_8_raises


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


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

        with pytest_8_raises(NotImplementedError, match=r"/X"):
            read_attr(group["X"])

        with pytest_8_raises(NotImplementedError, match=r"/group"):
            read_attr(group["group"])


def test_write_error_info(diskfmt, tmp_path):
    pth = tmp_path / f"failed_write.{diskfmt}"
    write = lambda x: getattr(x, f"write_{diskfmt}")(pth)

    # Assuming we don't define a writer for tuples
    a = ad.AnnData(uns={"a": {"b": {"c": (1, 2, 3)}}})

    with pytest_8_raises(IORegistryError, match=r"Error raised while writing key 'c'"):
        write(a)


def test_clean_uns():
    adata = ad.AnnData(
        uns=dict(species_categories=["a", "b"]),
        obs=pd.DataFrame({"species": [0, 1, 0]}, index=["a", "b", "c"]),
        var=pd.DataFrame({"species": [0, 1, 0, 2]}, index=["a", "b", "c", "d"]),
    )
    _clean_uns(adata)
    assert "species_categories" not in adata.uns
    assert isinstance(adata.obs["species"].dtype, pd.CategoricalDtype)
    assert adata.obs["species"].tolist() == ["a", "b", "a"]
    # var’s categories were overwritten by obs’s,
    # which we can detect here because var has too high codes
    assert pd.api.types.is_integer_dtype(adata.var["species"])


@pytest.mark.parametrize(
    "group_fn",
    [
        pytest.param(lambda _: zarr.group(), id="zarr"),
        pytest.param(lambda p: h5py.File(p / "test.h5", mode="a"), id="h5py"),
    ],
)
def test_only_child_key_reported_on_failure(tmp_path, group_fn):
    class Foo:
        pass

    group = group_fn(tmp_path)

    # This regex checks that the pattern inside the (?!...) group does not exist in the string
    # (?!...) is a negative lookahead
    # (?s) enables the dot to match newlines
    # https://stackoverflow.com/a/406408/130164 <- copilot suggested lol
    pattern = r"(?s)((?!Error raised while writing key '/?a').)*$"

    with pytest_8_raises(IORegistryError, match=pattern):
        write_elem(group, "/", {"a": {"b": Foo()}})

    write_elem(group, "/", {"a": {"b": [1, 2, 3]}})
    group["a/b"].attrs["encoding-type"] = "not a real encoding type"

    with pytest_8_raises(IORegistryError, match=pattern):
        read_elem(group)
