from __future__ import annotations

import importlib
import sys
import types
from contextlib import AbstractContextManager, nullcontext
from importlib.metadata import version
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from packaging.version import Version

import anndata as ad
from anndata._io.specs.registry import IORegistryError, to_writeable
from anndata._io.utils import report_read_key_on_error
from anndata._io.zarr import fast_zarr_context
from anndata.compat import _clean_uns
from anndata.tests.helpers import jnp

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal


@pytest.mark.parametrize(
    "group_fn",
    [
        pytest.param(lambda _: zarr.group(), id="zarr"),
        pytest.param(lambda p: h5py.File(p / "test.h5", mode="a"), id="h5py"),
    ],
)
@pytest.mark.parametrize("nested", [True, False], ids=["nested", "root"])
def test_key_error(
    *, tmp_path, group_fn: Callable[[Path], zarr.Group | h5py.Group], nested: bool
):
    @report_read_key_on_error
    def read_attr(_):
        raise NotImplementedError()

    group = group_fn(tmp_path)
    with group if isinstance(group, AbstractContextManager) else nullcontext():
        if nested:
            group = group.create_group("nested")
            path = "/nested"
        else:
            path = "/"
        group["X"] = np.array([1, 2, 3])
        group.create_group("group")

        with pytest.raises(
            NotImplementedError, match=rf"reading key 'X'.*from {path}$"
        ):
            read_attr(group["X"])

        with pytest.raises(
            NotImplementedError, match=rf"reading key 'group'.*from {path}$"
        ):
            read_attr(group["group"])


def test_write_error_info(diskfmt, diskfmt_store):
    write = lambda x: getattr(x, f"write_{diskfmt}")(diskfmt_store)

    # Assuming we don't define a writer for tuples
    a = ad.AnnData(uns={"a": {"b": {"c": (1, 2, 3)}}})
    assert a.unwriteable()

    with pytest.raises(
        IORegistryError, match=r"Error raised while writing key 'c'.*to /uns/a/b"
    ):
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
    pattern = r"(?s)^((?!Error raised while writing key '/?a').)*$"

    with pytest.raises(IORegistryError, match=pattern):
        ad.io.write_elem(group, "/", {"a": {"b": Foo()}})

    ad.io.write_elem(group, "/", {"a": {"b": [1, 2, 3]}})
    group["a/b"].attrs["encoding-type"] = "not a real encoding type"

    with pytest.raises(IORegistryError, match=pattern):
        ad.io.read_elem(group)


@pytest.mark.parametrize("to", [np.array, pd.Series, lambda x: x])
def test_to_writeable_passthrough(to: Callable[[list], np.ndarray | pd.Series | list]):
    x = to([1, 2, 3])
    result = to_writeable(x)
    assert result is x


def get_gpu_devices() -> list:
    assert jnp is not None
    return [
        d
        for d in jnp.__array_namespace_info__().devices()
        if d is not None and d.device_kind == "GPU"
    ]


skip_if_no_jax_gpu = pytest.mark.skipif(
    jnp is None or not any(get_gpu_devices()),
    reason="No GPU available",
)


@skip_if_no_jax_gpu
@pytest.mark.array_api
def test_to_writeable_jax_gpu_array() -> None:
    assert jnp is not None
    x = jnp.asarray([1.0, 2.0], device=get_gpu_devices()[0])
    result = to_writeable(x)
    assert isinstance(result, np.ndarray)
    np.testing.assert_array_equal(result, np.array([1.0, 2.0]))


@pytest.mark.array_api
def test_to_writeable_does_not_recurse() -> None:
    assert jnp is not None
    x = {"a": jnp.array([1.0, 2.0])}
    result = to_writeable(x)
    # since dict is not supported, it should return unchanged
    assert result is x


has_fused = Version(version("zarr")) >= Version("3.3")


@pytest.mark.parametrize(
    ("mk_group", "use_zarrs", "set_fused", "expected"),
    [
        pytest.param(
            zarr.create_group,
            True,
            True,
            "Fused",
            id="local_group_triggers_fused_when_set",
            marks=pytest.mark.skipif(
                not has_fused,
                reason="No fused pipeline to set",
            ),
        ),
        pytest.param(
            zarr.create_group,
            True,
            False,
            "Zarrs",
            id="local_group_triggers_zarrs",
        ),
        # zarrs cannot handle memory stores, so it falls back to zarr-python
        pytest.param(
            lambda x: zarr.create_group(zarr.storage.MemoryStore()),
            False,
            False,
            "Fused" if has_fused else "Batched",
            id="mem_group_triggers_fused",
        ),
    ],
)
def test_zarr_context(
    monkeypatch,
    tmp_path,
    mk_group,
    expected: Literal["Fused", "Batched", "Zarrs"],
    *,
    use_zarrs: bool,
    set_fused: bool,
):
    if use_zarrs:
        fake_module = types.ModuleType("zarrs")
        fake_module.__spec__ = importlib.machinery.ModuleSpec("zarrs", loader=None)

        monkeypatch.setitem(sys.modules, "zarrs", fake_module)
    outer_ctx = (
        zarr.config.set({
            "codec_pipeline.path": "zarr.core.codec_pipeline.FusedCodecPipeline"
        })
        if set_fused
        else nullcontext()
    )
    with outer_ctx:
        with fast_zarr_context():
            g = mk_group(tmp_path)
            if not use_zarrs:
                g["foo"] = np.ones((2, 3))
                pipeline = str(g["foo"]._async_array.codec_pipeline)
            else:
                # we don't install zarrs, so just check that the context sets it
                # instead of fallback behavior
                pipeline = zarr.config.get("codec_pipeline.path")
            assert expected in pipeline
        # `fast_zarr_context` context works and does not leak state i.e., old pipeline returns
        if set_fused:
            assert "Fused" in zarr.config.get("codec_pipeline.path")
        else:
            assert "Batched" in zarr.config.get("codec_pipeline.path")
