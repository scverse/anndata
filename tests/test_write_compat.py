"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING, cast, get_args

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr

import anndata as ad
from anndata._io.specs.registry import read_elem_partial
from anndata.acc import A, AdAcc, AdRef

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from zarr.storage import MemoryStore

    from anndata.experimental.backed import Dataset2D


def test_compares_by_features() -> None:
    assert ad.WriteConfig() == ad.WriteCompat.DEFAULT
    assert hash(ad.WriteConfig()) == hash(ad.WriteCompat.DEFAULT)
    assert ad.WriteConfig() <= ad.WriteCompat.V0_13
    assert not ad.WriteConfig() < ad.WriteCompat.V0_13
    assert ad.WriteCompat.V0_13 < ad.WriteCompat.V0_14
    assert ad.WriteConfig(accessor=True) < ad.WriteCompat.V0_14
    assert ad.WriteConfig(accessor=True, indexed_keys=True) == ad.WriteCompat.V0_14


def test_profiles_are_configs() -> None:
    assert isinstance(ad.WriteCompat.V0_13, ad.WriteConfig)
    assert not ad.WriteCompat.V0_13.accessor
    assert ad.WriteCompat.V0_14.accessor
    assert ad.WriteConfig().accessor == ad.WriteCompat.DEFAULT.accessor
    assert ad.WriteConfig(accessor=True).accessor


def test_literal_matches_enum() -> None:
    assert get_args(ad.WriteCompatStr.__value__) == tuple(
        str(c) for c in ad.WriteCompat
    )


def test_default_is_alias() -> None:
    assert ad.WriteCompat.DEFAULT in ad.WriteCompat
    assert "DEFAULT" not in {c.name for c in ad.WriteCompat}


def test_takes_str() -> None:
    assert ad.WriteCompat("0.13") is ad.WriteCompat.V0_13
    # this anndata version can’t write anything older than 0.13
    with pytest.raises(ValueError, match=r"0\.12"):
        ad.AnnData().unwriteable(compat="0.12")  # type: ignore[arg-type]


def test_accessor_needs_0_14(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """Accessors are refused before 0.14."""
    adata = ad.AnnData(shape=(3, 4), uns={"v": A.obsm["pca"][0]})
    assert adata.unwriteable()
    assert not adata.unwriteable(compat="0.14")
    assert not adata.unwriteable(compat=ad.WriteConfig(accessor=True))
    with pytest.raises(ValueError, match=r'Pass `compat="0.14"`'):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)


def test_accessor_subclass_unwriteable() -> None:
    """Subclasses from packages without an `anndata.acc` entry point can’t be written."""

    class MyRef(AdRef): ...

    adata = ad.AnnData(shape=(3, 4), uns={"r": AdAcc(ref_class=MyRef).obsm["pca"][0]})
    assert adata.unwriteable(compat="0.14")


def test_accessor_without_json_unwriteable(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """`RefAcc`s reach the accessor writer, but ones without JSON form (`A.obs`) fail."""
    adata = ad.AnnData(shape=(3, 4), uns={"r": A.obs})
    assert adata.unwriteable(compat="0.14")
    with pytest.raises(TypeError, match=r"Unsupported accessor"):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")


@pytest.mark.parametrize(
    "key",
    [
        pytest.param("a/b", id="slash"),
        pytest.param("CON", id="windows-device"),
        pytest.param("__x", id="zarr-reserved"),
        pytest.param("item_0", id="looks-like-child"),
    ],
)
def test_indexed_keys_roundtrip(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"], key: str
) -> None:
    """`compat="0.14"` allows any key everywhere keys can occur."""
    index = pd.Index(["c1", "c2"], name=f"idx{key}")
    adata = ad.AnnData(
        np.zeros((2, 3)),
        obs=pd.DataFrame({key: [1, 2]}, index=index),
        uns={key: 1},
        obsm={key: np.zeros((2, 2))},
        layers={key: np.zeros((2, 3))},
    )
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")
    back = getattr(ad, f"read_{diskfmt}")(diskfmt_store)

    assert list(back.obs.columns) == [key]
    assert back.obs.index.name == f"idx{key}"
    assert list(back.uns) == [key]
    assert list(back.obsm) == [key]
    assert set(back.layers) == {None, key}


def test_indexed_layout_on_disk(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """Names are stored in an array, children are numbered, dataframes follow the convention."""
    from anndata._io.utils import DF_CONVENTION

    obs = pd.DataFrame({"a/b": [1, 2], "x": [3, 4]}, index=["c1", "c2"])
    adata = ad.AnnData(np.zeros((2, 3)), obs=obs, uns={"c:d": 1})
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")

    def check(f: h5py.Group | zarr.Group) -> None:
        uns, obs_g = (cast("h5py.Group | zarr.Group", f[k]) for k in ("uns", "obs"))
        names = {
            k: cast("h5py.Dataset | zarr.Array", g[k])
            for g, k in [(uns, "keys"), (obs_g, "columns")]
        }

        assert dict(uns.attrs) == {"encoding-type": "dict", "encoding-version": "0.2.0"}
        assert set(uns) == {"keys", "item_0"}
        assert list(ad.io.read_elem(names["keys"])) == ["c:d"]  # type: ignore[arg-type]

        assert obs_g.attrs["encoding-version"] == "0.3.0"
        assert list(obs_g.attrs["df:shape"]) == [2, 3]  # type: ignore[arg-type]
        assert obs_g.attrs["df:index"] == 0
        assert set(obs_g) == {"columns", "column_0", "column_1", "column_2"}
        assert list(ad.io.read_elem(names["columns"])) == ["", "a/b", "x"]  # type: ignore[arg-type]
        if isinstance(f, zarr.Group):
            assert obs_g.attrs["zarr_conventions"] == [dict(DF_CONVENTION)]
        else:  # a zarr mechanism, and h5py can’t store it
            assert "zarr_conventions" not in obs_g.attrs

    if diskfmt == "zarr":
        check(zarr.open_group(diskfmt_store))
    else:
        with h5py.File(diskfmt_store) as f:
            check(f)


def test_column_named_index(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """With `compat="0.14"`, the index is found by position, so `_index` is a normal name."""
    df = pd.DataFrame({"_index": [1, 2]}, index=["c1", "c2"])
    adata = ad.AnnData(np.zeros((2, 3)), obs=df)
    with pytest.raises(ValueError, match=r"'_index' is a reserved name"):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)

    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")
    back = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    pd.testing.assert_frame_equal(back.obs, df)


def test_index_named_index(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """An unnamed index is stored as `""`, so `_index` keeps its name."""
    df = pd.DataFrame({"a": [1, 2]}, index=pd.Index(["c1", "c2"], name="_index"))
    adata = ad.AnnData(np.zeros((2, 3)), obs=df)
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")
    back = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert back.obs.index.name == "_index"


def test_dataframe_without_index() -> None:
    """The convention makes `df:index` optional, so we read a range index then."""
    g = zarr.open_group(zarr.storage.MemoryStore())
    ad.io.write_elem(g, "column_0", np.array([1, 2, 3]))
    ad.io.write_elem(g, "columns", np.array(["a"], dtype=object))
    g.attrs.update({
        "encoding-type": "dataframe",
        "encoding-version": "0.3.0",
        "df:shape": [3, 1],
    })
    expected = pd.DataFrame({"a": [1, 2, 3]})

    pd.testing.assert_frame_equal(cast("pd.DataFrame", ad.io.read_elem(g)), expected)
    partial = read_elem_partial(g, indices=(slice(1, 3), slice(None)))
    pd.testing.assert_frame_equal(partial, expected.iloc[1:3])
    if find_spec("xarray"):
        lazy = cast("Dataset2D", ad.experimental.read_elem_lazy(g))
        assert list(lazy.columns) == ["a"]
        np.testing.assert_array_equal(lazy["a"], [1, 2, 3])
        assert list(lazy.index) == ["0", "1", "2"]


def test_no_indexed_keys_before_0_14(tmp_path: Path) -> None:
    """The default profile still writes keys as child names and still rejects `/`."""
    ad.AnnData(np.zeros((2, 3)), uns={"a:b": 1}).write_h5ad(tmp_path / "t.h5ad")
    with h5py.File(tmp_path / "t.h5ad") as f:
        assert f["uns"].attrs["encoding-version"] == "0.1.0"
        assert f["obs"].attrs["encoding-version"] == "0.2.0"
        assert "a:b" in f["uns"]

    adata = ad.AnnData(np.zeros((2, 3)), uns={"a/b": 1})
    with pytest.raises(ValueError, match=r'compat="0.14"'):
        adata.write_h5ad(tmp_path / "slash.h5ad")


@pytest.mark.parametrize("attr", ["uns", "obsm"])
def test_empty_key(
    diskfmt_store: Path | MemoryStore,
    diskfmt: Literal["h5ad", "zarr"],
    attr: Literal["uns", "obsm"],
) -> None:
    """`""` names the parent group itself, so it needs `compat="0.14"`."""
    adata = ad.AnnData(np.zeros((2, 3)))
    getattr(adata, attr)[""] = np.zeros((2, 1))
    write = getattr(adata, f"write_{diskfmt}")
    with pytest.raises(ValueError, match=r"Empty keys are not allowed"):
        write(diskfmt_store)

    write(diskfmt_store, compat="0.14")
    back = getattr(ad, f"read_{diskfmt}")(diskfmt_store)
    assert list(getattr(back, attr)) == [""]


def test_indexed_keys_read_lazy(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    pytest.importorskip("xarray")
    index = pd.Index(["c1", "c2"], name="cell/id")
    obs = pd.DataFrame({"a/b": [1, 2], "c:d": [3, 4]}, index=index)
    adata = ad.AnnData(
        np.zeros((2, 3)), obs=obs, uns={"a/b": 1}, obsm={"c:d": np.zeros((2, 2))}
    )
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")
    lazy = ad.experimental.read_lazy(diskfmt_store)
    assert list(lazy.uns) == ["a/b"]
    assert list(lazy.obsm) == ["c:d"]
    assert set(lazy.obs.columns) == {"a/b", "c:d"}
    assert lazy.obs.index.name == "cell/id"
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", lazy.to_memory().obs), cast("pd.DataFrame", adata.obs)
    )


@pytest.mark.parametrize("axis", ["obs", "var"])
def test_indexed_keys_concat_on_disk(
    tmp_path: Path, diskfmt: Literal["h5ad", "zarr"], axis: Literal["obs", "var"]
) -> None:
    """Indexed keys survive on-disk concatenation, both along and across the axis."""
    paths: list[Path] = []
    for i in range(2):
        adata = ad.AnnData(
            np.ones((2, 2)),
            obs=pd.DataFrame({"g/h": ["x", "y"]}),
            obsm={"a/b": np.ones((2, 1))},
            varm={"c:d": np.ones((2, 1))},
            layers={"e/f": np.ones((2, 2))},
        )
        setattr(adata, f"{axis}_names", [f"{i}-{j}" for j in range(2)])
        paths.append(tmp_path / f"{i}.{diskfmt}")
        getattr(adata, f"write_{diskfmt}")(paths[-1], compat="0.14")
    out = tmp_path / f"out.{diskfmt}"
    ad.experimental.concat_on_disk(paths, out, axis=axis, merge="same")
    back = getattr(ad, f"read_{diskfmt}")(out)
    assert list(back.obsm) == ["a/b"]
    assert list(back.varm) == ["c:d"]
    assert set(back.layers) == {None, "e/f"}
    assert list(back.obs.columns) == ["g/h"]


def test_indexed_keys_roundtrip_full(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """Every element type survives `compat="0.14"`, e.g. nested and empty mappings."""
    from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata

    adata = gen_adata((10, 8), **GEN_ADATA_NO_XARRAY_ARGS)
    adata.uns["empty"] = {}
    getattr(adata, f"write_{diskfmt}")(diskfmt_store, compat="0.14")
    assert_equal(getattr(ad, f"read_{diskfmt}")(diskfmt_store), adata)


@pytest.mark.parametrize("compat", ["0.13", "0.14"])
def test_non_str_key(tmp_path: Path, compat: Literal["0.13", "0.14"]) -> None:
    adata = ad.AnnData(np.zeros((2, 3)), uns={1: 1})  # type: ignore[dict-item]
    with pytest.raises(TypeError, match=r"invalid key\. Should be str"):
        adata.write_h5ad(tmp_path / "t.h5ad", compat=compat)
