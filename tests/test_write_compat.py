"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from typing import TYPE_CHECKING, get_args

import h5py
import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import anndata as ad
from anndata.acc import A

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from zarr.storage import MemoryStore


def test_write_compat_compares_as_versions() -> None:
    assert Version("0.9") < ad.WriteCompat.V0_13
    assert Version("0.13") == ad.WriteCompat.V0_13


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


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(
            [1, "b", pd.DataFrame({"a": [1, 2]}, index=["x", "y"])], id="hetero"
        ),
        pytest.param([[1, 2], [3]], id="ragged"),
        pytest.param({"deep": [{"a": 1}, [2, "three"], None]}, id="nested"),
        pytest.param(A.obsm["pca"][0], id="accessor"),
    ],
)
def test_new_encodings_need_0_14(
    diskfmt_store: Path | MemoryStore,
    diskfmt: Literal["h5ad", "zarr"],
    value: object,
) -> None:
    """Values needing the `sequence`/`accessor` encodings are refused before 0.14."""
    adata = ad.AnnData(shape=(3, 4), uns={"v": value})
    assert adata.unwriteable()
    assert not adata.unwriteable(compat="0.14")
    with pytest.raises(ValueError, match=r'Pass `compat="0.14"`'):
        getattr(adata, f"write_{diskfmt}")(diskfmt_store)


def test_mixed_scalars_stringified_before_0_14(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]
) -> None:
    """Pre-0.14 numpy’s stringification of `[1, "b"]` is kept for compatibility."""
    adata = ad.AnnData(shape=(3, 4), uns={"v": [1, "b"]})
    getattr(adata, f"write_{diskfmt}")(diskfmt_store)
    assert getattr(ad, f"read_{diskfmt}")(diskfmt_store).uns["v"].tolist() == ["1", "b"]


@pytest.mark.parametrize(
    "value",
    [
        pytest.param([1, "b"], id="mixed-scalars"),
        pytest.param(["a", b"b"], id="mixed-str-bytes"),
    ],
)
def test_lossy_coercion_kept_before_0_14(value: list) -> None:
    """Sequences 0.14 refuses to mangle are still turned into arrays before it."""
    from anndata._io.specs.registry import sequence_is_arrayable

    assert sequence_is_arrayable(value, compat=ad.WriteCompat.V0_13)
    assert not sequence_is_arrayable(value, compat=ad.WriteCompat.V0_14)


@pytest.mark.parametrize(
    "key",
    [
        pytest.param("a/b", id="unsafe-char"),
        pytest.param("CON", id="reserved"),
        pytest.param("a%2Fb", id="looks-escaped-but-is-literal"),
    ],
)
def test_key_escaping_roundtrip(
    diskfmt_store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"], key: str
) -> None:
    """`compat="0.14"` escapes keys everywhere they can occur.

    One key per branch of `escape_key` is enough here – which keys need what is
    covered by `test_escape_key` and `test_escape_key_is_usable_and_reversible`.
    """
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


def test_reserved_keys_need_escaping_even_without_unsafe_chars(
    tmp_path: Path,
) -> None:
    """A group of otherwise-safe but reserved keys still switches to escaped mode."""
    adata = ad.AnnData(np.zeros((2, 2)), uns={"CON": 1, "plain": 2})
    adata.write_h5ad(tmp_path / "t.h5ad", compat="0.14")
    with h5py.File(tmp_path / "t.h5ad") as f:
        assert f["uns"].attrs["encoding-version"] == "0.2.0"
        assert sorted(f["uns"]) == ["%CON%", "plain"]
    assert set(ad.read_h5ad(tmp_path / "t.h5ad").uns) == {"CON", "plain"}


def test_escaped_keys_are_safe_on_disk(tmp_path: Path) -> None:
    """Child names carry no characters any mainstream file system forbids."""
    from anndata._io.utils import UNSAFE_KEY_CHARS

    adata = ad.AnnData(
        np.zeros((2, 3)),
        obs=pd.DataFrame({"a/b": [1, 2]}, index=["c1", "c2"]),
        uns={"c:d": 1, "100%": 2},
    )
    adata.write_h5ad(tmp_path / "t.h5ad", compat="0.14")
    with h5py.File(tmp_path / "t.h5ad") as f:
        assert sorted(f["uns"]) == ["100%25", "c%3Ad"]
        # `column-order`/`_index` name the child keys, so they are escaped too
        assert list(f["obs"].attrs["column-order"]) == ["a%2Fb"]
        assert "a%2Fb" in f["obs"]

        names: list[str] = []
        f.visit(names.append)
    # `%` introduces an escape sequence, everything else is gone
    unsafe = UNSAFE_KEY_CHARS - {"%"}
    assert not [n for n in names if unsafe & set(n.rsplit("/", 1)[-1])]


def test_version_bumped_only_where_needed(tmp_path: Path) -> None:
    """A group without special keys stays byte-compatible with older readers."""
    adata = ad.AnnData(
        np.zeros((2, 3)),
        obs=pd.DataFrame({"plain": [1, 2]}, index=["c1", "c2"]),
        uns={"plain": 1},
    )
    adata.write_h5ad(tmp_path / "t.h5ad", compat="0.14")
    with h5py.File(tmp_path / "t.h5ad") as f:
        assert f["uns"].attrs["encoding-version"] == "0.1.0"
        assert f["obs"].attrs["encoding-version"] == "0.2.0"


def test_no_escaping_before_0_14(tmp_path: Path) -> None:
    """The default profile still writes `:` literally and still rejects `/`."""
    ad.AnnData(np.zeros((2, 3)), uns={"a:b": 1}).write_h5ad(tmp_path / "t.h5ad")
    with h5py.File(tmp_path / "t.h5ad") as f:
        assert f["uns"].attrs["encoding-version"] == "0.1.0"
        assert "a:b" in f["uns"]

    adata = ad.AnnData(np.zeros((2, 3)), uns={"a/b": 1})
    with pytest.raises(ValueError, match=r'Pass `compat="0.14"`'):
        adata.write_h5ad(tmp_path / "slash.h5ad")
