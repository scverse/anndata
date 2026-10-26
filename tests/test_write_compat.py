"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from typing import TYPE_CHECKING, get_args

import pytest

import anndata as ad
from anndata.acc import A

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from zarr.storage import MemoryStore


def test_compares_by_features() -> None:
    assert ad.WriteConfig() == ad.WriteCompat.DEFAULT
    assert hash(ad.WriteConfig()) == hash(ad.WriteCompat.DEFAULT)
    assert ad.WriteConfig() <= ad.WriteCompat.V0_13
    assert not ad.WriteConfig() < ad.WriteCompat.V0_13
    assert ad.WriteCompat.V0_13 < ad.WriteCompat.V0_14
    assert ad.WriteConfig(accessor=True) == ad.WriteCompat.V0_14


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
