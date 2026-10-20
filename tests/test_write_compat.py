"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from typing import get_args

import pytest
from packaging.version import Version

import anndata as ad


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
