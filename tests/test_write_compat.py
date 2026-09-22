"""Tests for :class:`anndata.WriteCompat` and the `compat` write argument."""

from __future__ import annotations

from typing import TYPE_CHECKING, get_args

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
