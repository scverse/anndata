"""Tests for `.uns` accessors and accessors cast into nested types (`.AD`, `.DF`)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData
from anndata.acc import A, DFAcc, ItemAcc, UnsMapAcc

if TYPE_CHECKING:
    from pathlib import Path

    from anndata.acc import AdRef, MapAcc, RefAcc


@pytest.fixture
def adata(adata: AnnData) -> AnnData:
    inner = AnnData(
        np.arange(6.0).reshape(2, 3),
        obs=pd.DataFrame(dict(x=[1, 2]), index=["c0", "c1"]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    adata.uns.update(
        nested=inner,
        df=pd.DataFrame(dict(c=[1, 2, 3]), index=list("abc")),
        foo=dict(bar=dict(baz=2)),
    )
    return adata


@pytest.mark.parametrize(
    ("acc", "expected"),
    [
        pytest.param(A.uns, "A.uns", id="uns"),
        pytest.param(
            A.uns["foo"]["bar"]["baz"], "A.uns['foo']['bar']['baz']", id="item"
        ),
        pytest.param(
            A.uns["nested"].AD.obs["x"], "A.uns['nested'].AD.obs['x']", id="ad"
        ),
        pytest.param(
            A.uns["n"].AD.uns["m"].AD.X[:, "g"],
            "A.uns['n'].AD.uns['m'].AD.X[:, 'g']",
            id="ad-ad",
        ),
        pytest.param(A.uns["df"].DF["c"], "A.uns['df'].DF['c']", id="uns-df"),
        pytest.param(A.uns["df"].DF.index, "A.uns['df'].DF.index", id="uns-df-idx"),
        pytest.param(A.varm["df"].DF["a"], "A.varm['df'].DF['a']", id="varm-df"),
    ],
)
def test_repr(acc: MapAcc | RefAcc | AdRef, expected: str) -> None:
    assert repr(acc) == expected
    assert eval(expected) == acc


def test_types() -> None:
    assert isinstance(A.uns, UnsMapAcc)
    assert isinstance(A.uns["a"], ItemAcc)
    assert A.uns["a"]["b"].path == ("a", "b")
    assert isinstance(A.uns["a"].DF, DFAcc)
    assert isinstance(A.obsm["a"].DF, DFAcc)
    assert isinstance(A.obs, DFAcc)  # `MetaAcc` is a `DFAcc`
    assert A.uns["df"].DF[["a", "b"]] == [A.uns["df"].DF["a"], A.uns["df"].DF["b"]]
    # nesting makes accessors distinct
    assert A.uns["a"].AD.obs["x"] != A.obs["x"]


@pytest.mark.parametrize(
    ("spec", "expected"),
    [
        pytest.param("uns.a", A.uns["a"], id="top"),
        pytest.param("uns.a.b", A.uns["a"]["b"], id="nested"),
        pytest.param("uns.a.0", A.uns["a"]["0"], id="digits-are-keys"),
        pytest.param("uns.0", A.uns["0"], id="first-is-key"),
    ],
)
def test_resolve(spec: str, expected: ItemAcc) -> None:
    assert A.resolve(spec) == expected
    assert A.resolve(spec, vec=False) == expected


@pytest.mark.parametrize("spec", ["uns.", "uns.a..b", "uns.a."])
def test_resolve_empty_segment(spec: str) -> None:
    with pytest.raises(ValueError, match=r"Cannot parse uns accessor"):
        A.resolve(spec)


@pytest.mark.parametrize("k", [0, 1.0, True, None])
def test_item_invalid_key(k: object) -> None:
    with pytest.raises(TypeError, match=r"Unsupported uns key"):
        A.uns["a"][k]  # type: ignore[index]


@pytest.mark.parametrize(
    ("ref", "dims"),
    [
        pytest.param(A.varm["df"].DF["a"], {"var"}, id="varm-df"),
        pytest.param(A.uns["df"].DF["c"], set(), id="uns-df"),
        pytest.param(A.uns["nested"].AD.obs["x"], {"obs"}, id="ad"),
    ],
)
def test_dims(ref: AdRef, dims: set[str]) -> None:
    assert ref.dims == dims


@pytest.mark.parametrize(
    ("acc", "scope"),
    [
        pytest.param(A.obs, None, id="top"),
        pytest.param(A.varm["df"].DF, None, id="top-df"),
        pytest.param(A.uns["n"].AD.obs, A.uns["n"], id="ad"),
        pytest.param(A.uns["n"].AD.uns["df"].DF, A.uns["n"], id="ad-uns-df"),
        pytest.param(A.uns["n"].AD.obsm["df"].DF, A.uns["n"], id="ad-obsm-df"),
        pytest.param(A.uns["n"].AD.uns["m"].AD.X, A.uns["n"].AD.uns["m"], id="ad-ad"),
    ],
)
def test_scope(acc: RefAcc, scope: RefAcc | None) -> None:
    assert acc.scope == scope


@pytest.mark.parametrize(
    ("ref", "match"),
    [
        pytest.param(A.uns["foo"].AD.obs["x"], r"dict, not an AnnData", id="ad"),
        pytest.param(A.uns["nested"].DF["x"], r"AnnData, not a data frame", id="df"),
    ],
)
def test_cast_invalid(adata: AnnData, ref: AdRef, match: str) -> None:
    with pytest.raises(TypeError, match=match):
        adata[ref]


@pytest.mark.parametrize(
    "acc",
    [
        A.uns["foo"].AD.obs["x"],
        A.uns["foo"].AD.uns["x"],
        A.uns["foo"].AD.obsm,
        A.uns["nested"].DF["x"],
        A.uns["nested"].DF,
    ],
    ids=str,
)
def test_cast_invalid_not_in(adata: AnnData, acc: AdRef | RefAcc | MapAcc) -> None:
    """Nothing is in an element of the wrong type."""
    assert acc not in adata


def test_get(adata: AnnData) -> None:
    assert adata[A.uns["foo"]["bar"]["baz"]] == 2
    assert adata[A.uns["nested"]] is adata.uns["nested"]
    np.testing.assert_array_equal(adata[A.uns["nested"].AD.obs["x"]], [1, 2])
    np.testing.assert_array_equal(adata[A.uns["nested"].AD.X["c1", :]], [3, 4, 5])
    np.testing.assert_array_equal(adata[A.uns["df"].DF["c"]], [1, 2, 3])
    np.testing.assert_array_equal(adata[A.uns["df"].DF.index], list("abc"))
    np.testing.assert_array_equal(
        adata[A.varm["df"].DF["b"]], adata[A.varm["df"][:, 1]]
    )
    with pytest.raises(TypeError, match=r"can only be used to get the whole item"):
        A.uns["foo"].get(adata, 0)  # type: ignore[arg-type]


def test_set(adata: AnnData) -> None:
    inner: AnnData = adata.uns["nested"]
    adata[A.uns["new"]] = 5
    adata[A.uns["foo"]["bar"]["baz"]] = 7
    adata[A.uns["nested"].AD.obs["y"]] = [3, 4]
    adata[A.uns["nested"].AD.X["c0", :]] = 0
    adata[A.uns["df"].DF["d"]] = 0
    adata[A.varm["df"].DF["c"]] = 1
    assert adata.uns["new"] == 5
    assert adata.uns["foo"]["bar"] == dict(baz=7)
    assert inner.obs["y"].tolist() == [3, 4]
    np.testing.assert_array_equal(adata[A.uns["nested"].AD.X["c0", :]], [0, 0, 0])
    assert adata.uns["df"].columns.tolist() == ["c", "d"]
    assert adata[A.varm["df"].DF].columns.tolist() == ["a", "b", "c"]


@pytest.mark.parametrize(
    "acc",
    [
        A.uns,
        A.uns["foo"]["bar"]["baz"],
        A.uns["nested"].AD.obs["x"],
        A.uns["nested"].AD.X,
        A.uns["df"].DF["c"],
        A.varm["df"].DF["a"],
    ],
    ids=str,
)
def test_in(adata: AnnData, acc: MapAcc | RefAcc | AdRef) -> None:
    assert acc in adata


@pytest.mark.parametrize(
    "acc",
    [
        A.uns["missing"],
        A.uns["foo"]["nope"],
        A.uns["missing"].AD.obs["x"],
        A.uns["missing"].AD.obsm,
        A.uns["nested"].AD.obsm,
        A.uns["nested"].AD.obs["missing"],
        A.uns["df"].DF["missing"],
        A.varm["missing"].DF["a"],
    ],
    ids=str,
)
def test_not_in(adata: AnnData, acc: MapAcc | RefAcc | AdRef) -> None:
    assert acc not in adata


def test_dims_item_unsupported() -> None:
    with pytest.raises(TypeError, match=r"Cast it"):
        A.uns["a"].dims(None)  # type: ignore[arg-type]


@pytest.fixture
def written(adata: AnnData, tmp_path: Path) -> Path:
    """Write `adata`, then change things in memory, only one of which a test writes."""
    path = tmp_path / "a.h5ad"
    adata.write_h5ad(path)
    adata[A.uns["foo"]["bar"]["baz"]] = 7
    adata[A.uns["foo"]["bar"]["added"]] = 8
    adata[A.uns["foo"]["qux"]] = 9
    adata[A.uns["nested"].AD.obs["x"]] = [5, 6]
    adata[A.uns["nested"].AD.X] = np.ones((2, 3))
    adata[A.uns["df"].DF["c"]] = [4, 5, 6]
    adata[A.varm["df"].DF["a"]] = 0.0
    adata[A.varm["df"].DF["b"]] = 0.0
    # not in the store
    adata[A.uns["new"]] = adata[A.uns["nested"]].copy()
    adata[A.uns["new_df"]] = pd.DataFrame(dict(c=[1]))
    adata[A.uns["foo"]["new"]] = dict(x=1)
    return path


@pytest.mark.parametrize(
    ("use", "expected_foo"),
    [
        pytest.param(A.uns["foo"]["bar"]["baz"], dict(bar=dict(baz=7)), id="leaf"),
        pytest.param(
            A.uns["foo"]["bar"]["added"],
            dict(bar=dict(baz=2, added=8)),
            id="new-leaf",
        ),
    ],
)
def test_write_use(
    adata: AnnData, written: Path, use: ItemAcc, expected_foo: dict
) -> None:
    import anndata as ad

    adata.write_h5ad(written, use=use)
    assert ad.read_h5ad(written).uns["foo"] == expected_foo


@pytest.mark.parametrize(
    ("use", "group", "encoding"),
    [
        pytest.param(A.uns["foo"]["new"]["x"], "uns/foo/new", "dict", id="uns-dict"),
        pytest.param(
            A.uns["foo"]["bar"]["baz"]["x"], "uns/foo/bar/baz", "dict", id="uns-scalar"
        ),
        pytest.param(A.uns["new"].AD.obs["x"], "uns/new", "anndata", id="ad"),
        pytest.param(A.uns["new_df"].DF["c"], "uns/new_df", "dataframe", id="df"),
    ],
)
def test_write_use_missing(
    adata: AnnData, written: Path, use: AdRef | RefAcc, group: str, encoding: str
) -> None:
    """Only fixed attributes like `uns` get created, nothing that’s just in memory."""
    with pytest.raises(ValueError, match=rf"no '{encoding}' element at '{group}'"):
        adata.write_h5ad(written, use=use)


@pytest.mark.parametrize(
    ("use", "changed", "unchanged"),
    [
        pytest.param(
            A.uns["nested"].AD.obs["x"],
            [A.uns["nested"].AD.obs["x"]],
            [A.uns["nested"].AD.X],
            id="ad-col",
        ),
        pytest.param(
            A.uns["nested"].AD.X,
            [A.uns["nested"].AD.X],
            [A.uns["nested"].AD.obs["x"]],
            id="ad-x",
        ),
        pytest.param(
            A.varm["df"].DF["a"],
            [A.varm["df"].DF["a"]],
            [A.varm["df"].DF["b"]],
            id="df-col",
        ),
        pytest.param(
            A.varm["df"].DF,
            [A.varm["df"].DF["a"], A.varm["df"].DF["b"]],
            [],
            id="df-whole",
        ),
        pytest.param(A.uns["df"].DF["c"], [A.uns["df"].DF["c"]], [], id="uns-df-col"),
    ],
)
def test_write_use_cast(
    adata: AnnData,
    written: Path,
    use: AdRef | RefAcc,
    changed: list[AdRef],
    unchanged: list[AdRef],
) -> None:
    import anndata as ad

    before = ad.read_h5ad(written)
    adata.write_h5ad(written, use=use)
    after = ad.read_h5ad(written)
    for ref in changed:
        np.testing.assert_array_equal(after[ref], adata[ref])
    for ref in unchanged:
        np.testing.assert_array_equal(after[ref], before[ref])


def test_write_use_cast_index(adata: AnnData, written: Path) -> None:
    with pytest.raises(ValueError, match=r"name of the varm/df index"):
        adata.write_h5ad(written, use=A.varm["df"].DF["_index"])


def test_write_use_cast_names_mismatch(adata: AnnData, written: Path) -> None:
    adata.uns["nested"].obs_names = ["d0", "d1"]
    with pytest.raises(ValueError, match=r"obs_names in the store do not match"):
        adata.write_h5ad(written, use=A.uns["nested"].AD.obs["x"])
