from __future__ import annotations

from collections.abc import Callable
from string import ascii_lowercase
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData
from anndata.acc import A

if TYPE_CHECKING:
    from collections.abc import Callable, Collection
    from typing import Literal

    from anndata.acc import AdPath


type AdPathExpected = Callable[[AnnData], np.ndarray | sp.coo_array | pd.Series]

PATHS: list[tuple[AdPath, AdPathExpected]] = [
    (A[:, :], lambda ad: ad.X),
    (A[:, "gene-3"], lambda ad: ad[:, "gene-3"].X.flatten()),
    (A["cell-5", :], lambda ad: ad["cell-5"].X.flatten()),
    (A.obs["type"], lambda ad: ad.obs["type"]),
    (A.obs.index, lambda ad: ad.obs.index.values),
    (A.layers["a"][:, :], lambda ad: ad.layers["a"].copy().toarray()),
    (
        A.layers["a"][:, "gene-18"],
        lambda ad: ad[:, "gene-18"].layers["a"].copy().toarray().flatten(),
    ),
    (
        A.layers["a"]["cell-77", :],
        lambda ad: ad["cell-77"].layers["a"].copy().toarray().flatten(),
    ),
    (A.obsm["umap"][0], lambda ad: ad.obsm["umap"][:, 0]),
    (A.obsm["umap"][1], lambda ad: ad.obsm["umap"][:, 1]),
    (A.varp["cons"]["gene-46", :], lambda ad: ad.varp["cons"][46, :].toarray()),
    (A.varp["cons"][:, "gene-46"], lambda ad: ad.varp["cons"][:, 46].toarray()),
]


@pytest.fixture(scope="session", params=PATHS, ids=[str(p[0]) for p in PATHS])
def path_and_expected_fn(
    request: pytest.FixtureRequest,
) -> tuple[AdPath, AdPathExpected]:
    return request.param


@pytest.fixture(scope="session")
def ad_path(path_and_expected_fn: tuple[AdPath, AdPathExpected]) -> AdPath:
    return path_and_expected_fn[0]


@pytest.fixture(scope="session")
def ad_expected(path_and_expected_fn: tuple[AdPath, AdPathExpected]) -> AdPathExpected:
    return path_and_expected_fn[1]


@pytest.fixture
def adata() -> AnnData:
    gen = np.random.default_rng()
    x = gen.random((100, 50), dtype=np.float32)
    layers = dict(a=sp.random(100, 50, rng=gen, format="csr"))
    obs = pd.DataFrame(
        dict(type=gen.integers(0, 3, size=100)),
        index="cell-" + pd.array(range(100)).astype(str),
    )
    var_grp = pd.Categorical(
        gen.integers(0, 6, size=50), categories=list(ascii_lowercase[:5])
    )
    var = pd.DataFrame(
        dict(grp=var_grp),
        index="gene-" + pd.array(range(50)).astype(str),
    )
    obsm = dict(umap=gen.random((100, 2)))
    varp = dict(cons=sp.csr_array(sp.random(50, 50, rng=gen)))
    return AnnData(x, obs, var, layers=layers, obsm=obsm, varm={}, obsp={}, varp=varp)


def test_repr(ad_path: AdPath) -> None:
    from anndata.acc import A  # for eval

    assert repr(ad_path) == str(ad_path)
    assert repr(ad_path)[:2] in {"A.", "A["}
    assert eval(repr(ad_path)) == ad_path
    del A


@pytest.mark.parametrize(
    ("ad_path", "axes"),
    [
        pytest.param(A[:, :], {"obs", "var"}, id="x"),
        pytest.param(A.layers["y"][:, :], {"obs", "var"}, id="layer"),
        # selecting one obs gives a vector along the var axis:
        pytest.param(A.layers["y"]["c", :], {"var"}, id="layer-obs"),
        pytest.param(A.var["a"], {"var"}, id="var"),
        pytest.param(A.obsm["c"][:, 0], {"obs"}, id="obsm"),
        pytest.param(A.varp["d"][:, :], ("var", "var"), id="varp"),
        pytest.param(A.varp["d"][:, "c2"], {"var"}, id="varp-col"),
    ],
)
def test_axes(ad_path: AdPath, axes: Collection[Literal["obs", "var"]]) -> None:
    assert ad_path.axes == axes


@pytest.mark.parametrize(
    ("spec", "expected"),
    [
        # TODO: pytest.param("???", A[:, :], id="x"),
        pytest.param("layers.y[c,:]", A.layers["y"]["c", :], id="layer-obs"),
        pytest.param("layers.y[:,g]", A.layers["y"][:, "g"], id="layer-var"),
        pytest.param("obs.a", A.obs["a"], id="obs"),
        pytest.param("var.b", A.var["b"], id="var"),
        pytest.param("obsm.c.0", A.obsm["c"][:, 0], id="obsm"),
        pytest.param("varm.d.1", A.varm["d"][:, 1], id="varm"),
        pytest.param("obsp.g[c1,:]", A.obsp["d"]["c1", :], id="obsp"),
        pytest.param("obsp.g[:,c2]", A.obsp["d"][:, "c2"], id="varp"),
    ],
)
def test_resolve(spec: str, expected: AdPath) -> None:
    try:
        assert A.resolve(spec) == expected
        assert spec == expected
        assert expected != ""
    except NotImplementedError:
        pytest.xfail("not implemented")


@pytest.mark.parametrize(
    "mk_path",
    [
        pytest.param(lambda: A[:3, :], id="x-partslice"),
        pytest.param(lambda: A[:, b""], id="x-nostr"),
        pytest.param(lambda: A[["a"], ["b"]], id="x-twolists"),
        pytest.param(lambda: A.layers[1], id="layers-nostr"),
        pytest.param(lambda: A.layers["a"][:3, :], id="layer-partslice"),
        pytest.param(lambda: A.layers["a"][:, b""], id="layer-nostr"),
        pytest.param(lambda: A.layers["a"][["a"], ["b"]], id="layer-twolists"),
        pytest.param(lambda: A.obs[1], id="obs-nostr"),
        pytest.param(lambda: A.obsm[0], id="obsms-nostr"),
        pytest.param(lambda: A.obsm["a"][:3, 0], id="obsm-partslice"),
        pytest.param(lambda: A.obsm["a"]["b"], id="obsm-noint"),
        pytest.param(lambda: A.varp[0], id="varps-nostr"),
        pytest.param(lambda: A.varp["x"][0, :], id="varp-nostr-inner"),
        pytest.param(lambda: A.varp["x"]["a", "b"], id="varp-twostr"),
        pytest.param(lambda: A.varp["x"][["a"], ["b"]], id="varp-twolists"),
    ],
)
def test_invalid(mk_path: Callable[[], AdPath]) -> None:
    with pytest.raises((ValueError, TypeError)):
        mk_path()


@pytest.mark.parametrize(
    ("expr", "expanded"),
    [
        pytest.param(A.obs[["a", "b"]], [A.obs["a"], A.obs["b"]], id="obs"),
        pytest.param(A.var[["x", "y"]], [A.var["x"], A.var["y"]], id="var"),
        pytest.param(A[["a", "b"], :], [A["a", :], A["b", :]], id="x-obs"),
        pytest.param(A[:, ["x", "y"]], [A[:, "x"], A[:, "y"]], id="x-var"),
        pytest.param(
            A.layers["l"][["a", "b"], :],
            [A.layers["l"]["a", :], A.layers["l"]["b", :]],
            id="layer-obs",
        ),
        pytest.param(
            A.layers["l"][:, ["x", "y"]],
            [A.layers["l"][:, "x"], A.layers["l"][:, "y"]],
            id="layer-var",
        ),
        pytest.param(
            A.varm["p"][[0, 1]],
            [A.varm["p"][:, 0], A.varm["p"][:, 1]],
            id="varm-int",
        ),
        pytest.param(
            A.varm["p"][:, [3, 4]],
            [A.varm["p"][:, 3], A.varm["p"][:, 4]],
            id="varm-tuple",
        ),
        pytest.param(
            A.obsp["p"][["a", "b"], :],
            [A.obsp["p"]["a", :], A.obsp["p"]["b", :]],
            id="obsp-0",
        ),
        pytest.param(
            A.obsp["p"][:, ["a", "b"]],
            [A.obsp["p"][:, "a"], A.obsp["p"][:, "b"]],
            id="obsp-1",
        ),
    ],
)
def test_special(expr: AdPath, expanded: object) -> None:
    """Some paths have shortcuts for convenience."""
    assert expr == expanded


def test_get_values(
    adata: AnnData,
    ad_path: AdPath,
    ad_expected: Callable[[AnnData], np.ndarray | sp.coo_array | pd.Series],
) -> None:
    vals = ad_path(adata)
    np.testing.assert_array_equal(vals, ad_expected(adata), strict=True)
