from __future__ import annotations

import importlib
import importlib.resources
import json
from collections.abc import Callable, Sequence
from string import ascii_lowercase
from typing import TYPE_CHECKING

import jsonschema
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData
from anndata.acc import A

if TYPE_CHECKING:
    from collections.abc import Callable, Collection
    from typing import Literal

    from anndata.acc import AdRef


type AdRefExpected = Callable[[AnnData], np.ndarray | sp.coo_array | pd.Series]
type AdRefSer = Sequence[str | int | None]

with importlib.resources.open_text("anndata.acc", "acc-schema.json") as f:
    SCHEMA = json.load(f)

PATHS: list[tuple[AdRef, AdRefSer, AdRefExpected]] = [
    (A[:, :], ["layers", None, None, None], lambda ad: ad.X),
    (
        A[:, "gene-3"],
        ["layers", None, None, "gene-3"],
        lambda ad: ad[:, "gene-3"].X.flatten(),
    ),
    (
        A["cell-5", :],
        ["layers", None, "cell-5", None],
        lambda ad: ad["cell-5"].X.flatten(),
    ),
    (A.obs["type"], ["obs", "type"], lambda ad: ad.obs["type"]),
    (A.obs.index, ["obs", None], lambda ad: ad.obs.index.values),
    (
        A.layers["a"][:, :],
        ["layers", "a", None, None],
        lambda ad: ad.layers["a"].copy().toarray(),
    ),
    (
        A.layers["a"][:, "gene-18"],
        ["layers", "a", None, "gene-18"],
        lambda ad: ad[:, "gene-18"].layers["a"].copy().toarray().flatten(),
    ),
    (
        A.layers["a"]["cell-77", :],
        ["layers", "a", "cell-77", None],
        lambda ad: ad["cell-77"].layers["a"].copy().toarray().flatten(),
    ),
    (A.obsm["umap"][0], ["obsm", "umap", 0], lambda ad: ad.obsm["umap"][:, 0]),
    (A.obsm["umap"][1], ["obsm", "umap", 1], lambda ad: ad.obsm["umap"][:, 1]),
    (
        A.varp["cons"]["gene-46", :],
        ["varp", "cons", "gene-46", None],
        lambda ad: ad.varp["cons"][46, :].toarray(),
    ),
    (
        A.varp["cons"][:, "gene-46"],
        ["varp", "cons", None, "gene-46"],
        lambda ad: ad.varp["cons"][:, 46].toarray(),
    ),
]


@pytest.fixture(scope="session", params=PATHS, ids=[str(p[0]) for p in PATHS])
def path_and_expected_fn(
    request: pytest.FixtureRequest,
) -> tuple[AdRef, AdRefSer, AdRefExpected]:
    return request.param


@pytest.fixture(scope="session")
def ad_path(path_and_expected_fn: tuple[AdRef, AdRefSer, AdRefExpected]) -> AdRef:
    return path_and_expected_fn[0]


@pytest.fixture(scope="session")
def ad_serialized(
    path_and_expected_fn: tuple[AdRef, AdRefSer, AdRefExpected],
) -> AdRefSer:
    return path_and_expected_fn[1]


@pytest.fixture(scope="session")
def ad_expected(
    path_and_expected_fn: tuple[AdRef, AdRefSer, AdRefExpected],
) -> AdRefExpected:
    return path_and_expected_fn[2]


@pytest.fixture
def adata() -> AnnData:
    gen = np.random.default_rng()
    x = gen.random((100, 50), dtype=np.float32)
    layers = dict(a=sp.random(100, 50, format="csr"))
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
    varp = dict(cons=sp.csr_array(sp.random(50, 50)))
    return AnnData(x, obs, var, layers=layers, obsm=obsm, varm={}, obsp={}, varp=varp)


def test_repr(ad_path: AdRef) -> None:
    from anndata.acc import A  # for eval

    assert repr(ad_path) == str(ad_path)
    assert repr(ad_path)[:2] in {"A.", "A["}
    assert eval(repr(ad_path)) == ad_path
    del A


@pytest.mark.parametrize("compare", ["json", "obj"])
def test_serialization(
    ad_path: AdRef, ad_serialized: AdRefSer, compare: Literal["json", "obj"]
) -> None:
    if compare == "obj":
        assert A.from_json(ad_serialized) == ad_path
    else:
        assert A.to_json(ad_path) == ad_serialized


def test_serialization_schema(ad_serialized: AdRefSer) -> None:
    jsonschema.validate(ad_serialized, SCHEMA)


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
def test_axes(ad_path: AdRef, axes: Collection[Literal["obs", "var"]]) -> None:
    assert ad_path.axes == axes


@pytest.mark.parametrize(
    ("obj", "expected"),
    [
        pytest.param(A.layers["a"], "a", id="A.layers[]"),
        *[pytest.param(o, "obs", id=str(o)) for o in [A.obs, A.obsm, A.obsp]],
        *[pytest.param(v, "var", id=str(v)) for v in [A.var, A.varm, A.varp]],
        *[pytest.param(o["b"], ("obs", "b"), id=f"{o}[]") for o in [A.obsm, A.obsp]],
        *[pytest.param(v["c"], ("var", "c"), id=f"{v}[]") for v in [A.varm, A.varp]],
        pytest.param(A.obs["d"], (A.obs, "d"), id="path"),
    ],
)
def test_match(*, obj: object, expected: object) -> None:
    from anndata import acc

    match obj:
        case (
            acc.LayerAcc(arg)
            | acc.MetaAcc(arg)
            | acc.MultiMapAcc(arg)
            | acc.GraphMapAcc(arg)
        ):
            assert len(type(obj).__match_args__) == 1
            assert arg == expected
        case acc.MultiAcc(a0, a1) | acc.GraphAcc(a0, a1) | acc.AdRef(a0, a1):
            assert len(type(obj).__match_args__) == 2
            assert (a0, a1) == expected
        case _:
            pytest.fail(f"unhandled case: {obj}")


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
def test_invalid(mk_path: Callable[[], AdRef]) -> None:
    with pytest.raises((ValueError, TypeError)):
        mk_path()


@pytest.mark.parametrize(
    ("mk_expr", "expanded"),
    [
        pytest.param(lambda l: A.obs[l("a", "b")], [A.obs["a"], A.obs["b"]], id="obs"),
        pytest.param(lambda l: A.var[l("x", "y")], [A.var["x"], A.var["y"]], id="var"),
        pytest.param(lambda l: A[l("a", "b"), :], [A["a", :], A["b", :]], id="x-obs"),
        pytest.param(lambda l: A[:, l("x", "y")], [A[:, "x"], A[:, "y"]], id="x-var"),
        pytest.param(
            lambda l: A.layers["l"][l("a", "b"), :],
            [A.layers["l"]["a", :], A.layers["l"]["b", :]],
            id="layer-obs",
        ),
        pytest.param(
            lambda l: A.layers["l"][:, l("x", "y")],
            [A.layers["l"][:, "x"], A.layers["l"][:, "y"]],
            id="layer-var",
        ),
        pytest.param(
            lambda l: A.varm["p"][l(0, 1)],
            [A.varm["p"][:, 0], A.varm["p"][:, 1]],
            id="varm-int",
        ),
        pytest.param(
            lambda l: A.varm["p"][:, l(3, 4)],
            [A.varm["p"][:, 3], A.varm["p"][:, 4]],
            id="varm-tuple",
        ),
        pytest.param(
            lambda l: A.obsp["p"][l("a", "b"), :],
            [A.obsp["p"]["a", :], A.obsp["p"]["b", :]],
            id="obsp-0",
        ),
        pytest.param(
            lambda l: A.obsp["p"][:, l("a", "b")],
            [A.obsp["p"][:, "a"], A.obsp["p"][:, "b"]],
            id="obsp-1",
        ),
    ],
)
@pytest.mark.parametrize(
    "mk_list",
    [
        lambda *e: [*e],
        lambda *e: pd.Index([*e], dtype=type(e[0])),
    ],
    ids=["list", "index"],
)
def test_special[C](
    mk_expr: Callable[[C], list[AdRef]], mk_list: C, expanded: list[AdRef]
) -> None:
    """Some paths have shortcuts for convenience."""
    expr = mk_expr(mk_list)
    assert expr == expanded


@pytest.mark.parametrize(
    "ad_path",
    [
        *(p[0] for p in PATHS),
        A,
        A.obs,
        A.layers,
        A.layers["a"],
        A.varm,
        A.obsm["umap"],
        A.obsp,
        A.varp["cons"],
    ],
    ids=str,
)
def test_in(adata: AnnData, ad_path: AdRef) -> None:
    assert ad_path in adata


@pytest.mark.parametrize(
    "ad_path",
    [
        A.layers["a"]["gene-0", :],  # not an obs name
        A.layers["a"][:, "cell-3"],  # not a var name
        A.layers["b"],
        A.obsm["umap"][:, 3],
        A.obsm["b"],
        A.varp["cons"]["cell-1", :],  # not a var name
        A.varp["cons"][:, "cell-2"],  # not a var name either
        A.varp["b"],
    ],
    ids=str,
)
def test_not_in(adata: AnnData, ad_path: AdRef) -> None:
    assert ad_path not in adata


def test_not_in_empty(ad_path: AdRef) -> None:
    if ad_path in {A.obs.index, A.var.index}:
        pytest.xfail("obs.index and var.index are always in AnnData")
    assert ad_path not in AnnData()


def test_get_values(
    adata: AnnData,
    ad_path: AdRef,
    ad_expected: Callable[[AnnData], np.ndarray | sp.coo_array | pd.Series],
) -> None:
    vals = ad_path(adata)
    np.testing.assert_array_equal(vals, ad_expected(adata), strict=True)
