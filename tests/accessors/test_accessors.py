from __future__ import annotations

import importlib
import importlib.resources
import json
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING

import jsonschema
import pandas as pd
import pytest

from anndata import AnnData
from anndata.acc import A

if TYPE_CHECKING:
    from collections.abc import Callable, Collection
    from typing import Literal

    from anndata.acc import AdRef


type AdRefSer = Sequence[str | int | None]

with importlib.resources.open_text("anndata.acc", "acc-schema-v1.json") as f:
    SCHEMA = json.load(f)

PATHS: list[tuple[AdRef, AdRefSer]] = [
    (A.X[:, :], ["layers", None, None, None]),
    (A.X[:, "gene-3"], ["layers", None, None, "gene-3"]),
    (A.X["cell-5", :], ["layers", None, "cell-5", None]),
    (A.obs["type"], ["obs", "type"]),
    (A.obs.index, ["obs", None]),
    (A.layers["a"][:, :], ["layers", "a", None, None]),
    (A.layers["a"][:, "gene-18"], ["layers", "a", None, "gene-18"]),
    (A.layers["a"]["cell-77", :], ["layers", "a", "cell-77", None]),
    (A.obsm["umap"][0], ["obsm", "umap", 0]),
    (A.obsm["umap"][1], ["obsm", "umap", 1]),
    (A.varp["cons"]["gene-46", :], ["varp", "cons", "gene-46", None]),
    (A.varp["cons"][:, "gene-46"], ["varp", "cons", None, "gene-46"]),
]


@pytest.fixture(scope="session", params=PATHS, ids=[str(p[0]) for p in PATHS])
def path_and_expected_fn(
    request: pytest.FixtureRequest,
) -> tuple[AdRef, AdRefSer]:
    return request.param


@pytest.fixture(scope="session")
def ad_ref(path_and_expected_fn: tuple[AdRef, AdRefSer]) -> AdRef:
    return path_and_expected_fn[0]


@pytest.fixture(scope="session")
def ad_serialized(
    path_and_expected_fn: tuple[AdRef, AdRefSer],
) -> AdRefSer:
    return path_and_expected_fn[1]


def test_repr(ad_ref: AdRef) -> None:
    from anndata.acc import A  # for eval

    assert repr(ad_ref) == str(ad_ref)
    assert repr(ad_ref)[:2] == "A."
    assert eval(repr(ad_ref)) == ad_ref
    del A


@pytest.mark.parametrize("compare", ["json", "obj"])
def test_serialization(
    ad_ref: AdRef, ad_serialized: AdRefSer, compare: Literal["json", "obj"]
) -> None:
    if compare == "obj":
        assert A.from_json(ad_serialized) == ad_ref
    else:
        assert A.to_json(ad_ref) == ad_serialized


def test_serialization_schema(ad_serialized: AdRefSer) -> None:
    jsonschema.validate(ad_serialized, SCHEMA)


@pytest.mark.parametrize(
    ("ad_ref", "dims"),
    [
        pytest.param(A.X[:, :], {"obs", "var"}, id="x"),
        pytest.param(A.layers["y"][:, :], {"obs", "var"}, id="layer"),
        # selecting one obs gives a vector along the var dimension:
        pytest.param(A.layers["y"]["c", :], {"var"}, id="layer-obs"),
        pytest.param(A.var["a"], {"var"}, id="var"),
        pytest.param(A.obsm["c"][:, 0], {"obs"}, id="obsm"),
        pytest.param(A.varp["d"][:, :], ("var", "var"), id="varp"),
        pytest.param(A.varp["d"][:, "c2"], {"var"}, id="varp-col"),
    ],
)
def test_dims(ad_ref: AdRef, dims: Collection[Literal["obs", "var"]]) -> None:
    assert ad_ref.dims == dims


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
        pytest.param(lambda: A.X["c1"], id="x-notuple"),
        pytest.param(lambda: A.X[:3, :], id="x-partslice"),
        pytest.param(lambda: A.X[:, b""], id="x-nostr"),
        pytest.param(lambda: A.X[["a"], ["b"]], id="x-twolists"),
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
        pytest.param(
            lambda l: A.X[l("a", "b"), :], [A.X["a", :], A.X["b", :]], id="x-obs"
        ),
        pytest.param(
            lambda l: A.X[:, l("x", "y")], [A.X[:, "x"], A.X[:, "y"]], id="x-var"
        ),
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
    "ad_ref",
    [
        *(p[0] for p in PATHS),
        A.X,
        A.obs,
        A.layers,
        A.layers["a"],
        A.obsm,
        A.obsm["umap"],
        A.varp,
        A.varp["cons"],
    ],
    ids=str,
)
def test_in(adata: AnnData, ad_ref: AdRef) -> None:
    assert ad_ref in adata


@pytest.mark.parametrize(
    "ad_ref",
    [
        A.layers["a"]["gene-0", :],  # not an obs name
        A.layers["a"][:, "cell-3"],  # not a var name
        A.layers["b"],
        A.obsm["umap"][:, 3],
        A.obsm["b"],
        A.varm,
        A.obsp,
        A.varp["cons"]["cell-1", :],  # not a var name
        A.varp["cons"][:, "cell-2"],  # not a var name either
        A.varp["b"],
    ],
    ids=str,
)
def test_not_in(adata: AnnData, ad_ref: AdRef) -> None:
    assert ad_ref not in adata


def test_not_in_empty(ad_ref: AdRef) -> None:
    if ad_ref in {A.obs.index, A.var.index}:
        pytest.xfail("obs.index and var.index are always in AnnData")
    assert ad_ref not in AnnData()
