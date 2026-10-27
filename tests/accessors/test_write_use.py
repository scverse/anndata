from __future__ import annotations

from typing import TYPE_CHECKING, cast, overload

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse

import anndata as ad
from anndata.acc import A
from anndata.tests.helpers import assert_equal

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal

    from zarr.storage import MemoryStore

    from anndata.acc import AdRef, MapAcc, RefAcc


def _write(
    adata: ad.AnnData,
    store: Path | MemoryStore,
    diskfmt: Literal["h5ad", "zarr"],
    **kwargs,
) -> None:
    getattr(adata, f"write_{diskfmt}")(store, **kwargs)


def _read(store: Path | MemoryStore, diskfmt: Literal["h5ad", "zarr"]) -> ad.AnnData:
    return getattr(ad.io, f"read_{diskfmt}")(store)


@pytest.fixture
def adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    obs_names = [f"cell-{i}" for i in range(5)]
    var_names = [f"gene-{i}" for i in range(3)]
    return ad.AnnData(
        X=rng.random((5, 3)),
        obs=pd.DataFrame({"a": np.arange(5), "b": list("abcde")}, index=obs_names),
        var=pd.DataFrame({"c": np.arange(3.0)}, index=var_names),
        layers={"l": rng.random((5, 3))},
        obsm={"o": rng.random((5, 2))},
        obsp={"p": sparse.random(5, 5, format="csr", rng=rng)},
    )


@pytest.mark.parametrize("col", ["new", "a"], ids=["append", "overwrite"])
def test_write_obs_column(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    col: str,
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    adata.obs[col] = np.arange(5) * 10
    adata.obs["unwritten"] = 1
    _write(adata, diskfmt_store, diskfmt, use=A.obs[col])

    read = _read(diskfmt_store, diskfmt)
    expected_cols = ["a", "b"] if col == "a" else ["a", "b", col]
    assert list(read.obs.columns) == expected_cols
    assert_equal(read.obs, adata.obs[expected_cols])
    # unchanged otherwise
    read.obs["unwritten"] = cast("pd.DataFrame", adata.obs)["unwritten"]
    assert_equal(read, adata)


@overload
def _setter(
    elem: Literal["X", "obs", "obsm", "obsp", "layers"], sub_elem: str
) -> Callable[[ad.AnnData, ad.AnnData], None]: ...
@overload
def _setter(
    elem: Literal["X", "obs", "obsm", "obsp"],
    sub_elem: None = None,
) -> Callable[[ad.AnnData, ad.AnnData], None]: ...
def _setter(
    elem: Literal["X", "obs", "obsm", "obsp", "layers"], sub_elem: str | None = None
) -> Callable[[ad.AnnData, ad.AnnData], None]:
    def _set(dst: ad.AnnData, src: ad.AnnData) -> None:
        if sub_elem is not None:
            getattr(dst, elem).update({sub_elem: getattr(src, elem)[sub_elem]})
        else:
            setattr(dst, elem, getattr(src, elem))

    return _set


def _set_layers(dst: ad.AnnData, src: ad.AnnData) -> None:
    for k, v in src.layers.items():
        if k is not None:
            dst.layers[k] = v


@pytest.mark.parametrize(
    ("use", "apply"),
    [
        pytest.param(A.X, _setter("X"), id="X"),
        pytest.param(A.obs, _setter("obs"), id="obs"),
        pytest.param(
            A.layers["l"],
            _setter("layers", "l"),
            id="layer",
        ),
        pytest.param(
            A.layers["new"],
            _setter("layers", "new"),
            id="new-layer",
        ),
        pytest.param(A.layers, _set_layers, id="layers"),
        pytest.param(
            A.obsm["o"],
            _setter("obsm", "o"),
            id="obsm-key",
        ),
        pytest.param(A.obsm, _setter("obsm"), id="obsm"),
        pytest.param(
            A.obsp["p"],
            _setter("obsp", "p"),
            id="obsp-key",
        ),
        pytest.param(A.obsp, _setter("obsp"), id="obsp"),
    ],
)
def test_write_elem(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    use: RefAcc | MapAcc,
    apply: Callable[[ad.AnnData, ad.AnnData], None],
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    rng = np.random.default_rng(1)
    modified = adata.copy()
    modified.X = rng.random((5, 3))
    modified.obs["a"] = np.arange(5) + 1
    modified.layers["l"] = rng.random((5, 3))
    modified.layers["new"] = rng.random((5, 3))
    modified.obsm["o"] = rng.random((5, 2))
    modified.obsp["p"] = sparse.random(5, 5, format="csr", rng=rng)
    _write(modified, diskfmt_store, diskfmt, use=use)

    expected = adata.copy()
    apply(expected, modified)
    assert_equal(_read(diskfmt_store, diskfmt), expected)


@pytest.mark.parametrize(
    ("use", "match"),
    [
        pytest.param(A.obsm["o"][:, 0], r"only `obs`/`var` columns", id="vector"),
        pytest.param(A.obs.index, r"only `obs`/`var` columns", id="index"),
        pytest.param(A.obs["_index"], r"name of the obs index", id="reserved"),
    ],
)
def test_write_invalid_ref(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    use: AdRef,
    match: str,
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    with pytest.raises(ValueError, match=match):
        _write(adata, diskfmt_store, diskfmt, use=use)


def _rename_var(adata: ad.AnnData) -> ad.AnnData:
    adata.var_names = [f"other-{i}" for i in range(adata.n_vars)]
    return adata


@pytest.mark.parametrize(
    ("modify", "dim"),
    [
        pytest.param(lambda a: a[:2].copy(), "obs", id="obs-subset"),
        pytest.param(lambda a: a[::-1].copy(), "obs", id="obs-reordered"),
        pytest.param(_rename_var, "var", id="var-renamed"),
    ],
)
@pytest.mark.parametrize("use", [A.obs["a"], A.X], ids=["column", "X"])
def test_write_names_mismatch(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
    modify: Callable[[ad.AnnData], ad.AnnData],
    dim: Literal["obs", "var"],
    use: AdRef | RefAcc,
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    with pytest.raises(ValueError, match=rf"{dim}_names in the store do not match"):
        _write(modify(adata), diskfmt_store, diskfmt, use=use)


def test_write_not_anndata_store(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
) -> None:
    if diskfmt == "h5ad":
        h5py.File(diskfmt_store, "w").close()
    else:
        zarr.open_group(diskfmt_store, mode="w")
    with pytest.raises(ValueError, match=r"not an AnnData store"):
        _write(adata, diskfmt_store, diskfmt, use=A.X)


def test_write_list(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    rng = np.random.default_rng(1)
    adata.X = rng.random((5, 3))
    adata.obs["new"] = np.arange(5) * 10
    adata.obs["a"] = np.arange(5) + 1
    adata.obsm["o"] = rng.random((5, 2))
    _write(adata, diskfmt_store, diskfmt, use=[A.X, *A.obs[["new", "a"]], A.obsm])
    assert_equal(_read(diskfmt_store, diskfmt), adata)


def test_write_list_invalid_writes_nothing(
    adata: ad.AnnData,
    diskfmt: Literal["h5ad", "zarr"],
    diskfmt_store: Path | MemoryStore,
) -> None:
    _write(adata, diskfmt_store, diskfmt)
    expected = adata.copy()
    adata.X = np.zeros((5, 3))
    with pytest.raises(ValueError, match=r"only `obs`/`var` columns"):
        _write(adata, diskfmt_store, diskfmt, use=[A.X, A.obsm["o"][:, 0]])
    assert_equal(_read(diskfmt_store, diskfmt), expected)
