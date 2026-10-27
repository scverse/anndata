from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from pandas.api.extensions import ExtensionArray
from scipy import sparse as sp

from anndata import read_h5ad
from anndata.acc import A
from anndata.compat import (
    CSArray,
    CSMatrix,
    CupyArray,
    CupyCSRMatrix,
    CupySparseMatrix,
    DaskArray,
    XDataArray,
    XDataset,
    XVariable,
    is_scipy_sparse,
)
from anndata.tests.helpers import (
    DASK_CAN_SPARRAY,
    as_cupy,
    as_dense_dask_array,
    as_sparse_dask_array,
    as_sparse_dask_matrix,
)
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal

    from anndata import AnnData
    from anndata.acc import AdRef, InMemoryArray, RefAcc
    from anndata.typing import AlignedArray


needs_dask = pytest.mark.skipif(not find_spec("dask"), reason="dask not installed.")


ND_PATHS: list[tuple[AdRef, Callable[[AnnData], object]]] = [
    (A.X[:, :], lambda ad: ad.X),
    (A.X[:, "gene-3"], lambda ad: ad[:, "gene-3"].X),
    (A.X["cell-5", :], lambda ad: ad["cell-5"].X),
    (A.layers["a"][:, :], lambda ad: ad.layers["a"]),
    (A.layers["a"][:, "gene-18"], lambda ad: ad[:, "gene-18"].layers["a"]),
    (A.layers["a"]["cell-77", :], lambda ad: ad["cell-77"].layers["a"]),
    (A.obsm["umap"][0], lambda ad: asarray(ad.obsm["umap"])[:, 0]),
    (A.obsm["umap"][1], lambda ad: asarray(ad.obsm["umap"])[:, 1]),
    (A.obsm["umap"][-1], lambda ad: asarray(ad.obsm["umap"])[:, -1]),
    (A.obsm["umap"][-2], lambda ad: asarray(ad.obsm["umap"])[:, -2]),
    (A.varp["cons"][:, :], lambda ad: asarray(ad.varp["cons"])),
    (A.varp["cons"]["gene-46", :], lambda ad: asarray(ad.varp["cons"])[46, :]),
    (A.varp["cons"][:, "gene-46"], lambda ad: asarray(ad.varp["cons"])[:, 46]),
]

DF_PATHS: list[tuple[AdRef, Callable[[AnnData], object]]] = [
    (A.obs["type"], lambda ad: ad.obs["type"]),
    (A.obs.index, lambda ad: ad.obs.index),
    # TODO: copy entries from above that can also hold dataframes?
]

DASK_TYPES = [
    # TODO: check inner type
    pytest.param(as_dense_dask_array, DaskArray, DaskArray, id="da.array[np.ndarray]"),
    pytest.param(
        *(as_sparse_dask_array, DaskArray, DaskArray),
        marks=pytest.mark.skipif(
            not DASK_CAN_SPARRAY, reason="Dask does not support sparrays"
        ),
        id="da.array[sp.csr_array]",
    ),
    pytest.param(
        as_sparse_dask_matrix, DaskArray, DaskArray, id="da.array[sp.csc_matrix]"
    ),
    pytest.param(
        *(lambda x: as_cupy(as_dense_dask_array(x)), DaskArray, DaskArray),
        marks=pytest.mark.gpu,
        id="da.array[cp.ndarray]",
    ),
    pytest.param(
        *(lambda x: as_cupy(as_sparse_dask_matrix(x)), DaskArray, DaskArray),
        marks=pytest.mark.gpu,
        id="da.array[cpx.csr_matrix]",
    ),
]

# (container, 1D type, 2D type)
ND_TYPES = [
    pytest.param(np.asarray, np.ndarray, np.ndarray, id="np.ndarray"),
    # TODO: return 1D sparse array instead of 1D ndarray?
    pytest.param(sp.csr_array, np.ndarray, sp.csr_array, id="sp.csr_array"),
    pytest.param(sp.csc_matrix, np.ndarray, sp.csc_matrix, id="sp.csc_matrix"),
    pytest.param(as_cupy, CupyArray, CupyArray, marks=pytest.mark.gpu, id="cp.ndarray"),
    pytest.param(
        *(lambda x: as_cupy(sp.csr_array(x)), CupyArray, CupyCSRMatrix),
        marks=pytest.mark.gpu,
        id="cpx.csr_matrix",
    ),
    *[
        pytest.param(*t.values, marks=[*t.marks, needs_dask], id=t.id)  # noqa: PD011
        for t in DASK_TYPES
    ],
]

# (container, column type, index type)
DF_TYPES = [
    pytest.param(pd.DataFrame, ExtensionArray, ExtensionArray, id="pd.DataFrame"),
    pytest.param(
        lambda df: XDataset.from_dataframe(df),  # noqa: PLW0108
        XVariable,
        ExtensionArray,
        marks=pytest.mark.skipif(
            not find_spec("xarray"), reason="xarray not installed."
        ),
        id="xr.Dataset",
    ),
]


def convert_ndarrays(
    adata: AnnData, array_conv: Callable[[InMemoryArray], InMemoryArray], /
) -> None:
    def conv(v: AlignedArray) -> InMemoryArray:
        assert isinstance(v, np.ndarray | CSMatrix | CSArray)
        if isinstance(v, CSMatrix | CSArray):
            v = v.toarray()
        return array_conv(v)

    for k, v in adata.layers.items():
        adata.layers[k] = conv(v)
    for k, v in adata.varm.items():
        adata.varm[k] = conv(v)
    for k, v in adata.obsm.items():
        adata.obsm[k] = conv(v)
    for k, v in adata.obsp.items():
        adata.obsp[k] = conv(v)
    for k, v in adata.varp.items():
        adata.varp[k] = conv(v)


def convert_dataframes(
    adata: AnnData, df_conv: Callable[[pd.DataFrame], pd.DataFrame | XDataset], /
) -> None:
    assert isinstance(adata.obs, pd.DataFrame)
    adata.obs = df_conv(adata.obs)
    assert isinstance(adata.var, pd.DataFrame)
    adata.var = df_conv(adata.var)


def _expected2np(expected: object, ad_ref: AdRef, /) -> np.ndarray:
    ndim = len(ad_ref.dims)
    match expected:
        case np.ndarray():
            return expected.flatten() if ndim == 1 else expected
        case _ if is_scipy_sparse(expected):
            return _expected2np(expected.toarray(), ad_ref)
        case pd.Series() | pd.Index() | XDataArray():
            return expected.to_numpy()
        case DaskArray():
            return _expected2np(expected.compute(), ad_ref)
        case CupyArray() | CupySparseMatrix():
            return _expected2np(expected.get(), ad_ref)
        case _:
            pytest.fail(f"unhandled expected type {type(expected)}")


@pytest.fixture(
    params=[
        pytest.param(
            (ad_ref, ad_expected, *typ.values, convert),  # noqa: PD011
            marks=typ.marks,
            id=f"{ad_ref}-{typ.id}",
        )
        for paths, types, convert in (
            (ND_PATHS, ND_TYPES, convert_ndarrays),
            (DF_PATHS, DF_TYPES, convert_dataframes),
        )
        for ad_ref, ad_expected in paths
        for typ in types
    ]
)
def get_test_params(
    adata: AnnData, request: pytest.FixtureRequest
) -> tuple[AnnData, AdRef, type[InMemoryArray], InMemoryArray]:
    ad_ref, ad_expected, convert_array, *types_expected, convert = request.param
    convert(adata, convert_array)
    if convert is convert_ndarrays:  # (1D type, 2D type)
        type_expected = types_expected[len(ad_ref.dims) > 1]
    else:  # (column type, index type)
        type_expected = types_expected[ad_ref.idx is None]
    return adata, ad_ref, type_expected, _expected2np(ad_expected(adata), ad_ref)


def test_get_values(
    get_test_params: tuple[AnnData, AdRef, type[InMemoryArray], InMemoryArray],
) -> None:
    adata, ad_ref, type_expected, expected = get_test_params

    vals = adata[ad_ref]

    assert isinstance(vals, type_expected)
    vals_np = asarray(vals)
    np.testing.assert_array_equal(vals_np, expected, strict=True)


@pytest.mark.parametrize("i", [2, 99, -3], ids=str)
@pytest.mark.parametrize("dim", ["obs", "var"], ids=str)
def test_get_multi_out_of_range(
    adata: AnnData, dim: Literal["obs", "var"], i: int
) -> None:
    """Out-of-range column indices raise instead of yielding an empty array."""
    adata.varm["pcs"] = np.zeros((adata.n_vars, 2))  # `obsm["umap"]` also has 2 columns
    ad_ref = getattr(A, f"{dim}m")["umap" if dim == "obs" else "pcs"][i]
    assert ad_ref not in adata
    msg = rf"Column index `{i}` is out of range for A\.{dim}m\[.+\] with 2 columns\."
    with pytest.raises(IndexError, match=msg):
        adata[ad_ref]


@pytest.mark.parametrize(
    "ad_ref",
    [
        A.X[:, "gene-3"],
        A.X["cell-5", :],
        A.layers["a"][:, "gene-18"],
        A.obsm["umap"][1],
        A.obsm["umap"][-1],
        A.varp["cons"]["gene-46", :],
        A.varp["cons"][:, "gene-46"],
        A.obs["type"],
        A.obs["new"],
        A.var.index,
    ],
    ids=str,
)
@pytest.mark.parametrize("view", [False, True], ids=["actual", "view"])
def test_set_values(adata: AnnData, ad_ref: AdRef, *, view: bool) -> None:
    if view:
        adata = adata[:, :]
    (dim,) = ad_ref.dims
    value = np.arange(adata.n_obs if dim == "obs" else adata.n_vars)
    value = value.astype(str) if ad_ref.idx is None else value.astype(np.float32)

    adata[ad_ref] = value

    assert not adata.is_view
    np.testing.assert_array_equal(asarray(adata[ad_ref]), value)


@pytest.mark.parametrize(
    ("acc", "attr"),
    [
        (A.X, lambda ad: ad.X),
        (A.layers["new"], lambda ad: ad.layers["new"]),
        (A.obsm["new"], lambda ad: ad.obsm["new"]),
        (A.obsp["new"], lambda ad: ad.obsp["new"]),
    ],
    ids=str,
)
def test_set_full(
    adata: AnnData, acc: RefAcc, attr: Callable[[AnnData], object]
) -> None:
    n_cols = adata.n_obs if acc == A.obsp["new"] else adata.n_vars
    value = np.ones((adata.n_obs, n_cols), dtype=np.float32)
    adata[acc] = value
    assert attr(adata) is value


def test_set_full_df(adata: AnnData) -> None:
    df = pd.DataFrame(dict(x=np.arange(adata.n_obs)), index=adata.obs_names)
    adata[A.obs] = df
    obs = adata[A.obs]
    assert isinstance(obs, pd.DataFrame)
    pd.testing.assert_frame_equal(obs, df)


@pytest.mark.parametrize(
    ("index", "exc"), [(A.obsm, IndexError), (0, TypeError)], ids=str
)
def test_set_unsupported(adata: AnnData, index: object, exc: type[Exception]) -> None:
    with pytest.raises(exc):
        adata[index] = 0  # type: ignore[index]


def test_set_unsettable(adata: AnnData, tmp_path: Path) -> None:
    adata.X = adata.layers["a"]
    adata.write_h5ad(tmp_path / "test.h5ad")
    adata = read_h5ad(tmp_path / "test.h5ad", backed="r")
    with pytest.raises(TypeError, match=r"Cannot assign into _?CSRDataset"):
        adata[A.X[:, "gene-3"]] = np.zeros(adata.n_obs)


def test_subclass_without_set_warns() -> None:
    from anndata.acc import RefAcc

    with pytest.warns(FutureWarning, match=r"does not implement `set`"):

        class _NoSet(RefAcc):  # pyright: ignore[reportUnusedClass]
            pass


def test_set_multi_out_of_range(adata: AnnData) -> None:
    with pytest.raises(IndexError, match=r"out of range"):
        adata[A.obsm["umap"][2]] = np.zeros(adata.n_obs)
