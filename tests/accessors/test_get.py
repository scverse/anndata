from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from pandas.api.extensions import ExtensionArray
from scipy import sparse as sp

from anndata.acc import A
from anndata.compat import (
    CupyArray,
    CupyCSRMatrix,
    CupySparseMatrix,
    DaskArray,
    XDataArray,
    XDataset,
    XVariable,
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

    from anndata import AnnData
    from anndata.acc import AdRef, InMemoryArray


needs_dask = pytest.mark.skipif(not find_spec("dask"), reason="dask not installed.")


ND_PATHS: list[tuple[AdRef, Callable[[AnnData], InMemoryArray]]] = [
    (A.X[:, :], lambda ad: ad.X),
    (A.X[:, "gene-3"], lambda ad: ad[:, "gene-3"].X),
    (A.X["cell-5", :], lambda ad: ad["cell-5"].X),
    (A.layers["a"][:, :], lambda ad: ad.layers["a"]),
    (A.layers["a"][:, "gene-18"], lambda ad: ad[:, "gene-18"].layers["a"]),
    (A.layers["a"]["cell-77", :], lambda ad: ad["cell-77"].layers["a"]),
    (A.obsm["umap"][0], lambda ad: asarray(ad.obsm["umap"])[:, 0]),
    (A.obsm["umap"][1], lambda ad: asarray(ad.obsm["umap"])[:, 1]),
    (A.varp["cons"]["gene-46", :], lambda ad: asarray(ad.varp["cons"])[46, :]),
    (A.varp["cons"][:, "gene-46"], lambda ad: asarray(ad.varp["cons"])[:, 46]),
]

DF_PATHS: list[tuple[AdRef, Callable[[AnnData], InMemoryArray]]] = [
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
        pytest.param(*t.values, marks=[*t.marks, needs_dask], id=t.id)
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
    def conv(v: InMemoryArray) -> InMemoryArray:
        if isinstance(v, sp.sparray | sp.spmatrix):
            v = v.toarray()
        return array_conv(v)

    adata.X = conv(adata.X)
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
    adata: AnnData, df_conv: Callable[[pd.Series], InMemoryArray], /
) -> None:
    adata.obs = df_conv(adata.obs)
    adata.var = df_conv(adata.var)


def _expected2np(expected: InMemoryArray, ad_ref: AdRef, /) -> np.ndarray:
    ndim = len(ad_ref.dims)
    match expected:
        case np.ndarray():
            return expected.flatten() if ndim == 1 else expected
        case sp.sparray() | sp.spmatrix():
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
            (ad_ref, ad_expected, *typ.values, convert),
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
