from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from pandas.api.extensions import ExtensionArray
from scipy import sparse as sp

from anndata.acc import A
from anndata.compat import XDataArray, XDataset, XVariable
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata import AnnData
    from anndata.acc import AdRef, Array


ND_PATHS: list[tuple[AdRef, Callable[[AnnData], Array]]] = [
    (A[:, :], lambda ad: ad.X),
    (A[:, "gene-3"], lambda ad: ad[:, "gene-3"].X),
    (A["cell-5", :], lambda ad: ad["cell-5"].X),
    (A.layers["a"][:, :], lambda ad: ad.layers["a"]),
    (A.layers["a"][:, "gene-18"], lambda ad: ad[:, "gene-18"].layers["a"]),
    (A.layers["a"]["cell-77", :], lambda ad: ad["cell-77"].layers["a"]),
    (A.obsm["umap"][0], lambda ad: ad.obsm["umap"][:, 0]),
    (A.obsm["umap"][1], lambda ad: ad.obsm["umap"][:, 1]),
    (A.varp["cons"]["gene-46", :], lambda ad: ad.varp["cons"][46, :]),
    (A.varp["cons"][:, "gene-46"], lambda ad: ad.varp["cons"][:, 46]),
]

DF_PATHS: list[tuple[AdRef, Callable[[AnnData], Array]]] = [
    (A.obs["type"], lambda ad: ad.obs["type"]),
    (A.obs.index, lambda ad: ad.obs.index.values),
    # TODO: copy entries from above that can also hold dataframes.
]

ND_TYPES = [
    pytest.param(np.asarray, np.ndarray, np.ndarray, id="np.ndarray"),
    # TODO: return 1D sparse array instead of 1D ndarray?
    pytest.param(sp.csr_array, np.ndarray, sp.csr_array, id="sp.csr_array"),
    pytest.param(sp.csc_matrix, np.ndarray, sp.csc_matrix, id="sp.csc_matrix"),
]

DF_TYPES = [
    pytest.param(pd.DataFrame, ExtensionArray, None, id="pd.DataFrame"),
    pytest.param(
        lambda df: XDataset.from_dataframe(df),  # noqa: PLW0108
        XVariable,
        None,
        marks=pytest.mark.skipif(
            not find_spec("xarray"), reason="xarray not installed."
        ),
        id="xr.Dataset",
    ),
]


def convert_ndarrays(adata: AnnData, array_conv: Callable[[Array], Array], /) -> None:
    def conv(v: Array) -> Array:
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
    adata: AnnData, df_conv: Callable[[pd.Series], Array], /
) -> None:
    adata.obs = df_conv(adata.obs)
    adata.var = df_conv(adata.var)


def _expected2np(expected: Array, ad_path: AdRef) -> np.ndarray:
    ndim = len(ad_path.dims)
    match expected:
        case np.ndarray():
            return expected.flatten() if ndim == 1 else expected
        case sp.sparray() | sp.spmatrix():
            return _expected2np(expected.toarray(), ad_path)
        case pd.Series() | XDataArray():
            return expected.to_numpy()
        case _:
            pytest.fail(f"unhandled expected type {type(expected)}")


@pytest.fixture(
    params=[
        pytest.param(
            (ad_path, ad_expected, *typ.values, convert),
            marks=typ.marks,
            id=f"{ad_path}-{typ.id}",
        )
        for paths, types, convert in (
            (ND_PATHS, ND_TYPES, convert_ndarrays),
            (DF_PATHS, DF_TYPES, convert_dataframes),
        )
        for ad_path, ad_expected in paths
        for typ in types
    ]
)
def check_ad_path(
    adata: AnnData, request: pytest.FixtureRequest
) -> tuple[AnnData, AdRef, type[Array], Array]:
    ad_path, ad_expected, convert_array, type_expected_1d, type_expected_2d, convert = (
        request.param
    )
    convert(adata, convert_array)
    return (
        adata,
        ad_path,
        type_expected_1d if len(ad_path.dims) == 1 else type_expected_2d,
        _expected2np(ad_expected(adata), ad_path),
    )


def test_get_values(check_ad_path: tuple[AnnData, AdRef, type[Array], Array]) -> None:
    adata, ad_path, type_expected, expected = check_ad_path

    vals = ad_path(adata)  # TODO: allow returning sparse array?

    assert isinstance(vals, type_expected)
    np.testing.assert_array_equal(asarray(vals), expected, strict=True)
