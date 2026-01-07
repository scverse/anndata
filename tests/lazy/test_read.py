from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import zarr

from anndata import AnnData
from anndata.compat import DaskArray
from anndata.experimental import read_elem_lazy, read_lazy
from anndata.experimental.backed._io import ANNDATA_ELEMS
from anndata.io import read_zarr, write_elem
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    AccessTrackingStore,
    assert_equal,
    gen_adata,
    gen_typed_df,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from anndata._types import AnnDataElem


pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.mark.parametrize(
    ("elem_key", "sub_key"),
    [
        ("raw", "X"),
        ("obs", "cat"),
        ("obs", "int64"),
        *((elem_name, None) for elem_name in ANNDATA_ELEMS),
    ],
)
def test_access_count_elem_access(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
    elem_key: AnnDataElem,
    sub_key: str,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    full_path = f"{elem_key}/{sub_key}" if sub_key is not None else elem_key
    remote_store_tall_skinny.initialize_key_trackers({full_path, "X"})
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(adata_remote_tall_skinny), elem_key)
    if sub_key is not None:
        if elem_key in {"obs", "var"}:
            elem[sub_key]
        else:
            getattr(elem, sub_key)
    remote_store_tall_skinny.assert_access_count(full_path, 0)
    remote_store_tall_skinny.assert_access_count("X", 0)


def test_access_count_subset(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
):
    non_obs_elem_names = filter(lambda e: e != "obs", ANNDATA_ELEMS)
    remote_store_tall_skinny.initialize_key_trackers([
        "obs/cat/codes",
        *non_obs_elem_names,
    ])
    adata_remote_tall_skinny[adata_remote_tall_skinny.obs["cat"] == "a", :]
    # all codes read in for subset (from 4 chunks as set in the fixture)
    remote_store_tall_skinny.assert_access_count("obs/cat/codes", 4)
    for elem_name in non_obs_elem_names:
        remote_store_tall_skinny.assert_access_count(elem_name, 0)


def test_access_count_subset_column_compute(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
):
    remote_store_tall_skinny.initialize_key_trackers(["obs/int64"])
    adata_remote_tall_skinny[adata_remote_tall_skinny.shape[0] // 2, :].obs[
        "int64"
    ].compute()
    # two chunks needed for 0:10 subset
    remote_store_tall_skinny.assert_access_count("obs/int64", 1)


def test_access_count_index(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_with_store_tall_skinny_path: Path,
) -> None:
    adata_orig = read_zarr(adata_remote_with_store_tall_skinny_path)

    remote_store_tall_skinny.initialize_key_trackers(["obs/_index"])
    read_lazy(remote_store_tall_skinny, load_annotation_index=False)
    remote_store_tall_skinny.assert_access_count("obs/_index", 0)

    read_lazy(remote_store_tall_skinny)
    n_chunks = 4
    count_expected = (  # *2 when mask exists
        n_chunks * 2 if adata_orig.obs.index.dtype == "string" else n_chunks
    )
    remote_store_tall_skinny.assert_access_count("obs/_index", count_expected)


def test_access_count_dtype(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
    adata_remote_with_store_tall_skinny_path: Path,
) -> None:
    adata_orig = read_zarr(adata_remote_with_store_tall_skinny_path)

    remote_store_tall_skinny.initialize_key_trackers(["obs/cat/categories"])
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 0)

    count_expected = 2 if adata_orig.obs["cat"].cat.categories.dtype == "string" else 1
    # This should only cause categories to be read in once (and their mask if applicable)
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", count_expected)
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", count_expected)


def test_uns_uses_dask(adata_remote: AnnData):
    assert isinstance(adata_remote.uns["nested"]["nested_further"]["array"], DaskArray)


def test_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    remote_to_memory = adata_remote.to_memory()
    assert_equal(remote_to_memory, adata_orig)


def test_access_counts_obsm_df(tmp_path: Path):
    adata = AnnData(
        X=np.array(np.random.rand(100, 20)),
    )
    adata.obsm["df"] = pd.DataFrame(
        {"col1": np.random.rand(100), "col2": np.random.rand(100)},
        index=adata.obs_names,
    )
    adata.write_zarr(tmp_path)
    store = AccessTrackingStore(tmp_path)
    store.initialize_key_trackers(["obsm/df"])
    read_lazy(store, load_annotation_index=False)
    store.assert_access_count("obsm/df", 0)


def test_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    obs_cats = adata_orig.obs["obs_cat"].cat.categories
    subset_obs = adata_orig.obs["obs_cat"] == obs_cats[0]
    assert_equal(adata_orig[subset_obs, :], adata_remote[subset_obs, :].to_memory())

    var_cats = adata_orig.var["var_cat"].cat.categories
    subset_var = adata_orig.var["var_cat"] == var_cats[0]
    assert_equal(adata_orig[:, subset_var], adata_remote[:, subset_var].to_memory())


def test_view_of_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    cats_obs = adata_orig.obs["obs_cat"].cat.categories
    subset_obs = (adata_orig.obs["obs_cat"] == cats_obs[0]) | (
        adata_orig.obs["obs_cat"] == cats_obs[1]
    )
    subsetted_adata = adata_orig[subset_obs, :]
    subset_subset_obs = subsetted_adata.obs["obs_cat"] == cats_obs[1]
    subsetted_subsetted_adata = subsetted_adata[subset_subset_obs, :]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[subset_obs, :][subset_subset_obs, :].to_memory(),
    )

    cats_var = adata_orig.var["var_cat"].cat.categories
    subset_var = (adata_orig.var["var_cat"] == cats_var[0]) | (
        adata_orig.var["var_cat"] == cats_var[1]
    )
    subsetted_adata = adata_orig[:, subset_var]
    subset_subset_var = subsetted_adata.var["var_cat"] == cats_var[1]
    subsetted_subsetted_adata = subsetted_adata[:, subset_subset_var]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[:, subset_var][:, subset_subset_var].to_memory(),
    )


@pytest.mark.zarr_io
def test_unconsolidated(tmp_path: Path, mtx_format):
    adata = gen_adata((10, 10), mtx_format, **GEN_ADATA_NO_XARRAY_ARGS)
    orig_pth = tmp_path / "orig.zarr"
    adata.write_zarr(orig_pth)
    (orig_pth / ".zmetadata").unlink()
    store = AccessTrackingStore(orig_pth)
    store.initialize_key_trackers(["obs/.zgroup", ".zgroup"])
    with pytest.warns(UserWarning, match=r"Did not read zarr as consolidated"):
        remote = read_lazy(store)
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, adata)
    store.assert_access_count("obs/.zgroup", 1)


def test_h5_file_obj(tmp_path: Path):
    adata = gen_adata((10, 10), **GEN_ADATA_NO_XARRAY_ARGS)
    orig_pth = tmp_path / "adata.h5ad"
    adata.write_h5ad(orig_pth)
    remote = read_lazy(orig_pth)
    assert remote.file.is_open
    assert remote.filename == orig_pth
    assert_equal(remote.to_memory(), adata)


@pytest.fixture(scope="session")
def df_group(tmp_path_factory) -> zarr.Group:
    df = gen_typed_df(120)
    path = tmp_path_factory.mktemp("foo.zarr")
    g = zarr.open_group(path, mode="w", zarr_format=2)
    write_elem(g, "foo", df, dataset_kwargs={"chunks": 25})
    return zarr.open(path, mode="r")["foo"]


@pytest.mark.parametrize(
    ("chunks", "expected_chunks"),
    [((1,), (1,)), ((-1,), (120,)), (None, (25,))],
    ids=["small", "minus_one_uses_full", "none_uses_ondisk_chunking"],
)
def test_chunks_df(
    tmp_path: Path,
    chunks: tuple[int] | None,
    expected_chunks: tuple[int],
    df_group: zarr.Group,
):
    ds = read_elem_lazy(df_group, chunks=chunks)
    for k in ds:
        if isinstance(arr := ds[k].data, DaskArray):
            assert arr.chunksize == expected_chunks


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categories(tmp_path: Path, diskfmt: str):
    """Test LazyCategories slicing via .cat.categories accessor."""
    from anndata.experimental.backed._lazy_arrays import (
        CategoricalArray,
        LazyCategories,
    )

    n_cats = 100
    # Use zero-padded numbers to ensure consistent ordering
    categories = [f"Type_{i:03d}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    col = lazy.obs["cell_type"]

    # Verify underlying CategoricalArray exists
    cat_arr = col.variable._data.array
    assert isinstance(cat_arr, CategoricalArray)

    # Test .cat.categories returns LazyCategories
    cats = col.cat.categories
    assert isinstance(cats, LazyCategories)

    # Test len() is cheap (n_categories)
    assert len(cats) == n_cats

    # Test slicing - head
    head5 = cats[:5]
    assert len(head5) == 5
    assert list(head5) == [f"Type_{i:03d}" for i in range(5)]

    # Test slicing - tail
    tail5 = cats[-5:]
    assert len(tail5) == 5
    assert list(tail5) == [f"Type_{i:03d}" for i in range(95, 100)]

    # Test slicing - middle
    middle = cats[10:20]
    assert len(middle) == 10
    assert list(middle) == [f"Type_{i:03d}" for i in range(10, 20)]

    # Test single item access
    assert cats[0] == "Type_000"
    assert cats[-1] == "Type_099"

    # Test full array conversion
    all_cats = np.array(cats)
    assert len(all_cats) == n_cats


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_cat_accessor(tmp_path: Path, diskfmt: str):
    """Test the .cat accessor for user-facing API."""
    from anndata.experimental.backed._lazy_arrays import LazyCategories

    n_cats = 50
    # Use zero-padded numbers to ensure consistent ordering
    categories = [f"Cell_{i:02d}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    col = lazy.obs["cell_type"]

    # Test .cat accessor exists
    assert hasattr(col, "cat")

    # Test .categories returns LazyCategories
    cats = col.cat.categories
    assert isinstance(cats, LazyCategories)

    # Test len() via accessor
    assert len(cats) == n_cats

    # Test slicing via accessor - head
    head5 = cats[:5]
    assert len(head5) == 5
    assert list(head5) == [f"Cell_{i:02d}" for i in range(5)]

    # Test slicing via accessor - tail
    tail3 = cats[-3:]
    assert len(tail3) == 3
    assert list(tail3) == [f"Cell_{i:02d}" for i in range(47, 50)]

    # Test slicing via accessor - middle
    middle = cats[10:15]
    assert len(middle) == 5
    assert list(middle) == [f"Cell_{i:02d}" for i in range(10, 15)]

    # Test full array conversion
    all_cats = np.array(cats)
    assert len(all_cats) == n_cats


def test_cat_accessor_non_categorical(tmp_path: Path):
    """Test that .cat accessor returns None for non-categorical columns."""
    adata = AnnData(
        X=np.zeros((10, 2)),
        obs=pd.DataFrame({"numeric": np.arange(10)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    col = lazy.obs["numeric"]

    # Accessor returns None for non-categorical columns
    assert col.cat.categories is None
