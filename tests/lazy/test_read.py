from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING
from unittest.mock import patch

import h5py
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
    remote_store_tall_skinny.initialize_key_trackers(["obs/cat/categories"])
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 0)

    # Accessing dtype alone should NOT load categories (lazy loading)
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 0)

    # n_categories should also be cheap (metadata only)
    _ = adata_remote_tall_skinny.obs["cat"].dtype.n_categories
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 0)

    # Accessing categories should trigger loading (once, then cached)
    count_before = remote_store_tall_skinny.get_access_count("obs/cat/categories")
    _ = adata_remote_tall_skinny.obs["cat"].dtype.categories
    count_after = remote_store_tall_skinny.get_access_count("obs/cat/categories")
    assert count_after > count_before, "categories access should trigger read"

    # Subsequent accesses should use cache (no additional reads)
    _ = adata_remote_tall_skinny.obs["cat"].dtype.categories
    _ = adata_remote_tall_skinny.obs["cat"].dtype.categories
    assert (
        remote_store_tall_skinny.get_access_count("obs/cat/categories") == count_after
    ), "cached categories should not trigger additional reads"


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


# Session-scoped fixtures for categorical data (write once, read many)
# Each category type has zarr and h5ad path fixtures, plus a parametrized store fixture


def _write_categorical_zarr(tmp_path_factory, name: str, cat: pd.Categorical) -> Path:
    """Helper to write categorical to zarr and return path."""
    path = tmp_path_factory.mktemp(f"{name}.zarr")
    store = zarr.open(path, mode="w")
    write_elem(store, "cat", cat)
    return path


def _write_categorical_h5ad(tmp_path_factory, name: str, cat: pd.Categorical) -> Path:
    """Helper to write categorical to h5ad and return path."""
    path = tmp_path_factory.mktemp(name) / "cat.h5ad"
    with h5py.File(path, mode="w") as f:
        write_elem(f, "cat", cat)
    return path


def _open_categorical_store(path: Path, backend: str):
    """Helper to open categorical store for either backend."""
    if backend == "zarr":
        return zarr.open(path, mode="r")["cat"]
    else:
        return h5py.File(path, mode="r")["cat"]


# Small categorical ['a', 'b', 'c']
@pytest.fixture(scope="session")
def cat_small_path_zarr(tmp_path_factory) -> Path:
    return _write_categorical_zarr(
        tmp_path_factory, "cat_small", pd.Categorical(["a", "b", "c"])
    )


@pytest.fixture(scope="session")
def cat_small_path_h5ad(tmp_path_factory) -> Path:
    return _write_categorical_h5ad(
        tmp_path_factory, "cat_small", pd.Categorical(["a", "b", "c"])
    )


@pytest.fixture(params=["zarr", "h5ad"])
def cat_small_store(request, cat_small_path_zarr: Path, cat_small_path_h5ad: Path):
    """Parametrized fixture: small categorical ['a', 'b', 'c'] for both backends."""
    path = cat_small_path_zarr if request.param == "zarr" else cat_small_path_h5ad
    store = _open_categorical_store(path, request.param)
    yield store
    if request.param == "h5ad":
        store.file.close()


# Medium categorical ['a', 'b', 'c', 'd', 'e']
@pytest.fixture(scope="session")
def cat_medium_path_zarr(tmp_path_factory) -> Path:
    return _write_categorical_zarr(
        tmp_path_factory, "cat_medium", pd.Categorical(["a", "b", "c", "d", "e"])
    )


@pytest.fixture(scope="session")
def cat_medium_path_h5ad(tmp_path_factory) -> Path:
    return _write_categorical_h5ad(
        tmp_path_factory, "cat_medium", pd.Categorical(["a", "b", "c", "d", "e"])
    )


@pytest.fixture(params=["zarr", "h5ad"])
def cat_medium_store(request, cat_medium_path_zarr: Path, cat_medium_path_h5ad: Path):
    """Parametrized fixture: medium categorical for both backends."""
    path = cat_medium_path_zarr if request.param == "zarr" else cat_medium_path_h5ad
    store = _open_categorical_store(path, request.param)
    yield store
    if request.param == "h5ad":
        store.file.close()


# Large categorical with 100 categories
@pytest.fixture(scope="session")
def cat_large_path_zarr(tmp_path_factory) -> Path:
    categories = [f"cat_{i}" for i in range(100)]
    return _write_categorical_zarr(
        tmp_path_factory, "cat_large", pd.Categorical(categories)
    )


@pytest.fixture(scope="session")
def cat_large_path_h5ad(tmp_path_factory) -> Path:
    categories = [f"cat_{i}" for i in range(100)]
    return _write_categorical_h5ad(
        tmp_path_factory, "cat_large", pd.Categorical(categories)
    )


@pytest.fixture(params=["zarr", "h5ad"])
def cat_large_store(request, cat_large_path_zarr: Path, cat_large_path_h5ad: Path):
    """Parametrized fixture: large categorical (100 categories) for both backends."""
    path = cat_large_path_zarr if request.param == "zarr" else cat_large_path_h5ad
    store = _open_categorical_store(path, request.param)
    yield store
    if request.param == "h5ad":
        store.file.close()


# Ordered categorical ['low', 'medium', 'high']
@pytest.fixture(scope="session")
def cat_ordered_path_zarr(tmp_path_factory) -> Path:
    cat = pd.Categorical(
        ["low", "medium", "high"] * 3 + ["low"],
        categories=["low", "medium", "high"],
        ordered=True,
    )
    return _write_categorical_zarr(tmp_path_factory, "cat_ordered", cat)


@pytest.fixture(scope="session")
def cat_ordered_path_h5ad(tmp_path_factory) -> Path:
    cat = pd.Categorical(
        ["low", "medium", "high"] * 3 + ["low"],
        categories=["low", "medium", "high"],
        ordered=True,
    )
    return _write_categorical_h5ad(tmp_path_factory, "cat_ordered", cat)


@pytest.fixture(params=["zarr", "h5ad"])
def cat_ordered_store(
    request, cat_ordered_path_zarr: Path, cat_ordered_path_h5ad: Path
):
    """Parametrized fixture: ordered categorical for both backends."""
    path = cat_ordered_path_zarr if request.param == "zarr" else cat_ordered_path_h5ad
    store = _open_categorical_store(path, request.param)
    yield store
    if request.param == "h5ad":
        store.file.close()


# 50 categories for head/tail testing
@pytest.fixture(scope="session")
def cat_fifty_path_zarr(tmp_path_factory) -> Path:
    categories = [f"Type_{i:02d}" for i in range(50)]
    return _write_categorical_zarr(
        tmp_path_factory, "cat_fifty", pd.Categorical(categories)
    )


@pytest.fixture(scope="session")
def cat_fifty_path_h5ad(tmp_path_factory) -> Path:
    categories = [f"Type_{i:02d}" for i in range(50)]
    return _write_categorical_h5ad(
        tmp_path_factory, "cat_fifty", pd.Categorical(categories)
    )


@pytest.fixture(params=["zarr", "h5ad"])
def cat_fifty_store(request, cat_fifty_path_zarr: Path, cat_fifty_path_h5ad: Path):
    """Parametrized fixture: 50 categories for head/tail testing, both backends."""
    path = cat_fifty_path_zarr if request.param == "zarr" else cat_fifty_path_h5ad
    store = _open_categorical_store(path, request.param)
    yield store
    if request.param == "h5ad":
        store.file.close()


def test_lazy_categorical_dtype_n_categories(cat_large_store):
    """Test n_categories is cheap (metadata only) and uses cache when loaded."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_large_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Before loading: n_categories should work without loading categories
    assert "categories" not in dtype.__dict__
    assert dtype.n_categories == 100
    assert "categories" not in dtype.__dict__  # Still not loaded - proves metadata-only
    assert dtype.ordered is False

    # After loading: n_categories should use cache
    _ = dtype.categories  # Force load
    assert "categories" in dtype.__dict__
    assert dtype.n_categories == 100  # Uses cache now

    # Verify cache is used by modifying cached value
    dtype.__dict__["categories"] = pd.Index(["x", "y", "z"])
    assert dtype.n_categories == 3  # Returns cached length, not disk length


def test_lazy_categorical_dtype_head_tail_categories(cat_fifty_store):
    """Test head_categories and tail_categories perform partial reads without loading all."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_fifty_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Verify categories not loaded initially
    assert "categories" not in dtype.__dict__

    # Test head_categories (first n) - should NOT load all categories
    first5 = dtype.head_categories(5)
    assert len(first5) == 5
    assert list(first5) == [f"Type_{i:02d}" for i in range(5)]
    assert "categories" not in dtype.__dict__  # Still not fully loaded

    # Test head_categories default (first 5)
    default_head = dtype.head_categories()
    assert len(default_head) == 5
    assert list(default_head) == [f"Type_{i:02d}" for i in range(5)]
    assert "categories" not in dtype.__dict__  # Still not fully loaded

    # Test tail_categories (last n) - should NOT load all categories
    last3 = dtype.tail_categories(3)
    assert len(last3) == 3
    assert list(last3) == [f"Type_{i:02d}" for i in range(47, 50)]
    assert "categories" not in dtype.__dict__  # Still not fully loaded

    # Test tail_categories default (last 5)
    default_tail = dtype.tail_categories()
    assert len(default_tail) == 5
    assert list(default_tail) == [f"Type_{i:02d}" for i in range(45, 50)]
    assert "categories" not in dtype.__dict__  # Still not fully loaded

    # Test requesting more than available
    all_head = dtype.head_categories(100)
    assert len(all_head) == 50
    assert list(all_head) == [f"Type_{i:02d}" for i in range(50)]

    all_tail = dtype.tail_categories(100)
    assert len(all_tail) == 50
    assert list(all_tail) == [f"Type_{i:02d}" for i in range(50)]


def test_lazy_categorical_dtype_categories_caching(cat_medium_store):
    """Test that categories are cached after full load."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_medium_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Before loading, categories should not be cached (uses @cached_property)
    assert "categories" not in dtype.__dict__

    # Load categories
    cats = dtype.categories
    assert cats is not None
    assert list(cats) == ["a", "b", "c", "d", "e"]

    # After loading, should be cached in __dict__ (cached_property pattern)
    assert "categories" in dtype.__dict__

    # Verify head/tail_categories use cache by modifying the cached value
    dtype.__dict__["categories"] = pd.Index(["x", "y", "z", "w", "v"])
    head = dtype.head_categories(3)
    assert list(head) == ["x", "y", "z"]  # Returns cached values, not disk values
    tail = dtype.tail_categories(3)
    assert list(tail) == ["z", "w", "v"]  # Returns cached values, not disk values


def test_lazy_categorical_dtype_ordered(cat_ordered_store):
    """Test LazyCategoricalDtype with ordered categories."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_ordered_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    assert dtype.ordered is True
    assert dtype.n_categories == 3
    assert list(dtype.categories) == ["low", "medium", "high"]


def test_lazy_categorical_dtype_repr(cat_large_store, cat_small_store):
    """Test LazyCategoricalDtype repr shows truncated categories."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    # Test large number of categories (truncated repr)
    lazy_cat = read_elem_lazy(cat_large_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    repr_str = repr(dtype)
    assert "LazyCategoricalDtype" in repr_str
    assert "n=100" in repr_str
    assert "..." in repr_str  # Truncation indicator
    assert "cat_0" in repr_str  # Head category
    assert "cat_99" in repr_str  # Tail category

    # Test small number of categories (full repr)
    small_lazy_cat = read_elem_lazy(cat_small_store)
    small_dtype = small_lazy_cat.dtype

    small_repr = repr(small_dtype)
    assert "LazyCategoricalDtype" in small_repr
    assert "..." not in small_repr  # No truncation for small categories
    assert "'a'" in small_repr
    assert "'b'" in small_repr
    assert "'c'" in small_repr


def test_lazy_categorical_dtype_equality(cat_small_store):
    """Test LazyCategoricalDtype equality comparisons and basic properties."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_small_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Test name property (inherited from CategoricalDtype)
    assert dtype.name == "category"

    # Test string comparison (dtype == "category")
    assert dtype == "category"
    assert dtype != "int64"

    # Test comparison with regular CategoricalDtype
    regular_dtype = pd.CategoricalDtype(categories=["a", "b", "c"], ordered=False)
    assert dtype == regular_dtype

    # Test comparison with different categories
    different_dtype = pd.CategoricalDtype(categories=["x", "y", "z"], ordered=False)
    assert dtype != different_dtype

    # Test comparison with different ordered flag
    ordered_dtype = pd.CategoricalDtype(categories=["a", "b", "c"], ordered=True)
    assert dtype != ordered_dtype

    # Test comparison with CategoricalDtype with None categories
    # LazyCategoricalDtype always has categories, so should not equal None-categories dtype
    dtype_none = pd.CategoricalDtype(categories=None, ordered=False)
    assert dtype != dtype_none

    # Test comparison with non-CategoricalDtype
    assert dtype != np.dtype("int64")
    assert dtype != 123
    assert dtype is not None


@pytest.mark.parametrize("backend", ["zarr", "h5ad"])
def test_lazy_categorical_dtype_equality_no_load(
    cat_small_path_zarr: Path, cat_small_path_h5ad: Path, backend: str
):
    """Test same-location equality doesn't load category data.

    Both h5py (HDF5 object ID comparison) and zarr 3.x (StorePath comparison) use
    location-based equality that doesn't read array contents. This test verifies
    that behavior by patching __getitem__ to raise if called.
    """
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    if backend == "zarr":
        path = cat_small_path_zarr

        def open_store(p):
            return zarr.open(p, mode="r")["cat"]

    else:
        path = cat_small_path_h5ad
        # Keep h5py files open for the duration of the test
        open_store = lambda p: h5py.File(p, mode="r")["cat"]

    # Open the same file twice to get different Python objects pointing to same location
    store1 = open_store(path)
    store2 = open_store(path)
    dtype1 = read_elem_lazy(store1).dtype
    dtype2 = read_elem_lazy(store2).dtype

    assert isinstance(dtype1, LazyCategoricalDtype)
    assert isinstance(dtype2, LazyCategoricalDtype)
    # Verify these are different Python objects
    assert dtype1._categories_elem is not dtype2._categories_elem

    # Patch __getitem__ to raise if data is loaded during comparison
    cat_arr1 = dtype1._get_categories_array()
    cat_arr2 = dtype2._get_categories_array()

    with (
        patch.object(
            cat_arr1,
            "__getitem__",
            side_effect=AssertionError("Data was loaded from arr1"),
        ),
        patch.object(
            cat_arr2,
            "__getitem__",
            side_effect=AssertionError("Data was loaded from arr2"),
        ),
    ):
        # This should use location-based comparison without triggering __getitem__
        assert dtype1 == dtype2

    # Also verify our cache wasn't populated
    assert "categories" not in dtype1.__dict__
    assert "categories" not in dtype2.__dict__

    # Clean up h5py file handles
    if backend == "h5ad":
        store1.file.close()
        store2.file.close()


def test_lazy_categorical_roundtrip_via_anndata(tmp_path: Path):
    """Integration test: lazy categorical through full AnnData workflow.

    This test uses the full AnnData read/write path rather than write_elem/read_elem_lazy
    to verify end-to-end integration including dtype caching and equality.
    """
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["type_a", "type_b", "type_c"]
    adata = AnnData(
        X=np.zeros((6, 2)),
        obs=pd.DataFrame({
            "cat": pd.Categorical(categories * 2),
            "ordered_cat": pd.Categorical(
                ["low", "high"] * 3,
                categories=["low", "high"],
                ordered=True,
            ),
        }),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    # Read lazy and verify dtype
    lazy = read_lazy(path)
    dtype1 = lazy.obs["cat"].dtype
    dtype2 = lazy.obs["cat"].dtype  # Same underlying array

    assert isinstance(dtype1, LazyCategoricalDtype)
    assert dtype1 is dtype2  # Same instance (cached)
    assert dtype1 == dtype2

    # Verify ordered categorical
    ordered_dtype = lazy.obs["ordered_cat"].dtype
    assert isinstance(ordered_dtype, LazyCategoricalDtype)
    assert ordered_dtype.ordered is True

    # Round-trip: lazy -> memory should equal original
    loaded = lazy.to_memory()
    assert loaded.obs["cat"].equals(adata.obs["cat"])
    assert loaded.obs["ordered_cat"].equals(adata.obs["ordered_cat"])


def test_lazy_categorical_dtype_hash(cat_small_store):
    """Test LazyCategoricalDtype is hashable."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    lazy_cat = read_elem_lazy(cat_small_store)
    dtype = lazy_cat.dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Should be hashable (useful for collecting unique dtypes in sets/dicts)
    h = hash(dtype)
    assert isinstance(h, int)

    # Can be used in a set
    s = {dtype}
    assert dtype in s


def test_nullable_string_index_decoding(tmp_path: Path):
    """Test that nullable string indices are properly decoded from bytes.

    HDF5 stores strings as bytes. When reading nullable-string-array indices,
    they should be decoded to proper strings, not converted using str() which
    would produce "b'...'" representations.

    Regression test for https://github.com/scverse/anndata/issues/2271
    """
    expected_obs = ["cell_A", "cell_B", "cell_C", "cell_D", "cell_E"]
    expected_var = ["gene_X", "gene_Y", "gene_Z"]

    adata = AnnData(
        np.zeros((len(expected_obs), len(expected_var))),
        obs=dict(obs_names=expected_obs),
        var=dict(var_names=expected_var),
    )

    path = tmp_path / "test.h5ad"
    adata.write_h5ad(path)

    lazy = read_lazy(path)
    obs_names = list(lazy.obs_names)
    var_names = list(lazy.var_names)

    assert obs_names == expected_obs
    assert var_names == expected_var
