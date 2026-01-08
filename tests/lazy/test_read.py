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


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categorical_dtype_n_categories(tmp_path: Path, diskfmt: str):
    """Test LazyCategoricalDtype.n_categories is cheap (metadata only)."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    n_cats = 100
    categories = [f"Cat_{i:03d}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cell_type"].dtype

    # dtype should be LazyCategoricalDtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # n_categories should work without loading all categories
    assert dtype.n_categories == n_cats

    # ordered should be accessible
    assert dtype.ordered is False


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categorical_dtype_head_tail_categories(tmp_path: Path, diskfmt: str):
    """Test LazyCategoricalDtype.head_categories and tail_categories for partial reads."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    n_cats = 50
    categories = [f"Type_{i:02d}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cell_type"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Test head_categories (first n)
    first5 = dtype.head_categories(5)
    assert len(first5) == 5
    assert list(first5) == [f"Type_{i:02d}" for i in range(5)]

    # Test head_categories default (first 5)
    default_head = dtype.head_categories()
    assert len(default_head) == 5
    assert list(default_head) == [f"Type_{i:02d}" for i in range(5)]

    # Test tail_categories (last n)
    last3 = dtype.tail_categories(3)
    assert len(last3) == 3
    assert list(last3) == [f"Type_{i:02d}" for i in range(47, 50)]

    # Test tail_categories default (last 5)
    default_tail = dtype.tail_categories()
    assert len(default_tail) == 5
    assert list(default_tail) == [f"Type_{i:02d}" for i in range(45, 50)]

    # Test requesting more than available
    all_head = dtype.head_categories(100)
    assert len(all_head) == n_cats
    assert list(all_head) == categories

    all_tail = dtype.tail_categories(100)
    assert len(all_tail) == n_cats
    assert list(all_tail) == categories


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categorical_dtype_categories_caching(tmp_path: Path, diskfmt: str):
    """Test that categories are cached after full load."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["a", "b", "c", "d", "e"]
    adata = AnnData(
        X=np.zeros((5, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Before loading, categories should not be cached
    # (accessing internal state for testing)
    assert dtype._LazyCategoricalDtype__categories is None

    # Load categories
    cats = dtype.categories
    assert cats is not None
    assert list(cats) == categories

    # After loading, should be cached
    assert dtype._LazyCategoricalDtype__categories is not None

    # Verify head/tail_categories use cache by modifying cache
    dtype._LazyCategoricalDtype__categories = pd.Index(["x", "y", "z", "w", "v"])
    head = dtype.head_categories(3)
    assert list(head) == ["x", "y", "z"]  # Returns cached values, not disk values
    tail = dtype.tail_categories(3)
    assert list(tail) == ["z", "w", "v"]  # Returns cached values, not disk values


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categorical_dtype_ordered(tmp_path: Path, diskfmt: str):
    """Test LazyCategoricalDtype with ordered categories."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    adata = AnnData(
        X=np.zeros((10, 2)),
        obs=pd.DataFrame({
            "ordered_cat": pd.Categorical(
                ["low", "medium", "high"] * 3 + ["low"],
                categories=["low", "medium", "high"],
                ordered=True,
            )
        }),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["ordered_cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    assert dtype.ordered is True
    assert dtype.n_categories == 3
    assert list(dtype.categories) == ["low", "medium", "high"]


def test_lazy_categorical_dtype_repr(tmp_path: Path):
    """Test LazyCategoricalDtype repr before and after loading."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = [f"cat_{i}" for i in range(100)]
    adata = AnnData(
        X=np.zeros((100, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Before loading: lazy repr
    repr_before = repr(dtype)
    assert "LazyCategoricalDtype" in repr_before
    assert "n_categories=100" in repr_before

    # Load categories
    _ = dtype.categories

    # After loading: standard CategoricalDtype repr
    repr_after = repr(dtype)
    assert "CategoricalDtype" in repr_after


def test_lazy_categorical_dtype_equality(tmp_path: Path):
    """Test LazyCategoricalDtype equality comparisons."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["a", "b", "c"]
    adata = AnnData(
        X=np.zeros((3, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

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

    # Test comparison with non-CategoricalDtype
    assert dtype != np.dtype("int64")
    assert dtype != 123
    assert dtype is not None


def test_lazy_categorical_dtype_equality_same_array(tmp_path: Path):
    """Test LazyCategoricalDtype equality between instances with same underlying array."""

    categories = ["x", "y", "z"]
    adata = AnnData(
        X=np.zeros((3, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype1 = lazy.obs["cat"].dtype
    dtype2 = lazy.obs["cat"].dtype  # Same underlying array

    # Same object should be equal
    assert dtype1 is dtype2  # They are the same instance
    assert dtype1 == dtype2


def test_lazy_categorical_dtype_hash(tmp_path: Path):
    """Test LazyCategoricalDtype is hashable."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["a", "b", "c"]
    adata = AnnData(
        X=np.zeros((3, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Should be hashable (required for pandas internals)
    h = hash(dtype)
    assert isinstance(h, int)

    # Can be used in a set
    s = {dtype}
    assert dtype in s


def test_lazy_categorical_dtype_n_categories_from_cache(tmp_path: Path):
    """Test n_categories returns from cache when categories already loaded."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["a", "b", "c", "d", "e"]
    adata = AnnData(
        X=np.zeros((5, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # Load categories first
    cats = dtype.categories
    assert cats is not None

    # Verify n_categories uses cache by modifying cache
    dtype._LazyCategoricalDtype__categories = pd.Index(["x", "y", "z"])
    assert dtype.n_categories == 3  # Returns cached length, not disk length


def test_lazy_categorical_dtype_empty_array():
    """Test LazyCategoricalDtype with None categories_array."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    # Create dtype with None categories_array
    dtype = LazyCategoricalDtype(categories_array=None, ordered=False)

    # Properties should handle None gracefully
    assert dtype.n_categories == 0
    assert dtype.categories is None

    # head_categories and tail_categories should return empty arrays
    head = dtype.head_categories(5)
    assert len(head) == 0

    tail = dtype.tail_categories(5)
    assert len(tail) == 0

    # repr should still work
    r = repr(dtype)
    assert "LazyCategoricalDtype" in r
    assert "n_categories=0" in r


def test_lazy_categorical_dtype_name(tmp_path: Path):
    """Test LazyCategoricalDtype.name property."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    categories = ["a", "b"]
    adata = AnnData(
        X=np.zeros((2, 2)),
        obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
    )

    path = tmp_path / "test.zarr"
    adata.write_zarr(path)

    lazy = read_lazy(path)
    dtype = lazy.obs["cat"].dtype
    assert isinstance(dtype, LazyCategoricalDtype)

    # name should be "category"
    assert dtype.name == "category"


def test_lazy_categorical_dtype_equality_with_none_categories(tmp_path: Path):
    """Test LazyCategoricalDtype equality when comparing dtypes with None categories."""
    from anndata.experimental.backed._lazy_arrays import LazyCategoricalDtype

    # Create dtype with None categories
    dtype1 = LazyCategoricalDtype(categories_array=None, ordered=False)

    # Regular CategoricalDtype without categories set
    dtype2 = pd.CategoricalDtype(categories=None, ordered=False)

    # Both have None categories, should be equal
    assert dtype1 == dtype2


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
