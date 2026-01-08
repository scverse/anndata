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

    # Test len(), shape, ndim, size, dtype are cheap (metadata only)
    assert len(cats) == n_cats
    assert cats.shape == (n_cats,)
    assert cats.shape[0] == n_cats
    assert cats.ndim == 1
    assert cats.size == n_cats
    assert cats.dtype is not None  # dtype accessible without loading

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

    # Test .codes returns the underlying integer codes
    codes = col.cat.codes
    assert codes is not None
    assert codes.shape[0] == n_cats
    assert list(codes[:5]) == [0, 1, 2, 3, 4]  # integer codes, not category values

    # Test .ordered
    assert col.cat.ordered is False  # default is unordered


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_cat_accessor_ordered(tmp_path: Path, diskfmt: str):
    """Test .cat.ordered for ordered categorical."""
    adata = AnnData(
        X=np.zeros((10, 2)),
        obs=pd.DataFrame({
            "ordered_cat": pd.Categorical(["a", "b", "c"] * 3 + ["a"], ordered=True)
        }),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    assert lazy.obs["ordered_cat"].cat.ordered is True


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
    assert col.cat.codes is None
    assert col.cat.ordered is None


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categories_strided_slicing(tmp_path: Path, diskfmt: str):
    """Test LazyCategories with step != 1 slicing."""
    from anndata.experimental.backed._lazy_arrays import LazyCategories

    n_cats = 20
    categories = [f"Cat_{i:02d}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    cats = lazy.obs["cell_type"].cat.categories
    assert isinstance(cats, LazyCategories)

    # Test positive step > 1
    step2 = cats[::2]
    assert list(step2) == [f"Cat_{i:02d}" for i in range(0, 20, 2)]

    step3 = cats[1::3]
    assert list(step3) == [f"Cat_{i:02d}" for i in range(1, 20, 3)]

    # Test negative step (reverse)
    reversed_cats = cats[::-1]
    assert list(reversed_cats) == [f"Cat_{i:02d}" for i in range(19, -1, -1)]

    # Test negative step with stride
    rev_step2 = cats[::-2]
    assert list(rev_step2) == [f"Cat_{i:02d}" for i in range(19, -1, -2)]

    # Test partial reverse
    partial_rev = cats[10:5:-1]
    assert list(partial_rev) == [f"Cat_{i:02d}" for i in range(10, 5, -1)]


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categories_edge_cases(tmp_path: Path, diskfmt: str):
    """Test LazyCategories edge cases: empty slices, out of bounds, etc."""
    from anndata.experimental.backed._lazy_arrays import LazyCategories

    n_cats = 10
    categories = [f"Cat_{i}" for i in range(n_cats)]
    adata = AnnData(
        X=np.zeros((n_cats, 2)),
        obs=pd.DataFrame({"cell_type": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    cats = lazy.obs["cell_type"].cat.categories
    assert isinstance(cats, LazyCategories)

    # Empty slice
    empty = cats[5:5]
    assert len(empty) == 0

    # Out of bounds slice (should clamp like numpy)
    oob = cats[100:200]
    assert len(oob) == 0

    oob2 = cats[5:100]
    assert list(oob2) == [f"Cat_{i}" for i in range(5, 10)]

    # Single element slice
    single = cats[3:4]
    assert len(single) == 1
    assert single[0] == "Cat_3"

    # Out of bounds single index should raise
    with pytest.raises(IndexError):
        _ = cats[100]

    with pytest.raises(IndexError):
        _ = cats[-100]


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categories_methods(tmp_path: Path, diskfmt: str):
    """Test LazyCategories helper methods: __contains__, tolist, __array__."""
    from anndata.experimental.backed._lazy_arrays import LazyCategories

    categories = ["apple", "banana", "cherry"]
    adata = AnnData(
        X=np.zeros((3, 2)),
        obs=pd.DataFrame({"fruit": pd.Categorical(categories)}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    cats = lazy.obs["fruit"].cat.categories
    assert isinstance(cats, LazyCategories)

    # Test __contains__
    assert "apple" in cats
    assert "banana" in cats
    assert "grape" not in cats

    # Test tolist
    cat_list = cats.tolist()
    assert isinstance(cat_list, list)
    assert cat_list == categories

    # Test __array__ with copy parameter
    arr1 = np.array(cats)
    arr2 = np.array(cats, copy=True)
    assert list(arr1) == categories
    assert list(arr2) == categories

    # Test __array__ with dtype
    arr_obj = np.array(cats, dtype=object)
    assert arr_obj.dtype == object


@pytest.mark.parametrize("diskfmt", ["zarr", "h5ad"])
def test_lazy_categories_nullable_strings(tmp_path: Path, diskfmt: str):
    """Test LazyCategories with nullable string categories (ArrowStringArray)."""
    pytest.importorskip("pyarrow")

    # Create categorical with nullable string dtype
    categories = pd.array(["type_a", "type_b", "type_c"], dtype="string[pyarrow]")
    cat_data = pd.Categorical.from_codes([0, 1, 2, 0, 1], categories=categories)

    adata = AnnData(
        X=np.zeros((5, 2)),
        obs=pd.DataFrame({"cell_type": cat_data}),
    )

    path = tmp_path / f"test.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)

    lazy = read_lazy(path)
    cats = lazy.obs["cell_type"].cat.categories

    # Should still work - may return ArrowStringArray or converted array
    assert len(cats) == 3
    # Check we can access elements
    assert cats[0] in ["type_a", "type_b", "type_c"]
    assert cats[-1] in ["type_a", "type_b", "type_c"]
