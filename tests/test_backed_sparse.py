from __future__ import annotations

from functools import partial
from itertools import product
from typing import TYPE_CHECKING, Literal, get_args

import h5py
import numpy as np
import pytest
import zarr
from scipy import sparse

import anndata as ad
from anndata._core.anndata import AnnData
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs.registry import read_elem_lazy
from anndata._io.zarr import open_write_group
from anndata.compat import (
    CSArray,
    CSMatrix,
    DaskArray,
    ZarrGroup,
    is_zarr_v2,
)
from anndata.experimental import read_dispatched
from anndata.tests.helpers import AccessTrackingStore, assert_equal, subset_func

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Sequence
    from pathlib import Path
    from types import EllipsisType

    from _pytest.mark import ParameterSet
    from numpy.typing import ArrayLike, NDArray
    from pytest_mock import MockerFixture

    from anndata.abc import CSCDataset, CSRDataset

    Idx = slice | int | NDArray[np.integer] | NDArray[np.bool_]


subset_func2 = subset_func


M = 50
N = 50


@pytest.fixture
def zarr_metadata_key():
    return ".zarray" if ad.settings.zarr_write_format == 2 else "zarr.json"


@pytest.fixture
def zarr_separator():
    return "" if ad.settings.zarr_write_format == 2 else "/c"


@pytest.fixture
def ondisk_equivalent_adata(
    tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]
) -> tuple[AnnData, AnnData, AnnData, AnnData]:
    csr_path = tmp_path / f"csr.{diskfmt}"
    csc_path = tmp_path / f"csc.{diskfmt}"
    dense_path = tmp_path / f"dense.{diskfmt}"

    write = lambda x, pth, **kwargs: getattr(x, f"write_{diskfmt}")(pth, **kwargs)

    csr_mem = ad.AnnData(X=sparse.random(M, N, format="csr", density=0.1))
    csc_mem = ad.AnnData(X=csr_mem.X.tocsc())
    dense_mem = ad.AnnData(X=csr_mem.X.toarray())

    write(csr_mem, csr_path)
    write(csc_mem, csc_path)
    # write(csr_mem, dense_path, as_dense="X")
    write(dense_mem, dense_path)
    if diskfmt == "h5ad":
        csr_disk = ad.read_h5ad(csr_path, backed="r")
        csc_disk = ad.read_h5ad(csc_path, backed="r")
        dense_disk = ad.read_h5ad(dense_path, backed="r")
    else:

        def read_zarr_backed(path):
            path = str(path)

            f = zarr.open(path, mode="r")

            # Read with handling for backwards compat
            def callback(func, elem_name, elem, iospec):
                if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
                    return AnnData(
                        **{
                            k: read_dispatched(v, callback)
                            for k, v in dict(elem).items()
                        }
                    )
                if iospec.encoding_type in {"csc_matrix", "csr_matrix"}:
                    return sparse_dataset(elem)
                return func(elem)

            adata = read_dispatched(f, callback=callback)

            return adata

        csr_disk = read_zarr_backed(csr_path)
        csc_disk = read_zarr_backed(csc_path)
        dense_disk = read_zarr_backed(dense_path)

    return csr_mem, csr_disk, csc_disk, dense_disk


@pytest.mark.parametrize(
    "empty_mask", [[], np.zeros(M, dtype=bool)], ids=["empty_list", "empty_bool_mask"]
)
def test_empty_backed_indexing(
    ondisk_equivalent_adata: tuple[AnnData, AnnData, AnnData, AnnData],
    empty_mask,
):
    csr_mem, csr_disk, csc_disk, _ = ondisk_equivalent_adata

    assert_equal(csr_mem.X[empty_mask], csr_disk.X[empty_mask])
    assert_equal(csr_mem.X[:, empty_mask], csc_disk.X[:, empty_mask])

    # The following do not work because of https://github.com/scipy/scipy/issues/19919
    # Our implementation returns a (0,0) sized matrix but scipy does (1,0).

    # assert_equal(csr_mem.X[empty_mask, empty_mask], csr_disk.X[empty_mask, empty_mask])
    # assert_equal(csr_mem.X[empty_mask, empty_mask], csc_disk.X[empty_mask, empty_mask])


def test_backed_indexing(
    ondisk_equivalent_adata: tuple[AnnData, AnnData, AnnData, AnnData],
    subset_func,
    subset_func2,
):
    csr_mem, csr_disk, csc_disk, dense_disk = ondisk_equivalent_adata

    obs_idx = subset_func(csr_mem.obs_names)
    var_idx = subset_func2(csr_mem.var_names)

    assert_equal(csr_mem[obs_idx, var_idx].X, csr_disk[obs_idx, var_idx].X)
    assert_equal(csr_mem[obs_idx, var_idx].X, csc_disk[obs_idx, var_idx].X)
    assert_equal(csr_mem.X[...], csc_disk.X[...])
    assert_equal(csr_mem[obs_idx, :].X, dense_disk[obs_idx, :].X)
    assert_equal(csr_mem[obs_idx].X, csr_disk[obs_idx].X)
    assert_equal(csr_mem[:, var_idx].X, dense_disk[:, var_idx].X)


def test_backed_ellipsis_indexing(
    ondisk_equivalent_adata: tuple[AnnData, AnnData, AnnData, AnnData],
    ellipsis_index: tuple[EllipsisType | slice, ...] | EllipsisType,
    equivalent_ellipsis_index: tuple[slice, slice],
):
    csr_mem, csr_disk, csc_disk, _ = ondisk_equivalent_adata

    assert_equal(csr_mem.X[equivalent_ellipsis_index], csr_disk.X[ellipsis_index])
    assert_equal(csr_mem.X[equivalent_ellipsis_index], csc_disk.X[ellipsis_index])


def make_randomized_mask(size: int) -> np.ndarray:
    randomized_mask = np.zeros(size, dtype=bool)
    inds = np.random.choice(size, 20, replace=False)
    inds.sort()
    for i in range(0, len(inds) - 1, 2):
        randomized_mask[inds[i] : inds[i + 1]] = True
    return randomized_mask


def make_alternating_mask(size: int, step: int) -> np.ndarray:
    mask_alternating = np.ones(size, dtype=bool)
    for i in range(0, size, step):  # 5 is too low to trigger new behavior
        mask_alternating[i] = False
    return mask_alternating


# non-random indices, with alternating one false and n true
make_alternating_mask_5 = partial(make_alternating_mask, step=5)
make_alternating_mask_15 = partial(make_alternating_mask, step=15)


def make_one_group_mask(size: int) -> np.ndarray:
    one_group_mask = np.zeros(size, dtype=bool)
    one_group_mask[1 : size // 2] = True
    return one_group_mask


def make_one_elem_mask(size: int) -> np.ndarray:
    one_elem_mask = np.zeros(size, dtype=bool)
    one_elem_mask[size // 4] = True
    return one_elem_mask


# test behavior from https://github.com/scverse/anndata/pull/1233
@pytest.mark.parametrize(
    ("make_bool_mask", "should_trigger_optimization"),
    [
        (make_randomized_mask, None),
        (make_alternating_mask_15, True),
        (make_alternating_mask_5, False),
        (make_one_group_mask, True),
        (make_one_elem_mask, False),
    ],
    ids=["randomized", "alternating_15", "alternating_5", "one_group", "one_elem"],
)
def test_consecutive_bool(
    mocker: MockerFixture,
    ondisk_equivalent_adata: tuple[AnnData, AnnData, AnnData, AnnData],
    make_bool_mask: Callable[[int], np.ndarray],
    should_trigger_optimization: bool | None,
):
    """Tests for optimization from https://github.com/scverse/anndata/pull/1233

    Parameters
    ----------
    mocker
        Mocker object
    ondisk_equivalent_adata
        AnnData objects with sparse X for testing
    make_bool_mask
        Function for creating a boolean mask.
    should_trigger_optimization
        Whether or not a given mask should trigger the optimized behavior.
    """
    _, csr_disk, csc_disk, _ = ondisk_equivalent_adata
    mask = make_bool_mask(csr_disk.shape[0])

    # indexing needs to be on `X` directly to trigger the optimization.

    # `_normalize_indices`, which is used by `AnnData`, converts bools to ints with `np.where`
    from anndata._core import sparse_dataset

    spy = mocker.spy(sparse_dataset, "get_compressed_vectors_for_slices")
    assert_equal(csr_disk.X[mask, :], csr_disk.X[np.where(mask)])
    if should_trigger_optimization is not None:
        assert (
            spy.call_count == 1 if should_trigger_optimization else not spy.call_count
        )
    assert_equal(csc_disk.X[:, mask], csc_disk.X[:, np.where(mask)[0]])
    if should_trigger_optimization is not None:
        assert (
            spy.call_count == 2 if should_trigger_optimization else not spy.call_count
        )
    assert_equal(csr_disk[mask, :], csr_disk[np.where(mask)])
    if should_trigger_optimization is not None:
        assert (
            spy.call_count == 3 if should_trigger_optimization else not spy.call_count
        )
    subset = csc_disk[:, mask]
    assert_equal(subset, csc_disk[:, np.where(mask)[0]])
    if should_trigger_optimization is not None:
        assert (
            spy.call_count == 4 if should_trigger_optimization else not spy.call_count
        )
    if should_trigger_optimization is not None and not csc_disk.isbacked:
        size = subset.shape[1]
        if should_trigger_optimization:
            subset_subset_mask = np.ones(size).astype("bool")
            subset_subset_mask[size // 2] = False
        else:
            subset_subset_mask = make_one_elem_mask(size)
        assert_equal(
            subset[:, subset_subset_mask], subset[:, np.where(subset_subset_mask)[0]]
        )
        assert (
            spy.call_count == 5 if should_trigger_optimization else not spy.call_count
        ), f"Actual count: {spy.call_count}"


@pytest.mark.parametrize(
    ("sparse_format", "append_method"),
    [
        pytest.param(sparse.csr_matrix, sparse.vstack),
        pytest.param(sparse.csc_matrix, sparse.hstack),
        pytest.param(sparse.csr_array, sparse.vstack),
        pytest.param(sparse.csc_array, sparse.hstack),
    ],
)
def test_dataset_append_memory(
    tmp_path: Path,
    sparse_format: Callable[[ArrayLike], CSMatrix],
    append_method: Callable[[list[CSMatrix]], CSMatrix],
    diskfmt: Literal["h5ad", "zarr"],
):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a = sparse_format(sparse.random(100, 100))
    b = sparse_format(sparse.random(100, 100))
    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")
    ad.io.write_elem(f, "mtx", a)
    diskmtx = sparse_dataset(f["mtx"])

    diskmtx.append(b)
    fromdisk = diskmtx.to_memory()

    frommem = append_method([a, b])

    assert_equal(fromdisk, frommem)


def test_append_array_cache_bust(tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a = sparse.random(100, 100, format="csr")
    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")
    ad.io.write_elem(f, "mtx", a)
    ad.io.write_elem(f, "mtx_2", a)
    diskmtx = sparse_dataset(f["mtx"])
    old_array_shapes = {}
    array_names = ["indptr", "indices", "data"]
    for name in array_names:
        old_array_shapes[name] = getattr(diskmtx, f"_{name}").shape
    diskmtx.append(sparse_dataset(f["mtx_2"]))
    for name in array_names:
        assert old_array_shapes[name] != getattr(diskmtx, f"_{name}").shape


@pytest.mark.parametrize("sparse_format", [sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize(
    ("subset_func", "subset_func2"),
    product(
        [
            ad.tests.helpers.array_subset,
            ad.tests.helpers.slice_subset,
            ad.tests.helpers.array_int_subset,
            ad.tests.helpers.array_bool_subset,
        ],
        repeat=2,
    ),
)
def test_read_array(
    tmp_path: Path,
    sparse_format: Callable[[ArrayLike], CSMatrix],
    diskfmt: Literal["h5ad", "zarr"],
    subset_func,
    subset_func2,
):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a = sparse_format(sparse.random(100, 100))
    obs_idx = subset_func(np.arange(100))
    var_idx = subset_func2(np.arange(100))
    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")
    ad.io.write_elem(f, "mtx", a)
    diskmtx = sparse_dataset(f["mtx"])
    ad.settings.use_sparse_array_on_read = True
    assert issubclass(type(diskmtx[obs_idx, var_idx]), CSArray)
    ad.settings.use_sparse_array_on_read = False
    assert issubclass(type(diskmtx[obs_idx, var_idx]), CSMatrix)


@pytest.mark.parametrize(
    ("sparse_format", "append_method"),
    [
        pytest.param(sparse.csr_matrix, sparse.vstack),
        pytest.param(sparse.csc_matrix, sparse.hstack),
    ],
)
def test_dataset_append_disk(
    tmp_path: Path,
    sparse_format: Callable[[ArrayLike], CSMatrix],
    append_method: Callable[[list[CSMatrix]], CSMatrix],
    diskfmt: Literal["h5ad", "zarr"],
):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a = sparse_format(sparse.random(10, 10))
    b = sparse_format(sparse.random(10, 10))

    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")
    ad.io.write_elem(f, "a", a)
    ad.io.write_elem(f, "b", b)
    a_disk = sparse_dataset(f["a"])
    b_disk = sparse_dataset(f["b"])

    a_disk.append(b_disk)
    fromdisk = a_disk.to_memory()

    frommem = append_method([a, b])

    assert_equal(fromdisk, frommem)


@pytest.mark.parametrize("sparse_format", [sparse.csr_matrix, sparse.csc_matrix])
def test_lazy_array_cache(
    tmp_path: Path, sparse_format: Callable[[ArrayLike], CSMatrix], zarr_metadata_key
):
    elems = {"indptr", "indices", "data"}
    path = tmp_path / "test.zarr"
    a = sparse_format(sparse.random(10, 10))
    f = open_write_group(path, mode="a")
    ad.io.write_elem(f, "X", a)
    store = AccessTrackingStore(path)
    for elem in elems:
        store.initialize_key_trackers([f"X/{elem}"])
    f = open_write_group(store, mode="a")
    a_disk = sparse_dataset(f["X"])
    a_disk[:1]
    a_disk[3:5]
    a_disk[6:7]
    a_disk[8:9]
    # one each for .zarray and actual access
    # see https://github.com/zarr-developers/zarr-python/discussions/2760 for why 4
    assert store.get_access_count("X/indptr") == 2 if is_zarr_v2() else 4
    for elem_not_indptr in elems - {"indptr"}:
        assert (
            sum(
                zarr_metadata_key in key_accessed
                for key_accessed in store.get_accessed_keys(f"X/{elem_not_indptr}")
            )
            == 1
        )


Kind = Literal["slice", "int", "array", "mask"]


def mk_idx_kind(idx: Sequence[int], *, kind: Kind, l: int) -> Idx | None:
    """Convert sequence of consecutive integers (e.g. range with step=1) into different kinds of indexing."""
    if kind == "slice":
        start = idx[0] if idx[0] > 0 else None
        if len(idx) == 1:
            return slice(start, idx[0] + 1)
        if all(np.diff(idx) == 1):
            stop = idx[-1] + 1 if idx[-1] < l - 1 else None
            return slice(start, stop)
    if kind == "int" and len(idx) == 1:
        return idx[0]
    if kind == "array":
        return np.asarray(idx)
    if kind == "mask":
        return np.isin(np.arange(l), idx)
    return None


def idify(x: object) -> str:
    if isinstance(x, slice):
        start, stop = ("" if s is None else str(s) for s in (x.start, x.stop))
        return f"{start}:{stop}" + (f":{x.step}" if x.step not in (1, None) else "")
    return str(x)


def width_idx_kinds(
    *idxs: tuple[Sequence[int], Idx, Sequence[str]], l: int
) -> Generator[ParameterSet, None, None]:
    """Convert major (first) index into various identical kinds of indexing."""
    for (idx_maj_raw, idx_min, exp), maj_kind in product(idxs, get_args(Kind)):
        if (idx_maj := mk_idx_kind(idx_maj_raw, kind=maj_kind, l=l)) is None:
            continue
        id_ = "-".join(map(idify, [idx_maj_raw, idx_min, maj_kind]))
        yield pytest.param(idx_maj, idx_min, exp, id=id_)


@pytest.mark.parametrize("sparse_format", [sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize(
    ("idx_maj", "idx_min", "exp"),
    width_idx_kinds(
        (
            [0],
            slice(None, None),
            ["X/data/{zarr_metadata_key}", "X/data{zarr_separator}/0"],
        ),
        (
            [0],
            slice(None, 3),
            ["X/data/{zarr_metadata_key}", "X/data{zarr_separator}/0"],
        ),
        (
            [3, 4, 5],
            slice(None, None),
            [
                "X/data/{zarr_metadata_key}",
                "X/data{zarr_separator}/3",
                "X/data{zarr_separator}/4",
                "X/data{zarr_separator}/5",
            ],
        ),
        l=10,
    ),
)
@pytest.mark.parametrize(
    "open_func",
    [
        sparse_dataset,
        lambda x: read_elem_lazy(
            x, chunks=(1, -1) if x.attrs["encoding-type"] == "csr_matrix" else (-1, 1)
        ),
    ],
    ids=["sparse_dataset", "read_elem_lazy"],
)
def test_data_access(
    tmp_path: Path,
    sparse_format: Callable[[ArrayLike], CSMatrix],
    idx_maj: Idx,
    idx_min: Idx,
    exp: list[str],
    open_func: Callable[[ZarrGroup], CSRDataset | CSCDataset | DaskArray],
    zarr_metadata_key,
    zarr_separator,
):
    exp = [
        e.format(zarr_metadata_key=zarr_metadata_key, zarr_separator=zarr_separator)
        for e in exp
    ]
    path = tmp_path / "test.zarr"
    a = sparse_format(np.eye(10, 10))
    f = open_write_group(path, mode="a")
    ad.io.write_elem(f, "X", a)
    data = f["X/data"][...]
    del f["X/data"]
    # chunk one at a time to count properly
    zarr.array(
        data,
        store=path / "X" / "data",
        chunks=(1,),
        zarr_format=ad.settings.zarr_write_format,
    )
    store = AccessTrackingStore(path)
    store.initialize_key_trackers(["X/data"])
    f = zarr.open_group(store)
    a_disk = AnnData(X=open_func(f["X"]))
    subset = a_disk[idx_maj, idx_min] if a.format == "csr" else a_disk[idx_min, idx_maj]
    if isinstance(subset.X, DaskArray):
        subset.X.compute(scheduler="single-threaded")
    # zarr v2 fetches all and not just metadata for that node in 3.X.X python package
    # TODO: https://github.com/zarr-developers/zarr-python/discussions/2760
    if ad.settings.zarr_write_format == 2 and not is_zarr_v2():
        exp = [*exp, "X/data/.zgroup", "X/data/.zattrs"]

    assert store.get_access_count("X/data") == len(exp), store.get_accessed_keys(
        "X/data"
    )
    # dask access order is not guaranteed so need to sort
    assert sorted(store.get_accessed_keys("X/data")) == sorted(exp)


@pytest.mark.parametrize(
    ("sparse_format", "a_shape", "b_shape"),
    [
        pytest.param("csr", (100, 100), (100, 200)),
        pytest.param("csc", (100, 100), (200, 100)),
    ],
)
def test_wrong_shape(
    tmp_path: Path,
    sparse_format: Literal["csr", "csc"],
    a_shape: tuple[int, int],
    b_shape: tuple[int, int],
    diskfmt: Literal["h5ad", "zarr"],
):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    a_mem = sparse.random(*a_shape, format=sparse_format)
    b_mem = sparse.random(*b_shape, format=sparse_format)

    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")

    ad.io.write_elem(f, "a", a_mem)
    ad.io.write_elem(f, "b", b_mem)
    a_disk = sparse_dataset(f["a"])
    b_disk = sparse_dataset(f["b"])

    with pytest.raises(AssertionError):
        a_disk.append(b_disk)


def test_reset_group(tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]):
    path = tmp_path / "test.zarr"
    base = sparse.random(100, 100, format="csr")

    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")

    ad.io.write_elem(f, "base", base)
    disk_mtx = sparse_dataset(f["base"])
    with pytest.raises(AttributeError):
        disk_mtx.group = f


def test_wrong_formats(tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    base = sparse.random(100, 100, format="csr")

    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")

    ad.io.write_elem(f, "base", base)
    disk_mtx = sparse_dataset(f["base"])
    pre_checks = disk_mtx.to_memory()

    with pytest.raises(ValueError, match="must have same format"):
        disk_mtx.append(sparse.random(100, 100, format="csc"))
    with pytest.raises(ValueError, match="must have same format"):
        disk_mtx.append(sparse.random(100, 100, format="coo"))
    with pytest.raises(NotImplementedError):
        disk_mtx.append(np.random.random((100, 100)))
    if isinstance(f, ZarrGroup) and not is_zarr_v2():
        data = np.random.random((100, 100))
        disk_dense = f.create_array("dense", shape=(100, 100), dtype=data.dtype)
        disk_dense[...] = data
    else:
        disk_dense = f.create_dataset(
            "dense", data=np.random.random((100, 100)), shape=(100, 100)
        )
    with pytest.raises(NotImplementedError):
        disk_mtx.append(disk_dense)

    post_checks = disk_mtx.to_memory()

    # Check nothing changed
    assert not np.any((pre_checks != post_checks).toarray())


def test_anndata_sparse_compat(tmp_path: Path, diskfmt: Literal["h5ad", "zarr"]):
    path = tmp_path / f"test.{diskfmt.replace('ad', '')}"
    base = sparse.random(100, 100, format="csr")

    f = open_write_group(path, mode="a") if diskfmt == "zarr" else h5py.File(path, "a")

    ad.io.write_elem(f, "/", base)
    adata = ad.AnnData(sparse_dataset(f["/"]))
    assert_equal(adata.X, base)


def test_backed_sizeof(
    ondisk_equivalent_adata: tuple[AnnData, AnnData, AnnData, AnnData],
):
    csr_mem, csr_disk, csc_disk, _ = ondisk_equivalent_adata

    assert csr_mem.__sizeof__() == csr_disk.__sizeof__(with_disk=True)
    assert csr_mem.__sizeof__() == csc_disk.__sizeof__(with_disk=True)
    assert csr_disk.__sizeof__(with_disk=True) == csc_disk.__sizeof__(with_disk=True)
    assert csr_mem.__sizeof__() > csr_disk.__sizeof__()
    assert csr_mem.__sizeof__() > csc_disk.__sizeof__()


@pytest.mark.parametrize(
    "group_fn",
    [
        pytest.param(lambda _: zarr.group(), id="zarr"),
        pytest.param(lambda p: h5py.File(p / "test.h5", mode="a"), id="h5py"),
    ],
)
@pytest.mark.parametrize(
    "sparse_class",
    [
        sparse.csr_matrix,
        pytest.param(
            sparse.csr_array,
            marks=[pytest.mark.skip(reason="scipy bug causes view to be allocated")],
        ),
    ],
)
def test_append_overflow_check(group_fn, sparse_class, tmp_path):
    group = group_fn(tmp_path)
    typemax_int32 = np.iinfo(np.int32).max
    orig_mtx = sparse_class(np.ones((1, 1), dtype=bool))
    # Minimally allocating new matrix
    new_mtx = sparse_class(
        (
            np.broadcast_to(True, typemax_int32 - 1),  # noqa: FBT003
            np.broadcast_to(np.int32(1), typemax_int32 - 1),
            [0, typemax_int32 - 1],
        ),
        shape=(1, 2),
    )

    ad.io.write_elem(group, "mtx", orig_mtx)
    backed = sparse_dataset(group["mtx"])

    # Checking for correct caching behaviour
    backed._indptr  # noqa: B018

    with pytest.raises(
        OverflowError,
        match=r"This array was written with a 32 bit intptr, but is now large.*",
    ):
        backed.append(new_mtx)

    # Check for any modification
    assert_equal(backed, orig_mtx)
