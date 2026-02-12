"""Tests for backing using the `.file` and `.isbacked` attributes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import joblib
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import (
    GEN_ADATA_DASK_ARGS,
    GEN_ADATA_NO_XARRAY_ARGS,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
    subset_func,
)
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal

    from anndata.compat import DaskArray

subset_func2 = subset_func


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------


@pytest.fixture
def adata() -> ad.AnnData:
    X_list = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]  # data matrix of shape n_obs x n_vars
    X = np.array(X_list)
    obs_dict = dict(  # annotation of observations / rows
        row_names=["name1", "name2", "name3"],  # row annotation
        oanno1=["cat1", "cat2", "cat2"],  # categorical annotation
        oanno2=["o1", "o2", "o3"],  # string annotation
        oanno3=[2.1, 2.2, 2.3],  # float annotation
    )
    var_dict = dict(vanno1=[3.1, 3.2, 3.3])  # annotation of variables / columns
    uns_dict = dict(  # unstructured annotation
        oanno1_colors=["#000000", "#FFFFFF"], uns2=["some annotation"]
    )
    return ad.AnnData(
        X,
        obs=obs_dict,
        var=var_dict,
        uns=uns_dict,
        obsm=dict(o1=np.zeros((X.shape[0], 10))),
        varm=dict(v1=np.ones((X.shape[1], 20))),
        layers=dict(float=X.astype(float), sparse=sparse.csr_matrix(X)),
    )


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
)
def mtx_format(
    request,
) -> Callable[
    [np.ndarray], DaskArray | np.ndarray | sparse.csr_array | sparse.csr_matrix
]:
    return request.param


@pytest.fixture(params=[sparse.csr_matrix, sparse.csc_matrix])
def sparse_format(request) -> type[sparse.csr_matrix | sparse.csc_matrix]:
    return request.param


@pytest.fixture(params=["r+", "r", False])
def backed_mode(request) -> Literal["r+", "r", False]:
    return request.param


@pytest.fixture(params=(("X",), ()))
def as_dense(request) -> tuple[str] | tuple:
    return request.param


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


# h5py internally calls `product` on min-versions
@pytest.mark.filterwarnings("ignore:`product` is deprecated as of NumPy 1.25.0")
# TODO: Check to make sure obs, obsm, layers, ... are written and read correctly as well
@pytest.mark.filterwarnings("error")
def test_read_write_X(
    tmp_path: Path,
    mtx_format: Callable[
        [np.ndarray], DaskArray | np.ndarray | sparse.csr_array | sparse.csr_matrix
    ],
    backed_mode: Literal["r+", "r", False],
    *,
    as_dense: tuple[str] | tuple,
):
    orig_pth = tmp_path / "orig.h5ad"
    backed_pth = tmp_path / "backed.h5ad"

    orig = ad.AnnData(mtx_format(asarray(sparse.random(10, 10, format="csr"))))
    orig.write(orig_pth)

    backed = ad.read_h5ad(orig_pth, backed=backed_mode)
    backed.write(backed_pth, as_dense=as_dense)
    backed.file.close()

    from_backed = ad.read_h5ad(backed_pth)
    assert np.all(asarray(orig.X) == asarray(from_backed.X))


def test_backed_view(tmp_path: Path, backed_mode: Literal["r+", "r", False]):
    orig_pth = tmp_path / "orig.h5ad"

    orig = ad.AnnData(sparse.random(100, 10, format="csr"))
    orig.write(orig_pth)

    adata = ad.read_h5ad(orig_pth, backed=backed_mode)

    for i in range(0, adata.shape[0], 10):
        chunk_path = tmp_path / f"chunk_{i}.h5ad"
        adata[i : i + 5].write_h5ad(tmp_path / f"chunk_{i}.h5ad")
        chunk = adata[i : i + 5]
        assert_equal(chunk, ad.read_h5ad(chunk_path))


# this is very similar to the views test
@pytest.mark.filterwarnings("ignore::anndata.ImplicitModificationWarning")
def test_backing(adata: ad.AnnData, tmp_path: Path, backing_h5ad: Path) -> None:
    assert not adata.isbacked

    adata.filename = backing_h5ad
    adata.write()
    assert not adata.file.is_open
    assert adata.isbacked
    assert adata[:, 0].is_view
    assert adata[:, 0].X.tolist() == np.reshape([1, 4, 7], (3, 1)).tolist()
    # this might give us a trouble as the user might not
    # know that the file is open again....
    assert adata.file.is_open

    adata[:2, 0].X = np.array([[0], [0]])
    assert adata[:, 0].X.tolist() == np.reshape([1, 4, 7], (3, 1)).tolist()

    adata_subset = adata[:2, [0, 1]]
    assert adata_subset.is_view
    subset_hash = joblib.hash(adata_subset)

    # cannot set view in backing mode...
    with pytest.raises(ValueError, match=r"pass a filename.*to_memory"):
        adata_subset.obs["foo"] = range(2)
    with pytest.raises(ValueError, match=r"pass a filename.*to_memory"):
        adata_subset.var["bar"] = -12
    with pytest.raises(ValueError, match=r"pass a filename.*to_memory"):
        adata_subset.obsm["o2"] = np.ones((2, 2))
    with pytest.raises(ValueError, match=r"pass a filename.*to_memory"):
        adata_subset.varm["v2"] = np.zeros((2, 2))
    with pytest.raises(ValueError, match=r"pass a filename.*to_memory"):
        adata_subset.layers["float2"] = adata_subset.layers["float"].copy()

    # Things should stay the same after failed operations
    assert subset_hash == joblib.hash(adata_subset)
    assert adata_subset.is_view

    # need to copy first
    adata_subset = adata_subset.copy(tmp_path / "test.subset.h5ad")
    # now transition to actual object
    assert not adata_subset.is_view
    adata_subset.obs["foo"] = range(2)
    assert not adata_subset.is_view
    assert adata_subset.isbacked
    assert adata_subset.obs["foo"].tolist() == list(range(2))

    # save
    adata_subset.write()


def test_backing_copy(adata, tmp_path: Path, backing_h5ad: Path):
    adata.filename = backing_h5ad
    adata.write()

    copypath = tmp_path / "test.copy.h5ad"
    copy = adata.copy(copypath)

    assert adata.filename == backing_h5ad
    assert copy.filename == copypath
    assert adata.isbacked
    assert copy.isbacked


# TODO: Also test updating the backing file inplace
def test_backed_raw(tmp_path: Path):
    backed_pth = tmp_path / "backed.h5ad"
    final_pth = tmp_path / "final.h5ad"
    mem_adata = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    mem_adata.raw = mem_adata
    mem_adata.write(backed_pth)

    backed_adata = ad.read_h5ad(backed_pth, backed="r")
    assert_equal(backed_adata, mem_adata)
    backed_adata.write_h5ad(final_pth)

    final_adata = ad.read_h5ad(final_pth)
    assert_equal(final_adata, mem_adata)


@pytest.mark.parametrize(
    "array_type",
    [
        pytest.param(asarray, id="dense_array"),
        pytest.param(sparse.csr_matrix, id="csr_matrix"),
        pytest.param(sparse.csr_array, id="csr_array"),
    ],
)
def test_backed_raw_subset(
    tmp_path: Path,
    array_type: Callable[
        [np.ndarray], np.ndarray | sparse.csr_array | sparse.csr_matrix
    ],
    subset_func: Callable[[ad.AnnData], ad.AnnData],
    subset_func2: Callable[[ad.AnnData], ad.AnnData],
):
    backed_pth = tmp_path / "backed.h5ad"
    final_pth = tmp_path / "final.h5ad"
    mem_adata = gen_adata((10, 10), X_type=array_type, **GEN_ADATA_NO_XARRAY_ARGS)
    mem_adata.raw = mem_adata
    obs_idx = subset_func(mem_adata.obs_names)
    var_idx = subset_func2(mem_adata.var_names)
    mem_adata.write(backed_pth)

    ### Backed view has same values as in memory view ###
    backed_adata = ad.read_h5ad(backed_pth, backed="r")
    backed_v = backed_adata[obs_idx, var_idx]
    assert backed_v.is_view
    mem_v = mem_adata[obs_idx, var_idx]

    # Value equivalent
    assert_equal(mem_v, backed_v)
    # Type and value equivalent
    assert_equal(mem_v.copy(), backed_v.to_memory(copy=True), exact=True)
    assert backed_v.is_view
    assert backed_v.isbacked

    ### Write from backed view ###
    backed_v.write_h5ad(final_pth)
    final_adata = ad.read_h5ad(final_pth)

    assert_equal(mem_v, final_adata)
    assert_equal(final_adata, backed_v.to_memory())  # assert loading into memory


@pytest.mark.parametrize(
    "array_type",
    [
        pytest.param(asarray, id="dense_array"),
        pytest.param(sparse.csr_matrix, id="csr_matrix"),
        pytest.param(as_dense_dask_array, id="dask_array"),
    ],
)
def test_to_memory_full(
    tmp_path: Path,
    array_type: Callable[[np.ndarray], np.ndarray | DaskArray | sparse.csr_matrix],
):
    backed_pth = tmp_path / "backed.h5ad"
    mem_adata = gen_adata((15, 10), X_type=array_type, **GEN_ADATA_DASK_ARGS)
    mem_adata.raw = gen_adata((15, 12), X_type=array_type, **GEN_ADATA_DASK_ARGS)
    mem_adata.write_h5ad(backed_pth, compression="lzf")

    backed_adata = ad.read_h5ad(backed_pth, backed="r")
    assert_equal(mem_adata, backed_adata.to_memory())

    # Test that raw can be removed
    del backed_adata.raw
    del mem_adata.raw
    assert_equal(mem_adata, backed_adata.to_memory())


def test_double_index(adata: ad.AnnData, backing_h5ad: Path):
    adata.filename = backing_h5ad
    with pytest.raises(ValueError, match=r"cannot make a view of a view"):
        # no view of view of backed object currently
        adata[:2][:, 0]

    # close backing file
    adata.write()


def test_return_to_memory_mode(adata: ad.AnnData, backing_h5ad: Path):
    bdata = adata.copy()
    adata.filename = backing_h5ad
    assert adata.isbacked

    adata.filename = None
    assert not adata.isbacked

    assert adata.X is not None

    # make sure the previous file had been properly closed
    # when setting `adata.filename = None`
    # if it hadn’t the following line would throw an error
    bdata.filename = backing_h5ad
    # close the file
    bdata.filename = None


def test_backed_modification(adata: ad.AnnData, backing_h5ad: Path):
    adata.X[:, 1] = 0  # Make it a little sparse
    adata.X = sparse.csr_matrix(adata.X)
    assert not adata.isbacked

    # While this currently makes the file backed, it doesn’t write it as sparse
    adata.filename = backing_h5ad
    adata.write()
    assert not adata.file.is_open
    assert adata.isbacked

    adata.X[0, [0, 2]] = 10
    adata.X[1, [0, 2]] = [11, 12]
    adata.X[2, 1] = 13  # If it were written as sparse, this should fail

    assert adata.isbacked

    assert np.all(adata.X[0, :] == np.array([10, 0, 10]))
    assert np.all(adata.X[1, :] == np.array([11, 0, 12]))
    assert np.all(adata.X[2, :] == np.array([7, 13, 9]))


def test_backed_modification_sparse(
    adata: ad.AnnData,
    backing_h5ad: Path,
    sparse_format: type[sparse.csr_matrix | sparse.csc_matrix],
):
    adata.X[:, 1] = 0  # Make it a little sparse
    adata.X = sparse_format(adata.X)
    orig = adata.X.copy()
    assert not adata.isbacked

    adata.write(backing_h5ad)
    adata = ad.read_h5ad(backing_h5ad, backed="r+")

    assert adata.filename == backing_h5ad
    assert adata.isbacked

    # Does not modify backed store
    adata[0, [0, 2]].X = np.array([[10, 10]])
    np.testing.assert_equal(orig.toarray(), adata.X[...].toarray())


# TODO: Work around h5py not supporting this
# def test_backed_view_modification(adata, backing_h5ad):
#     adata.write(backing_h5ad)
#     backed_adata = ad.read_h5ad(backing_h5ad, backed=True)

#     backed_view = backed_adata[[1, 2], :]
#     backed_view.X = 0

#     assert np.all(backed_adata.X[:3, :] == 0)


# TODO: Implement
# def test_backed_view_modification_sparse(adata, backing_h5ad, sparse_format):
#     adata[:, 1] = 0  # Make it a little sparse
#     adata.X = sparse_format(adata.X)
#     adata.write(backing_h5ad)
#     backed_adata = ad.read_h5ad(backing_h5ad, backed=True)

#     backed_view = backed_adata[[1,2], :]
#     backed_view.X = 0
#     assert np.all(backed_adata.X[[1,2], :] == 0)


@pytest.mark.parametrize(
    ("obs_idx", "var_idx"),
    [
        pytest.param(np.array([0, 1, 2]), np.array([1, 2]), id="no_dupes"),
        pytest.param(np.array([0, 1, 0, 2]), slice(None), id="1d_dupes"),
        pytest.param(np.array([0, 1, 0, 2]), np.array([1, 2, 1]), id="2d_dupes"),
    ],
)
def test_backed_duplicate_indices(tmp_path, obs_idx, var_idx):
    """Test that backed HDF5 datasets handle duplicate indices correctly."""
    backed_pth = tmp_path / "backed.h5ad"

    # Create test data
    mem_adata = gen_adata((6, 4), X_type=asarray, **GEN_ADATA_NO_XARRAY_ARGS)
    mem_adata.write(backed_pth)

    # Load backed data
    backed_adata = ad.read_h5ad(backed_pth, backed="r")

    # Test the indexing
    mem_result_multi = mem_adata[obs_idx, var_idx]
    backed_result_multi = backed_adata[obs_idx, var_idx]
    assert_equal(mem_result_multi, backed_result_multi)


@pytest.fixture
def h5py_test_data(tmp_path):
    """Create test HDF5 file with dataset for _safe_fancy_index_h5py tests."""
    import h5py

    h5_path = tmp_path / "test_dataset.h5"
    test_data = np.arange(24).reshape(6, 4)  # 6x4 matrix

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("test", data=test_data)

    return h5_path, test_data


@pytest.mark.parametrize(
    ("indices", "description"),
    [
        pytest.param((np.array([0, 1, 0, 2]),), "single_dimension_with_duplicates"),
        pytest.param(
            (np.array([0, 1, 2]), np.array([1, 2])), "multi_dimensional_no_duplicates"
        ),
        pytest.param(
            (np.array([0, 1, 0, 2]), np.array([1, 2])),
            "multi_dimensional_duplicates_first_dim",
        ),
        pytest.param(
            (np.array([0, 1, 2]), np.array([1, 2, 1])),
            "multi_dimensional_duplicates_second_dim",
        ),
        pytest.param(
            (np.array([0, 1, 0]), np.array([1, 2, 1])),
            "multi_dimensional_duplicates_both_dims",
        ),
        pytest.param(
            (np.array([True, False, True, False, False, True]),), "boolean_arrays"
        ),
        pytest.param((np.array([0, 1, 0]), slice(1, 3)), "mixed_indexing_with_slices"),
        pytest.param(
            (np.array([0, 1, 0]), [1, 2]), "mixed_indexing_with_slices_and_lists"
        ),
        pytest.param((np.array([3, 1, 3, 0, 1]),), "unsorted_indices_with_duplicates"),
    ],
)
def test_safe_fancy_index_h5py_function(h5py_test_data, indices, description):
    """Test the _safe_fancy_index_h5py function directly with various indexing patterns."""
    import h5py

    from anndata._core.index import _safe_fancy_index_h5py

    h5_path, test_data = h5py_test_data

    with h5py.File(h5_path, "r") as f:
        dataset = f["test"]

        # Get result from the function
        result = _safe_fancy_index_h5py(dataset, indices)

    # Calculate expected result using NumPy
    if isinstance(indices, tuple) and len(indices) > 1:
        # Multi-dimensional case - use np.ix_ for fancy indexing
        if isinstance(indices[1], slice):
            # Handle mixed case with slice
            expected = test_data[
                np.ix_(indices[0], np.arange(indices[1].start, indices[1].stop))
            ]
        else:
            expected = test_data[np.ix_(*indices)]
    else:
        # Single dimensional case
        expected = test_data[indices]

    # Assert arrays are equal
    np.testing.assert_array_equal(
        result, expected, err_msg=f"Failed for test case: {description}"
    )
