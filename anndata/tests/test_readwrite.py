from contextlib import contextmanager
from importlib.util import find_spec
from os import PathLike
from pathlib import Path
import re
from string import ascii_letters
import warnings

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix, csc_matrix
import zarr

import anndata as ad
from anndata._io.specs.registry import IORegistryError
from anndata.compat import _read_attr, DaskArray

from anndata.tests.helpers import (
    gen_adata,
    assert_equal,
    as_dense_dask_array,
    pytest_8_raises,
)

HERE = Path(__file__).parent


# ------------------------------------------------------------------------------
# Some test data
# ------------------------------------------------------------------------------


X_sp = csr_matrix([[1, 0, 0], [3, 0, 0], [5, 6, 0], [0, 0, 0], [0, 0, 0]])

X_list = [[1, 0], [3, 0], [5, 6]]  # data matrix of shape n_obs x n_vars

obs_dict = dict(  # annotation of observations / rows
    row_names=["name1", "name2", "name3"],  # row annotation
    oanno1=["cat1", "cat2", "cat2"],  # categorical annotation
    oanno1b=["cat1", "cat1", "cat1"],  # categorical annotation with one category
    oanno1c=["cat1", "cat1", np.nan],  # categorical annotation with a missing value
    oanno2=["o1", "o2", "o3"],  # string annotation
    oanno3=[2.1, 2.2, 2.3],  # float annotation
    oanno4=[3.3, 1.1, 2.2],  # float annotation
)

var_dict = dict(  # annotation of variables / columns
    vanno1=[3.1, 3.2],
    vanno2=["cat1", "cat1"],  # categorical annotation
    vanno3=[2.1, 2.2],  # float annotation
    vanno4=[3.3, 1.1],  # float annotation
)

uns_dict = dict(  # unstructured annotation
    oanno1_colors=["#000000", "#FFFFFF"],
    uns2=["some annotation"],
    uns3="another annotation",
    uns4=dict(
        a=1,
        b=[2, 3],
        c="4",
        d=["some", "strings"],
        e=np.ones(5),
        f=np.int32(7),
        g=[1, np.float32(2.5)],
    ),
)


@pytest.fixture(params=[{}, dict(compression="gzip")])
def dataset_kwargs(request):
    return request.param


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.fixture
def rw(backing_h5ad):
    M, N = 100, 101
    orig = gen_adata((M, N))
    orig.write(backing_h5ad)
    curr = ad.read(backing_h5ad)
    return curr, orig


@pytest.fixture(params=[np.uint8, np.int32, np.int64, np.float32, np.float64])
def dtype(request):
    return request.param


diskfmt2 = diskfmt


# ------------------------------------------------------------------------------
# The test functions
# ------------------------------------------------------------------------------


@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_roundtrip(typ, tmp_path, diskfmt, diskfmt2):
    tmpdir = Path(tmp_path)
    pth1 = tmpdir / f"first.{diskfmt}"
    write1 = lambda x: getattr(x, f"write_{diskfmt}")(pth1)
    read1 = lambda: getattr(ad, f"read_{diskfmt}")(pth1)
    pth2 = tmpdir / f"second.{diskfmt2}"
    write2 = lambda x: getattr(x, f"write_{diskfmt2}")(pth2)
    read2 = lambda: getattr(ad, f"read_{diskfmt2}")(pth2)

    adata1 = ad.AnnData(typ(X_list), obs=obs_dict, var=var_dict, uns=uns_dict)
    write1(adata1)
    adata2 = read1()
    write2(adata2)
    adata3 = read2()

    assert_equal(adata2, adata1)
    assert_equal(adata3, adata1)
    assert_equal(adata2, adata1)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_h5ad(tmp_path, typ, dataset_kwargs, backing_h5ad):
    mid_pth = tmp_path / "mid.h5ad"

    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    assert not isinstance(adata_src.obs["oanno1"].dtype, pd.CategoricalDtype)
    adata_src.raw = adata_src
    adata_src.write(backing_h5ad, **dataset_kwargs)

    adata_mid = ad.read(backing_h5ad)
    adata_mid.write(mid_pth, **dataset_kwargs)

    adata = ad.read_h5ad(mid_pth)
    assert isinstance(adata.obs["oanno1"].dtype, pd.CategoricalDtype)
    assert not isinstance(adata.obs["oanno2"].dtype, pd.CategoricalDtype)
    assert adata.obs.index.tolist() == ["name1", "name2", "name3"]
    assert adata.obs["oanno1"].cat.categories.tolist() == ["cat1", "cat2"]
    assert adata.obs["oanno1c"].cat.categories.tolist() == ["cat1"]
    assert isinstance(adata.raw.var["vanno2"].dtype, pd.CategoricalDtype)
    pd.testing.assert_frame_equal(adata.obs, adata_src.obs)
    pd.testing.assert_frame_equal(adata.var, adata_src.var)
    assert_equal(adata.var.index, adata_src.var.index)
    assert adata.var.index.dtype == adata_src.var.index.dtype

    assert isinstance(adata_src.raw.X, (type(adata.raw.X), DaskArray))
    assert isinstance(
        adata_src.uns["uns4"]["c"], (type(adata.uns["uns4"]["c"]), DaskArray)
    )
    assert isinstance(adata_src.varm, (type(adata.varm), DaskArray))

    assert_equal(adata.raw.X, adata_src.raw.X)
    pd.testing.assert_frame_equal(adata.raw.var, adata_src.raw.var)
    assert isinstance(adata.uns["uns4"]["a"], (int, np.integer))
    assert isinstance(adata_src.uns["uns4"]["a"], (int, np.integer))
    assert_equal(adata, adata_src)


@pytest.mark.skipif(not find_spec("zarr"), reason="Zarr is not installed")
@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_zarr(typ, tmp_path):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.raw = adata_src
    assert not isinstance(adata_src.obs["oanno1"].dtype, pd.CategoricalDtype)
    adata_src.write_zarr(tmp_path / "test_zarr_dir", chunks=True)

    adata = ad.read_zarr(tmp_path / "test_zarr_dir")
    assert isinstance(adata.obs["oanno1"].dtype, pd.CategoricalDtype)
    assert not isinstance(adata.obs["oanno2"].dtype, pd.CategoricalDtype)
    assert adata.obs.index.tolist() == ["name1", "name2", "name3"]
    assert adata.obs["oanno1"].cat.categories.tolist() == ["cat1", "cat2"]
    assert adata.obs["oanno1c"].cat.categories.tolist() == ["cat1"]
    assert isinstance(adata.raw.var["vanno2"].dtype, pd.CategoricalDtype)
    pd.testing.assert_frame_equal(adata.obs, adata_src.obs)
    pd.testing.assert_frame_equal(adata.var, adata_src.var)
    assert_equal(adata.var.index, adata_src.var.index)
    assert adata.var.index.dtype == adata_src.var.index.dtype

    # Dev. Note:
    # either load as same type or load the convert DaskArray to array
    # since we tested if assigned types and loaded types are DaskArray
    # this would also work if they work
    assert isinstance(adata_src.raw.X, (type(adata.raw.X), DaskArray))
    assert isinstance(
        adata_src.uns["uns4"]["c"], (type(adata.uns["uns4"]["c"]), DaskArray)
    )
    assert isinstance(adata_src.varm, (type(adata.varm), DaskArray))

    assert_equal(adata.raw.X, adata_src.raw.X)
    assert_equal(adata.raw.var, adata_src.raw.var)
    assert isinstance(adata.uns["uns4"]["a"], (int, np.integer))
    assert isinstance(adata_src.uns["uns4"]["a"], (int, np.integer))
    assert_equal(adata, adata_src)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_maintain_X_dtype(typ, backing_h5ad):
    X = typ(X_list).astype("int8")
    adata_src = ad.AnnData(X)
    adata_src.write(backing_h5ad)

    adata = ad.read(backing_h5ad)
    assert adata.X.dtype == adata_src.X.dtype


def test_read_write_maintain_obsmvarm_dtypes(rw):
    curr, orig = rw

    assert type(orig.obsm["array"]) is type(curr.obsm["array"])
    assert np.all(orig.obsm["array"] == curr.obsm["array"])
    assert np.all(orig.varm["array"] == curr.varm["array"])
    assert type(orig.obsm["sparse"]) is type(curr.obsm["sparse"])
    assert not np.any((orig.obsm["sparse"] != curr.obsm["sparse"]).toarray())
    assert not np.any((orig.varm["sparse"] != curr.varm["sparse"]).toarray())
    assert type(orig.obsm["df"]) is type(curr.obsm["df"])
    assert np.all(orig.obsm["df"] == curr.obsm["df"])
    assert np.all(orig.varm["df"] == curr.varm["df"])


def test_maintain_layers(rw):
    curr, orig = rw

    assert type(orig.layers["array"]) is type(curr.layers["array"])
    assert np.all(orig.layers["array"] == curr.layers["array"])
    assert type(orig.layers["sparse"]) is type(curr.layers["sparse"])
    assert not np.any((orig.layers["sparse"] != curr.layers["sparse"]).toarray())


@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_h5ad_one_dimension(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_one = adata_src[:, 0].copy()
    adata_one.write(backing_h5ad)
    adata = ad.read(backing_h5ad)
    assert adata.shape == (3, 1)
    assert_equal(adata, adata_one)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, as_dense_dask_array])
def test_readwrite_backed(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.filename = backing_h5ad  # change to backed mode
    adata_src.write()

    adata = ad.read(backing_h5ad)
    assert isinstance(adata.obs["oanno1"].dtype, pd.CategoricalDtype)
    assert not isinstance(adata.obs["oanno2"].dtype, pd.CategoricalDtype)
    assert adata.obs.index.tolist() == ["name1", "name2", "name3"]
    assert adata.obs["oanno1"].cat.categories.tolist() == ["cat1", "cat2"]
    assert_equal(adata, adata_src)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csc_matrix])
def test_readwrite_equivalent_h5ad_zarr(tmp_path, typ):
    h5ad_pth = tmp_path / "adata.h5ad"
    zarr_pth = tmp_path / "adata.zarr"

    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ)
    adata.raw = adata

    adata.write_h5ad(h5ad_pth)
    adata.write_zarr(zarr_pth)
    from_h5ad = ad.read_h5ad(h5ad_pth)
    from_zarr = ad.read_zarr(zarr_pth)

    assert_equal(from_h5ad, from_zarr, exact=True)


@contextmanager
def store_context(path: Path):
    if path.suffix == ".zarr":
        store = zarr.open(path, "r+")
    else:
        file = h5py.File(path, "r+")
        store = file["/"]
    yield store
    if "file" in locals():
        file.close()


@pytest.mark.parametrize(
    ["name", "read", "write"],
    [
        ("adata.h5ad", ad.read_h5ad, ad.AnnData.write_h5ad),
        ("adata.zarr", ad.read_zarr, ad.AnnData.write_zarr),
    ],
)
def test_read_full_io_error(tmp_path, name, read, write):
    adata = gen_adata((4, 3))
    path = tmp_path / name
    write(adata, path)
    with store_context(path) as store:
        store["obs"].attrs["encoding-type"] = "invalid"
    with pytest_8_raises(
        IORegistryError,
        match=r"raised while reading key '/obs'",
    ) as exc_info:
        read(path)
    assert re.search(
        r"No read method registered for IOSpec\(encoding_type='invalid', encoding_version='0.2.0'\)",
        str(exc_info.value),
    )


@pytest.mark.parametrize(
    "compression,compression_opts",
    [
        (None, None),
        ("lzf", None),
        ("gzip", None),
        ("gzip", 8),
    ],
)
def test_hdf5_compression_opts(tmp_path, compression, compression_opts):
    # https://github.com/scverse/anndata/issues/497
    pth = Path(tmp_path) / "adata.h5ad"
    adata = gen_adata((10, 8))
    kwargs = {}
    if compression is not None:
        kwargs["compression"] = compression
    if compression_opts is not None:
        kwargs["compression_opts"] = compression_opts
    not_compressed = []

    adata.write_h5ad(pth, **kwargs)

    def check_compressed(key, value):
        if isinstance(value, h5py.Dataset) and value.shape != ():
            if compression is not None and value.compression != compression:
                not_compressed.append(key)
            elif (
                compression_opts is not None
                and value.compression_opts != compression_opts
            ):
                not_compressed.append(key)

    with h5py.File(pth) as f:
        f.visititems(check_compressed)

    if not_compressed:
        msg = "\n\t".join(not_compressed)
        raise AssertionError(f"These elements were not compressed correctly:\n\t{msg}")

    assert_equal(adata, ad.read_h5ad(pth))


def test_zarr_compression(tmp_path):
    from numcodecs import Blosc

    pth = str(Path(tmp_path) / "adata.zarr")
    adata = gen_adata((10, 8))
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    not_compressed = []

    ad._io.write_zarr(pth, adata, compressor=compressor)

    def check_compressed(key, value):
        if isinstance(value, zarr.Array) and value.shape != ():
            if value.compressor != compressor:
                not_compressed.append(key)

    with zarr.open(str(pth), "r") as f:
        f.visititems(check_compressed)

    if not_compressed:
        msg = "\n\t".join(not_compressed)
        raise AssertionError(f"These elements were not compressed correctly:\n\t{msg}")

    assert_equal(adata, ad.read_zarr(pth))


def test_changed_obs_var_names(tmp_path, diskfmt):
    filepth = tmp_path / f"test.{diskfmt}"

    orig = gen_adata((10, 10))
    orig.obs_names.name = "obs"
    orig.var_names.name = "var"
    modified = orig.copy()
    modified.obs_names.name = "cells"
    modified.var_names.name = "genes"

    getattr(orig, f"write_{diskfmt}")(filepth)
    read = getattr(ad, f"read_{diskfmt}")(filepth)

    assert_equal(orig, read, exact=True)
    assert orig.var.index.name == "var"
    assert read.obs.index.name == "obs"
    with pytest.raises(AssertionError):
        assert_equal(orig, modified, exact=True)
    with pytest.raises(AssertionError):
        assert_equal(read, modified, exact=True)


@pytest.mark.skipif(not find_spec("loompy"), reason="Loompy is not installed")
@pytest.mark.parametrize("typ", [np.array, csr_matrix])
@pytest.mark.parametrize("obsm_mapping", [{}, dict(X_composed=["oanno3", "oanno4"])])
@pytest.mark.parametrize("varm_mapping", [{}, dict(X_composed2=["vanno3", "vanno4"])])
def test_readwrite_loom(typ, obsm_mapping, varm_mapping, tmp_path):
    X = typ(X_list)
    obs_dim = "meaningful_obs_dim_name"
    var_dim = "meaningful_var_dim_name"
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.obs_names.name = obs_dim
    adata_src.var_names.name = var_dim
    adata_src.obsm["X_a"] = np.zeros((adata_src.n_obs, 2))
    adata_src.varm["X_b"] = np.zeros((adata_src.n_vars, 3))

    adata_src.write_loom(tmp_path / "test.loom", write_obsm_varm=True)

    adata = ad.read_loom(
        tmp_path / "test.loom",
        sparse=typ is csr_matrix,
        obsm_mapping=obsm_mapping,
        obs_names=obs_dim,
        varm_mapping=varm_mapping,
        var_names=var_dim,
        cleanup=True,
    )

    if isinstance(X, np.ndarray):
        assert np.allclose(adata.X, X)
    else:
        # TODO: this should not be necessary
        assert np.allclose(adata.X.toarray(), X.toarray())
    assert "X_a" in adata.obsm_keys() and adata.obsm["X_a"].shape[1] == 2
    assert "X_b" in adata.varm_keys() and adata.varm["X_b"].shape[1] == 3
    # as we called with `cleanup=True`
    assert "oanno1b" in adata.uns["loom-obs"]
    assert "vanno2" in adata.uns["loom-var"]
    for k, v in obsm_mapping.items():
        assert k in adata.obsm_keys() and adata.obsm[k].shape[1] == len(v)
    for k, v in varm_mapping.items():
        assert k in adata.varm_keys() and adata.varm[k].shape[1] == len(v)
    assert adata.obs_names.name == obs_dim
    assert adata.var_names.name == var_dim


@pytest.mark.skipif(not find_spec("loompy"), reason="Loompy is not installed")
def test_readloom_deprecations(tmp_path):
    loom_pth = tmp_path / "test.loom"
    adata_src = gen_adata((5, 10), obsm_types=[np.ndarray], varm_types=[np.ndarray])
    adata_src.write_loom(loom_pth, write_obsm_varm=True)

    # obsm_names -> obsm_mapping
    obsm_mapping = {"df": adata_src.obs.columns}
    with pytest.warns(FutureWarning):
        depr_result = ad.read_loom(loom_pth, obsm_names=obsm_mapping)
    actual_result = ad.read_loom(loom_pth, obsm_mapping=obsm_mapping)
    assert_equal(actual_result, depr_result)
    with pytest.raises(ValueError, match="ambiguous"):
        ad.read_loom(loom_pth, obsm_mapping=obsm_mapping, obsm_names=obsm_mapping)

    # varm_names -> varm_mapping
    varm_mapping = {"df": adata_src.var.columns}
    with pytest.warns(FutureWarning):
        depr_result = ad.read_loom(loom_pth, varm_names=varm_mapping)
    actual_result = ad.read_loom(loom_pth, varm_mapping=varm_mapping)
    assert_equal(actual_result, depr_result)
    with pytest.raises(ValueError, match="ambiguous"):
        ad.read_loom(loom_pth, varm_mapping=varm_mapping, varm_names=varm_mapping)

    # positional -> keyword
    with pytest.warns(FutureWarning, match="sparse"):
        depr_result = ad.read_loom(loom_pth, True)
    actual_result = ad.read_loom(loom_pth, sparse=True)
    assert type(depr_result.X) == type(actual_result.X)


def test_read_csv():
    adata = ad.read_csv(HERE / "adata.csv")
    assert adata.obs_names.tolist() == ["r1", "r2", "r3"]
    assert adata.var_names.tolist() == ["c1", "c2"]
    assert adata.X.tolist() == X_list


def test_read_tsv_strpath():
    adata = ad.read_text(str(HERE / "adata-comments.tsv"), "\t")
    assert adata.obs_names.tolist() == ["r1", "r2", "r3"]
    assert adata.var_names.tolist() == ["c1", "c2"]
    assert adata.X.tolist() == X_list


def test_read_tsv_iter():
    with (HERE / "adata-comments.tsv").open() as f:
        adata = ad.read_text(f, "\t")
    assert adata.obs_names.tolist() == ["r1", "r2", "r3"]
    assert adata.var_names.tolist() == ["c1", "c2"]
    assert adata.X.tolist() == X_list


@pytest.mark.parametrize("typ", [np.array, csr_matrix])
def test_write_csv(typ, tmp_path):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata.write_csvs(tmp_path / "test_csv_dir", skip_data=False)


@pytest.mark.parametrize("typ", [np.array, csr_matrix])
def test_write_csv_view(typ, tmp_path):
    # https://github.com/scverse/anndata/issues/401
    import hashlib

    def md5_path(pth: PathLike) -> bytes:
        checksum = hashlib.md5()
        with open(pth, "rb") as f:
            while True:
                buf = f.read(checksum.block_size * 100)
                if not buf:
                    break
                checksum.update(buf)
        return checksum.digest()

    def hash_dir_contents(dir: Path) -> "dict[str, bytes]":
        root_pth = str(dir)
        return {
            str(k)[len(root_pth) :]: md5_path(k) for k in dir.rglob("*") if k.is_file()
        }

    adata = ad.AnnData(typ(X_list), obs=obs_dict, var=var_dict, uns=uns_dict)

    # Test writing a view
    view_pth = tmp_path / "test_view_csv_dir"
    copy_pth = tmp_path / "test_copy_csv_dir"
    adata[::2].write_csvs(view_pth, skip_data=False)
    adata[::2].copy().write_csvs(copy_pth, skip_data=False)

    assert hash_dir_contents(view_pth) == hash_dir_contents(copy_pth)


@pytest.mark.parametrize(
    ["read", "write", "name"],
    [
        pytest.param(ad.read_h5ad, ad._io.write_h5ad, "test_empty.h5ad"),
        pytest.param(
            ad.read_loom,
            ad._io.write_loom,
            "test_empty.loom",
            marks=pytest.mark.xfail(reason="Loom can’t handle 0×0 matrices"),
        ),
        pytest.param(ad.read_zarr, ad._io.write_zarr, "test_empty.zarr"),
        pytest.param(
            ad.read_zarr,
            ad._io.write_zarr,
            "test_empty.zip",
            marks=pytest.mark.xfail(reason="Zarr zip storage doesn’t seem to work…"),
        ),
    ],
)
def test_readwrite_hdf5_empty(read, write, name, tmp_path):
    adata = ad.AnnData(uns=dict(empty=np.array([], dtype=float)))
    write(tmp_path / name, adata)
    ad_read = read(tmp_path / name)
    assert ad_read.uns["empty"].shape == (0,)


def test_read_excel():
    adata = ad.read_excel(HERE / "data/excel.xlsx", "Sheet1", dtype=int)
    assert adata.X.tolist() == X_list


def test_read_umi_tools():
    adata = ad.read_umi_tools(HERE / "data/umi_tools.tsv.gz")
    assert adata.obs_names.name == "cell"
    assert adata.var_names.name == "gene"
    assert adata.shape == (2, 13)
    assert "ENSG00000070404.9" in adata.var_names
    assert set(adata.obs_names) == {"ACAAGG", "TTCACG"}


def test_write_categorical(tmp_path, diskfmt):
    adata_pth = tmp_path / f"adata.{diskfmt}"
    orig = ad.AnnData(
        obs=pd.DataFrame(
            dict(
                cat1=["a", "a", "b", np.nan, np.nan],
                cat2=pd.Categorical(["a", "a", "b", np.nan, np.nan]),
            )
        ),
    )
    getattr(orig, f"write_{diskfmt}")(adata_pth)
    curr = getattr(ad, f"read_{diskfmt}")(adata_pth)
    assert np.all(orig.obs.notna() == curr.obs.notna())
    assert np.all(orig.obs.stack().dropna() == curr.obs.stack().dropna())


def test_write_categorical_index(tmp_path, diskfmt):
    adata_pth = tmp_path / f"adata.{diskfmt}"
    orig = ad.AnnData(
        uns={"df": pd.DataFrame({}, index=pd.Categorical(list("aabcd")))},
    )
    getattr(orig, f"write_{diskfmt}")(adata_pth)
    curr = getattr(ad, f"read_{diskfmt}")(adata_pth)
    # Also covered by next assertion, but checking this value specifically
    pd.testing.assert_index_equal(
        orig.uns["df"].index, curr.uns["df"].index, exact=True
    )
    assert_equal(orig, curr, exact=True)


def test_dataframe_reserved_columns(tmp_path, diskfmt):
    reserved = ("_index",)
    adata_pth = tmp_path / f"adata.{diskfmt}"
    orig = ad.AnnData(
        obs=pd.DataFrame(index=np.arange(5)), var=pd.DataFrame(index=np.arange(5))
    )
    for colname in reserved:
        to_write = orig.copy()
        to_write.obs[colname] = np.ones(5)
        with pytest.raises(ValueError) as exc_info:
            getattr(to_write, f"write_{diskfmt}")(adata_pth)
        assert colname in str(exc_info.value)
    for colname in reserved:
        to_write = orig.copy()
        to_write.varm["df"] = pd.DataFrame(
            {colname: list("aabcd")}, index=to_write.var_names
        )
        with pytest.raises(ValueError) as exc_info:
            getattr(to_write, f"write_{diskfmt}")(adata_pth)
        assert colname in str(exc_info.value)


def test_write_large_categorical(tmp_path, diskfmt):
    M = 30_000
    N = 1000
    ls = np.array(list(ascii_letters))

    def random_cats(n):
        cats = {
            "".join(np.random.choice(ls, np.random.choice(range(5, 30))))
            for _ in range(n)
        }
        while len(cats) < n:  # For the rare case that there’s duplicates
            cats |= random_cats(n - len(cats))
        return cats

    cats = np.array(sorted(random_cats(10_000)))
    adata_pth = tmp_path / f"adata.{diskfmt}"
    n_cats = len(np.unique(cats))
    orig = ad.AnnData(
        csr_matrix(([1], ([0], [0])), shape=(M, N)),
        obs=dict(
            cat1=cats[np.random.choice(n_cats, M)],
            cat2=pd.Categorical.from_codes(np.random.choice(n_cats, M), cats),
        ),
    )
    getattr(orig, f"write_{diskfmt}")(adata_pth)
    curr = getattr(ad, f"read_{diskfmt}")(adata_pth)
    assert_equal(orig, curr)


def test_write_string_types(tmp_path, diskfmt):
    # https://github.com/scverse/anndata/issues/456
    adata_pth = tmp_path / f"adata.{diskfmt}"

    adata = ad.AnnData(
        obs=pd.DataFrame(
            np.ones((3, 2)),
            columns=["a", np.str_("b")],
            index=["a", "b", "c"],
        ),
    )

    write = getattr(adata, f"write_{diskfmt}")
    read = getattr(ad, f"read_{diskfmt}")

    write(adata_pth)
    from_disk = read(adata_pth)

    assert_equal(adata, from_disk)

    adata.obs[b"c"] = np.zeros(3)
    # This should error, and tell you which key is at fault
    with pytest_8_raises(TypeError, match=r"writing key 'obs'") as exc_info:
        write(adata_pth)

    assert str("b'c'") in str(exc_info.value)


@pytest.mark.parametrize(
    "teststring",
    ["teststring", np.asarray(["test1", "test2", "test3"], dtype="object")],
)
@pytest.mark.parametrize("encoding", ["ascii", "utf-8"])
@pytest.mark.parametrize("length", [None, 15])
def test_hdf5_attribute_conversion(tmp_path, teststring, encoding, length):
    with h5py.File(tmp_path / "attributes.h5", "w") as file:
        dset = file.create_dataset("dset", data=np.arange(10))
        attrs = dset.attrs
        attrs.create(
            "string",
            teststring,
            dtype=h5py.h5t.string_dtype(encoding=encoding, length=length),
        )

        assert_equal(teststring, _read_attr(attrs, "string"))


def test_zarr_chunk_X(tmp_path):
    import zarr

    zarr_pth = Path(tmp_path) / "test.zarr"
    adata = gen_adata((100, 100), X_type=np.array)
    adata.write_zarr(zarr_pth, chunks=(10, 10))

    z = zarr.open(str(zarr_pth))  # As of v2.3.2 zarr won’t take a Path
    assert z["X"].chunks == (10, 10)
    from_zarr = ad.read_zarr(zarr_pth)
    assert_equal(from_zarr, adata)


################################
# Round-tripping scanpy datasets
################################

diskfmt2 = diskfmt


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_scanpy_pbmc68k(tmp_path, diskfmt, diskfmt2):
    read1 = lambda pth: getattr(ad, f"read_{diskfmt}")(pth)
    write1 = lambda adata, pth: getattr(adata, f"write_{diskfmt}")(pth)
    read2 = lambda pth: getattr(ad, f"read_{diskfmt2}")(pth)
    write2 = lambda adata, pth: getattr(adata, f"write_{diskfmt2}")(pth)

    filepth1 = tmp_path / f"test1.{diskfmt}"
    filepth2 = tmp_path / f"test2.{diskfmt2}"

    import scanpy as sc

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ad.OldFormatWarning)
        pbmc = sc.datasets.pbmc68k_reduced()

    write1(pbmc, filepth1)
    from_disk1 = read1(filepth1)  # Do we read okay
    write2(from_disk1, filepth2)  # Can we round trip
    from_disk2 = read2(filepth2)

    assert_equal(pbmc, from_disk1)  # Not expected to be exact due to `nan`s
    assert_equal(pbmc, from_disk2)


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_scanpy_krumsiek11(tmp_path, diskfmt):
    filepth = tmp_path / f"test.{diskfmt}"
    import scanpy as sc

    orig = sc.datasets.krumsiek11()
    del orig.uns["highlights"]  # Can’t write int keys
    getattr(orig, f"write_{diskfmt}")(filepth)
    read = getattr(ad, f"read_{diskfmt}")(filepth)

    assert_equal(orig, read, exact=True)


# Checking if we can read legacy zarr files
# TODO: Check how I should add this file to the repo
@pytest.mark.filterwarnings("ignore::anndata.OldFormatWarning")
@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
@pytest.mark.skipif(
    not Path(HERE / "data/pbmc68k_reduced_legacy.zarr.zip").is_file(),
    reason="File not present.",
)
def test_backwards_compat_zarr():
    import scanpy as sc
    import zarr

    pbmc_orig = sc.datasets.pbmc68k_reduced()
    # Old zarr writer couldn’t do sparse arrays
    pbmc_orig.raw._X = pbmc_orig.raw.X.toarray()
    del pbmc_orig.uns["neighbors"]
    # Since these have moved, see PR #337
    del pbmc_orig.obsp["distances"]
    del pbmc_orig.obsp["connectivities"]

    # This was written out with anndata=0.6.22.post1
    zarrpth = HERE / "data/pbmc68k_reduced_legacy.zarr.zip"
    with zarr.ZipStore(zarrpth, mode="r") as z:
        pbmc_zarr = ad.read_zarr(z)

    assert_equal(pbmc_zarr, pbmc_orig)


# TODO: use diskfmt fixture once zarr backend implemented
def test_adata_in_uns(tmp_path, diskfmt):
    pth = tmp_path / f"adatas_in_uns.{diskfmt}"
    read = lambda pth: getattr(ad, f"read_{diskfmt}")(pth)
    write = lambda adata, pth: getattr(adata, f"write_{diskfmt}")(pth)

    orig = gen_adata((4, 5))
    orig.uns["adatas"] = {
        "a": gen_adata((1, 2)),
        "b": gen_adata((12, 8)),
    }
    another_one = gen_adata((2, 5))
    another_one.raw = gen_adata((2, 7))
    orig.uns["adatas"]["b"].uns["another_one"] = another_one

    write(orig, pth)
    curr = read(pth)

    assert_equal(orig, curr)


def test_io_dtype(tmp_path, diskfmt, dtype):
    pth = tmp_path / f"adata_dtype.{diskfmt}"
    read = lambda pth: getattr(ad, f"read_{diskfmt}")(pth)
    write = lambda adata, pth: getattr(adata, f"write_{diskfmt}")(pth)

    orig = ad.AnnData(np.ones((5, 8), dtype=dtype))
    write(orig, pth)
    curr = read(pth)

    assert curr.X.dtype == dtype
