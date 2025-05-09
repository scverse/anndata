from __future__ import annotations

import re
import warnings
from contextlib import contextmanager
from functools import partial
from importlib.util import find_spec
from pathlib import Path
from string import ascii_letters
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from numba.core.errors import NumbaDeprecationWarning
from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix

import anndata as ad
from anndata._io.specs.registry import IORegistryError
from anndata._io.zarr import open_write_group
from anndata.compat import (
    CSArray,
    CSMatrix,
    DaskArray,
    ZarrArray,
    ZarrGroup,
    _read_attr,
    is_zarr_v2,
)
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
)

if TYPE_CHECKING:
    from typing import Literal

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
    uns5=None,
)


@pytest.fixture(params=[{}, dict(compression="gzip")])
def dataset_kwargs(request):
    return request.param


@pytest.fixture
def rw(backing_h5ad):
    M, N = 100, 101
    orig = gen_adata((M, N), **GEN_ADATA_NO_XARRAY_ARGS)
    orig.write(backing_h5ad)
    curr = ad.read_h5ad(backing_h5ad)
    return curr, orig


@pytest.fixture(params=[np.uint8, np.int32, np.int64, np.float32, np.float64])
def dtype(request):
    return request.param


# ------------------------------------------------------------------------------
# The test functions
# ------------------------------------------------------------------------------


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csr_array, as_dense_dask_array])
def test_readwrite_roundtrip(typ, tmp_path, diskfmt, diskfmt2):
    pth1 = tmp_path / f"first.{diskfmt}"
    write1 = lambda x: getattr(x, f"write_{diskfmt}")(pth1)
    read1 = lambda: getattr(ad, f"read_{diskfmt}")(pth1)
    pth2 = tmp_path / f"second.{diskfmt2}"
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


def test_readwrite_roundtrip_async(tmp_path):
    import asyncio

    async def _do_test():
        zarr_path = tmp_path / "first.zarr"

        adata1 = ad.AnnData(
            csr_matrix(X_list), obs=obs_dict, var=var_dict, uns=uns_dict
        )
        adata1.write_zarr(zarr_path)
        adata2 = ad.read_zarr(zarr_path)

        assert_equal(adata2, adata1)

    # This test ensures our file i/o never calls `asyncio.run` internally
    asyncio.run(_do_test())


@pytest.mark.parametrize("storage", ["h5ad", "zarr"])
@pytest.mark.parametrize("typ", [np.array, csr_matrix, csr_array, as_dense_dask_array])
def test_readwrite_kitchensink(tmp_path, storage, typ, backing_h5ad, dataset_kwargs):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    assert not isinstance(adata_src.obs["oanno1"].dtype, pd.CategoricalDtype)
    adata_src.raw = adata_src.copy()

    if storage == "h5ad":
        adata_src.write(backing_h5ad, **dataset_kwargs)
        adata_mid = ad.read_h5ad(backing_h5ad)
        adata_mid.write(tmp_path / "mid.h5ad", **dataset_kwargs)
        adata = ad.read_h5ad(tmp_path / "mid.h5ad")
    else:
        adata_src.write_zarr(tmp_path / "test_zarr_dir")
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
    if isinstance(adata_src.raw.X, CSArray):
        assert isinstance(adata.raw.X, CSMatrix)
    else:
        assert isinstance(adata_src.raw.X, type(adata.raw.X) | DaskArray)
    assert isinstance(
        adata_src.uns["uns4"]["c"], type(adata.uns["uns4"]["c"]) | DaskArray
    )
    assert isinstance(adata_src.varm, type(adata.varm) | DaskArray)

    assert_equal(adata.raw.X, adata_src.raw.X)
    pd.testing.assert_frame_equal(adata.raw.var, adata_src.raw.var)
    assert isinstance(adata.uns["uns4"]["a"], int | np.integer)
    assert isinstance(adata_src.uns["uns4"]["a"], int | np.integer)
    assert_equal(adata, adata_src)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csr_array, as_dense_dask_array])
def test_readwrite_maintain_X_dtype(typ, backing_h5ad):
    X = typ(X_list).astype("int8")
    adata_src = ad.AnnData(X)
    adata_src.write(backing_h5ad)

    adata = ad.read_h5ad(backing_h5ad)
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


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csr_array, as_dense_dask_array])
def test_readwrite_h5ad_one_dimension(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_one = adata_src[:, 0].copy()
    adata_one.write(backing_h5ad)
    adata = ad.read_h5ad(backing_h5ad)
    assert adata.shape == (3, 1)
    assert_equal(adata, adata_one)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csr_array, as_dense_dask_array])
def test_readwrite_backed(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.filename = backing_h5ad  # change to backed mode
    adata_src.write()

    adata = ad.read_h5ad(backing_h5ad)
    assert isinstance(adata.obs["oanno1"].dtype, pd.CategoricalDtype)
    assert not isinstance(adata.obs["oanno2"].dtype, pd.CategoricalDtype)
    assert adata.obs.index.tolist() == ["name1", "name2", "name3"]
    assert adata.obs["oanno1"].cat.categories.tolist() == ["cat1", "cat2"]
    assert_equal(adata, adata_src)


@pytest.mark.parametrize(
    "typ", [np.array, csr_matrix, csc_matrix, csr_array, csc_array]
)
def test_readwrite_equivalent_h5ad_zarr(tmp_path, typ):
    h5ad_pth = tmp_path / "adata.h5ad"
    zarr_pth = tmp_path / "adata.zarr"

    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ, **GEN_ADATA_NO_XARRAY_ARGS)
    adata.raw = adata.copy()

    adata.write_h5ad(h5ad_pth)
    adata.write_zarr(zarr_pth)
    from_h5ad = ad.read_h5ad(h5ad_pth)
    from_zarr = ad.read_zarr(zarr_pth)

    assert_equal(from_h5ad, from_zarr, exact=True)


@contextmanager
def store_context(path: Path):
    if path.suffix == ".zarr":
        store = open_write_group(path, mode="r+")
    else:
        file = h5py.File(path, "r+")
        store = file["/"]
    yield store
    if "file" in locals():
        file.close()


@pytest.mark.parametrize(
    ("name", "read", "write"),
    [
        ("adata.h5ad", ad.read_h5ad, ad.AnnData.write_h5ad),
        ("adata.zarr", ad.read_zarr, ad.AnnData.write_zarr),
    ],
)
def test_read_full_io_error(tmp_path, name, read, write):
    adata = gen_adata((4, 3), **GEN_ADATA_NO_XARRAY_ARGS)
    path = tmp_path / name
    write(adata, path)
    with store_context(path) as store:
        if not is_zarr_v2() and isinstance(store, ZarrGroup):
            # see https://github.com/zarr-developers/zarr-python/issues/2716 for the issue
            # with re-opening without syncing attributes explicitly
            # TODO: Having to fully specify attributes to not override fixed in zarr v3.0.5
            # See https://github.com/zarr-developers/zarr-python/pull/2870
            store["obs"].update_attributes(
                {**dict(store["obs"].attrs), "encoding-type": "invalid"}
            )
            zarr.consolidate_metadata(store.store)
        else:
            store["obs"].attrs["encoding-type"] = "invalid"

    with pytest.raises(
        IORegistryError,
        match=r"raised while reading key 'obs'.*from /$",
    ) as exc_info:
        read(path)
    assert re.search(
        r"No read method registered for IOSpec\(encoding_type='invalid', encoding_version='0.2.0'\)",
        str(exc_info.value),
    )


@pytest.mark.parametrize(
    ("compression", "compression_opts"),
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
    adata = gen_adata((10, 8), **GEN_ADATA_NO_XARRAY_ARGS)
    kwargs = {}
    if compression is not None:
        kwargs["compression"] = compression
    if compression_opts is not None:
        kwargs["compression_opts"] = compression_opts
    not_compressed = []

    adata.write_h5ad(pth, **kwargs)

    def check_compressed(key, value):
        if not isinstance(value, h5py.Dataset) or value.shape == ():
            return
        if (compression is not None and value.compression != compression) or (
            compression_opts is not None and value.compression_opts != compression_opts
        ):
            not_compressed.append(key)

    with h5py.File(pth) as f:
        f.visititems(check_compressed)

    if not_compressed:
        sep = "\n\t"
        msg = (
            f"These elements were not compressed correctly:{sep}"
            f"{sep.join(not_compressed)}"
        )
        raise AssertionError(msg)

    expected = ad.read_h5ad(pth)
    assert_equal(adata, expected)


@pytest.mark.parametrize("zarr_write_format", [2, 3])
def test_zarr_compression(tmp_path, zarr_write_format):
    ad.settings.zarr_write_format = zarr_write_format
    pth = str(Path(tmp_path) / "adata.zarr")
    adata = gen_adata((10, 8), **GEN_ADATA_NO_XARRAY_ARGS)
    if zarr_write_format == 2 or is_zarr_v2():
        from numcodecs import Blosc

        compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    else:
        from zarr.codecs import BloscCodec

        # Typesize is forced to be 1 so that the codecs always match on the roundtrip.
        # Otherwise this value would vary depending on the datatype.
        # See github.com/zarr-developers/numcodecs/pull/713 for a related issue/explanation.
        # In practice, you would never want to set this parameter.
        compressor = BloscCodec(
            cname="zstd", clevel=3, shuffle="bitshuffle", typesize=1
        )
    not_compressed = []

    ad.io.write_zarr(pth, adata, compressor=compressor)

    def check_compressed(value, key):
        if not isinstance(value, ZarrArray) or value.shape == ():
            return None
        (read_compressor,) = value.compressors
        if zarr_write_format == 2:
            if read_compressor != compressor:
                not_compressed.append(key)
            return None
        if read_compressor.to_dict() != compressor.to_dict():
            not_compressed.append(key)

    if is_zarr_v2():
        with zarr.open(str(pth), "r") as f:
            f.visititems(check_compressed)
    else:
        f = zarr.open(str(pth), mode="r")
        for key, value in f.members(max_depth=None):
            check_compressed(value, key)

    if not_compressed:
        sep = "\n\t"
        msg = (
            f"These elements were not compressed correctly:{sep}"
            f"{sep.join(not_compressed)}"
        )
        raise AssertionError(msg)

    expected = ad.read_zarr(pth)
    assert_equal(adata, expected)


def test_changed_obs_var_names(tmp_path, diskfmt):
    filepth = tmp_path / f"test.{diskfmt}"

    orig = gen_adata((10, 10), **GEN_ADATA_NO_XARRAY_ARGS)
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

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
        # loompy uses “is” for ints
        warnings.filterwarnings("ignore", category=SyntaxWarning)
        warnings.filterwarnings(
            "ignore",
            message=r"datetime.datetime.utcnow\(\) is deprecated",
            category=DeprecationWarning,
        )
        adata_src.write_loom(tmp_path / "test.loom", write_obsm_varm=True)

    adata = ad.io.read_loom(
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
    assert "X_a" in adata.obsm_keys()
    assert adata.obsm["X_a"].shape[1] == 2
    assert "X_b" in adata.varm_keys()
    assert adata.varm["X_b"].shape[1] == 3
    # as we called with `cleanup=True`
    assert "oanno1b" in adata.uns["loom-obs"]
    assert "vanno2" in adata.uns["loom-var"]
    for k, v in obsm_mapping.items():
        assert k in adata.obsm_keys()
        assert adata.obsm[k].shape[1] == len(v)
    for k, v in varm_mapping.items():
        assert k in adata.varm_keys()
        assert adata.varm[k].shape[1] == len(v)
    assert adata.obs_names.name == obs_dim
    assert adata.var_names.name == var_dim


@pytest.mark.skipif(not find_spec("loompy"), reason="Loompy is not installed")
def test_readloom_deprecations(tmp_path):
    loom_pth = tmp_path / "test.loom"
    adata_src = gen_adata((5, 10), obsm_types=[np.ndarray], varm_types=[np.ndarray])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
        warnings.filterwarnings(
            "ignore",
            message=r"datetime.datetime.utcnow\(\) is deprecated",
            category=DeprecationWarning,
        )
        adata_src.write_loom(loom_pth, write_obsm_varm=True)

    # obsm_names -> obsm_mapping
    obsm_mapping = {"df": adata_src.obs.columns}
    with pytest.warns(FutureWarning):
        depr_result = ad.io.read_loom(loom_pth, obsm_names=obsm_mapping)
    actual_result = ad.io.read_loom(loom_pth, obsm_mapping=obsm_mapping)
    assert_equal(actual_result, depr_result)
    with pytest.raises(ValueError, match=r"ambiguous"), pytest.warns(FutureWarning):
        ad.io.read_loom(loom_pth, obsm_mapping=obsm_mapping, obsm_names=obsm_mapping)

    # varm_names -> varm_mapping
    varm_mapping = {"df": adata_src.var.columns}
    with pytest.warns(FutureWarning):
        depr_result = ad.io.read_loom(loom_pth, varm_names=varm_mapping)
    actual_result = ad.io.read_loom(loom_pth, varm_mapping=varm_mapping)
    assert_equal(actual_result, depr_result)
    with pytest.raises(ValueError, match=r"ambiguous"), pytest.warns(FutureWarning):
        ad.io.read_loom(loom_pth, varm_mapping=varm_mapping, varm_names=varm_mapping)

    # positional -> keyword
    with pytest.warns(FutureWarning, match=r"sparse"):
        depr_result = ad.io.read_loom(loom_pth, True)  # noqa: FBT003
    actual_result = ad.io.read_loom(loom_pth, sparse=True)
    assert type(depr_result.X) == type(actual_result.X)


def test_read_csv():
    adata = ad.io.read_csv(HERE / "data" / "adata.csv")
    assert adata.obs_names.tolist() == ["r1", "r2", "r3"]
    assert adata.var_names.tolist() == ["c1", "c2"]
    assert adata.X.tolist() == X_list


def test_read_tsv_strpath():
    adata = ad.io.read_text(str(HERE / "data" / "adata-comments.tsv"), "\t")
    assert adata.obs_names.tolist() == ["r1", "r2", "r3"]
    assert adata.var_names.tolist() == ["c1", "c2"]
    assert adata.X.tolist() == X_list


def test_read_tsv_iter():
    with (HERE / "data" / "adata-comments.tsv").open() as f:
        adata = ad.io.read_text(f, "\t")
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

    def md5_path(pth: Path) -> bytes:
        checksum = hashlib.md5()
        with pth.open("rb") as f:
            while True:
                buf = f.read(checksum.block_size * 100)
                if not buf:
                    break
                checksum.update(buf)
        return checksum.digest()

    def hash_dir_contents(dir: Path) -> dict[str, bytes]:
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
    ("read", "write", "name"),
    [
        pytest.param(ad.read_h5ad, ad.io.write_h5ad, "test_empty.h5ad"),
        pytest.param(
            ad.io.read_loom,
            ad.io.write_loom,
            "test_empty.loom",
            marks=pytest.mark.xfail(reason="Loom can’t handle 0×0 matrices"),
        ),
        pytest.param(ad.read_zarr, ad.io.write_zarr, "test_empty.zarr"),
    ],
)
def test_readwrite_empty(read, write, name, tmp_path):
    adata = ad.AnnData(uns=dict(empty=np.array([], dtype=float)))
    write(tmp_path / name, adata)
    ad_read = read(tmp_path / name)
    assert ad_read.uns["empty"].shape == (0,)


def test_read_excel():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"datetime.datetime.utcnow\(\) is deprecated",
            category=DeprecationWarning,
        )
        adata = ad.io.read_excel(HERE / "data/excel.xlsx", "Sheet1", dtype=int)
    assert adata.X.tolist() == X_list


def test_read_umi_tools():
    adata = ad.io.read_umi_tools(HERE / "data/umi_tools.tsv.gz")
    assert adata.obs_names.name == "cell"
    assert adata.var_names.name == "gene"
    assert adata.shape == (2, 13)
    assert "ENSG00000070404.9" in adata.var_names
    assert set(adata.obs_names) == {"ACAAGG", "TTCACG"}


@pytest.mark.parametrize("s2c", [True, False], ids=["str2cat", "preserve"])
def test_write_categorical(
    *, tmp_path: Path, diskfmt: Literal["h5ad", "zarr"], s2c: bool
) -> None:
    with ad.settings.override(allow_write_nullable_strings=True):
        adata_pth = tmp_path / f"adata.{diskfmt}"
        obs = dict(
            str=pd.array(["a", "a", "b", pd.NA, pd.NA], dtype="string"),
            cat=pd.Categorical(["a", "a", "b", np.nan, np.nan]),
            **(dict(obj=["a", "a", "b", np.nan, np.nan]) if s2c else {}),
        )
        orig = ad.AnnData(obs=pd.DataFrame(obs))
        getattr(orig, f"write_{diskfmt}")(
            adata_pth, convert_strings_to_categoricals=s2c
        )
        curr: ad.AnnData = getattr(ad, f"read_{diskfmt}")(adata_pth)
        assert np.all(orig.obs.notna() == curr.obs.notna())
        assert np.all(orig.obs.stack().dropna() == curr.obs.stack().dropna())
        assert curr.obs["str"].dtype == ("category" if s2c else "string")
        assert curr.obs["cat"].dtype == "category"


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


@pytest.mark.parametrize("colname", ["_index"])
@pytest.mark.parametrize("attr", ["obs", "varm_df"])
def test_dataframe_reserved_columns(tmp_path, diskfmt, colname, attr):
    adata_pth = tmp_path / f"adata.{diskfmt}"
    orig = ad.AnnData(
        obs=pd.DataFrame(index=np.arange(5)), var=pd.DataFrame(index=np.arange(5))
    )

    to_write = orig.copy()
    if attr == "obs":
        to_write.obs[colname] = np.ones(5)
    elif attr == "varm_df":
        to_write.varm["df"] = pd.DataFrame(
            {colname: list("aabcd")}, index=to_write.var_names
        )
    else:
        pytest.fail(f"Unexpected attr: {attr}")
    with pytest.raises(ValueError, match=rf"{colname}.*reserved name"):
        getattr(to_write, f"write_{diskfmt}")(adata_pth)


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


def test_write_string_type_error(tmp_path, diskfmt):
    adata = ad.AnnData(obs=dict(obs_names=list("abc")))
    adata.obs[b"c"] = np.zeros(3)

    # This should error, and tell you which key is at fault
    with pytest.raises(TypeError, match=r"writing key 'obs'") as exc_info:
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"adata.{diskfmt}")

    assert "b'c'" in str(exc_info.value)


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
    adata = gen_adata((100, 100), X_type=np.array, **GEN_ADATA_NO_XARRAY_ARGS)
    adata.write_zarr(zarr_pth, chunks=(10, 10))

    z = zarr.open(str(zarr_pth))  # As of v2.3.2 zarr won’t take a Path
    assert z["X"].chunks == (10, 10)
    from_zarr = ad.read_zarr(zarr_pth)
    assert_equal(from_zarr, adata)


################################
# Round-tripping scanpy datasets
################################


def _do_roundtrip(
    adata: ad.AnnData, pth: Path, diskfmt: Literal["h5ad", "zarr"]
) -> ad.AnnData:
    getattr(adata, f"write_{diskfmt}")(pth)
    return getattr(ad, f"read_{diskfmt}")(pth)


@pytest.fixture
def roundtrip(diskfmt):
    return partial(_do_roundtrip, diskfmt=diskfmt)


def test_write_string_types(tmp_path, diskfmt, roundtrip):
    # https://github.com/scverse/anndata/issues/456
    adata_pth = tmp_path / f"adata.{diskfmt}"

    adata = ad.AnnData(
        obs=pd.DataFrame(
            np.ones((3, 2)),
            columns=["a", np.str_("b")],
            index=["a", "b", "c"],
        ),
    )

    from_disk = roundtrip(adata, adata_pth)

    assert_equal(adata, from_disk)


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_scanpy_pbmc68k(tmp_path, diskfmt, roundtrip, diskfmt2):
    roundtrip2 = partial(_do_roundtrip, diskfmt=diskfmt2)

    filepth1 = tmp_path / f"test1.{diskfmt}"
    filepth2 = tmp_path / f"test2.{diskfmt2}"

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r"Importing read_.* from `anndata` is deprecated"
        )
        import scanpy as sc

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ad.OldFormatWarning)
        pbmc = sc.datasets.pbmc68k_reduced()
        # zarr v3 can't write recarray
        # https://github.com/zarr-developers/zarr-python/issues/2134
        if ad.settings.zarr_write_format == 3:
            del pbmc.uns["rank_genes_groups"]["names"]
            del pbmc.uns["rank_genes_groups"]["scores"]

    from_disk1 = roundtrip(pbmc, filepth1)  # Do we read okay
    from_disk2 = roundtrip2(from_disk1, filepth2)  # Can we round trip

    assert_equal(pbmc, from_disk1)  # Not expected to be exact due to `nan`s
    assert_equal(pbmc, from_disk2)


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_scanpy_krumsiek11(tmp_path, diskfmt, roundtrip):
    filepth = tmp_path / f"test.{diskfmt}"
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r"Importing read_.* from `anndata` is deprecated"
        )
        import scanpy as sc

    # TODO: this should be fixed in scanpy instead
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        orig = sc.datasets.krumsiek11()
    del orig.uns["highlights"]  # Can’t write int keys
    # Can’t write "string" dtype: https://github.com/scverse/anndata/issues/679
    orig.obs["cell_type"] = orig.obs["cell_type"].astype(str)
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        curr = roundtrip(orig, filepth)

    assert_equal(orig, curr, exact=True)


# Checking if we can read legacy zarr files
# TODO: Check how I should add this file to the repo
@pytest.mark.filterwarnings("ignore::anndata.OldFormatWarning")
@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
@pytest.mark.skipif(
    not Path(HERE / "data/pbmc68k_reduced_legacy.zarr.zip").is_file(),
    reason="File not present.",
)
def test_backwards_compat_zarr():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r"Importing read_.* from `anndata` is deprecated"
        )
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


def test_adata_in_uns(tmp_path, diskfmt, roundtrip):
    pth = tmp_path / f"adatas_in_uns.{diskfmt}"

    orig = gen_adata((4, 5), **GEN_ADATA_NO_XARRAY_ARGS)
    orig.uns["adatas"] = {
        "a": gen_adata((1, 2), **GEN_ADATA_NO_XARRAY_ARGS),
        "b": gen_adata((12, 8), **GEN_ADATA_NO_XARRAY_ARGS),
    }
    another_one = gen_adata((2, 5), **GEN_ADATA_NO_XARRAY_ARGS)
    another_one.raw = gen_adata((2, 7), **GEN_ADATA_NO_XARRAY_ARGS)
    orig.uns["adatas"]["b"].uns["another_one"] = another_one

    curr = roundtrip(orig, pth)

    assert_equal(orig, curr)


@pytest.mark.parametrize(
    "uns_val",
    [
        pytest.param(dict(base=None), id="dict_val"),
        pytest.param(
            pd.DataFrame(dict(col_0=["string", None])).convert_dtypes(), id="df"
        ),
    ],
)
def test_none_dict_value_in_uns(diskfmt, tmp_path, roundtrip, uns_val):
    pth = tmp_path / f"adata_dtype.{diskfmt}"

    orig = ad.AnnData(np.ones((3, 4)), uns=dict(val=uns_val))
    with ad.settings.override(allow_write_nullable_strings=True):
        curr = roundtrip(orig, pth)

    if isinstance(orig.uns["val"], pd.DataFrame):
        pd.testing.assert_frame_equal(curr.uns["val"], orig.uns["val"])
    else:
        assert curr.uns["val"] == orig.uns["val"]


def test_io_dtype(tmp_path, diskfmt, dtype, roundtrip):
    pth = tmp_path / f"adata_dtype.{diskfmt}"

    orig = ad.AnnData(np.ones((5, 8), dtype=dtype))
    curr = roundtrip(orig, pth)

    assert curr.X.dtype == dtype


def test_h5py_attr_limit(tmp_path):
    N = 10_000
    a = ad.AnnData(np.ones((5, 10)))
    a.obsm["df"] = pd.DataFrame(
        np.ones((5, N)), index=a.obs_names, columns=[str(i) for i in range(N)]
    )
    a.write(tmp_path / "tmp.h5ad")


@pytest.mark.parametrize(
    "elem_key", ["obs", "var", "obsm", "varm", "layers", "obsp", "varp", "uns"]
)
def test_forward_slash_key(elem_key, tmp_path):
    a = ad.AnnData(np.ones((10, 10)))
    getattr(a, elem_key)["bad/key"] = np.ones(
        (10,) if elem_key in ["obs", "var"] else (10, 10)
    )
    with pytest.raises(ValueError, match="Forward slashes"):
        a.write_h5ad(tmp_path / "does_not_matter_the_path.h5ad")


@pytest.mark.skipif(
    find_spec("xarray"), reason="Xarray is installed so `read_lazy` will not error"
)
def test_read_lazy_import_error():
    with pytest.raises(ImportError, match="xarray"):
        ad.experimental.read_lazy("test.zarr")
