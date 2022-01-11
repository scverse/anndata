"""
Tests that each element in an anndata is written correctly
"""
from tempfile import TemporaryDirectory
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import zarr


import anndata as ad
from anndata.compat import _read_attr
from anndata._io.specs import write_elem, read_elem
from anndata.tests.helpers import assert_equal, gen_adata


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.fixture(scope="function", params=["h5", "zarr"])
def store(request):
    with TemporaryDirectory() as tmpdir:

        if request.param == "h5":
            file = h5py.File(Path(tmpdir) / "test.h5", "w")
            store = file["/"]
        elif request.param == "zarr":
            store = zarr.open(Path(tmpdir) / "test.zarr")

        yield store

        if request.param == "h5":
            file.close()


@pytest.mark.parametrize(
    "value,encoding_type",
    [
        ("hello world", "string"),
        (np.str_("hello world"), "string"),
        (np.array([1, 2, 3]), "array"),
        (np.array(["hello", "world"], dtype=object), "string-array"),
        (1, "numeric-scalar"),
        (True, "numeric-scalar"),
        (1.0, "numeric-scalar"),
        ({"a": 1}, "dict"),
        (gen_adata((3, 2)), "anndata"),
        (sparse.random(5, 3, format="csr", density=0.5), "csr_matrix"),
        (sparse.random(5, 3, format="csc", density=0.5), "csc_matrix"),
        (pd.DataFrame({"a": [1, 2, 3]}), "dataframe"),
        (pd.Categorical(list("aabccedd")), "categorical"),
        (pd.Categorical(list("aabccedd"), ordered=True), "categorical"),
        (pd.Categorical([1, 2, 1, 3], ordered=True), "categorical"),
        (
            pd.arrays.IntegerArray(
                np.ones(5, dtype=int), mask=np.array([True, False, True, False, True])
            ),
            "nullable-integer",
        ),
        (pd.array([1, 2, 3]), "nullable-integer"),
        (
            pd.arrays.BooleanArray(
                np.random.randint(0, 2, size=5, dtype=bool),
                mask=np.random.randint(0, 2, size=5, dtype=bool),
            ),
            "nullable-boolean",
        ),
        (pd.array([True, False, True, True]), "nullable-boolean"),
        # (bytes, b"some bytes", "bytes"), # Does not work for zarr
        # TODO consider how specific encodings should be. Should we be fully describing the written type?
        # Currently the info we add is: "what you wouldn't be able to figure out yourself"
        # but that's not really a solid rule.
        # (bool, True, "bool"),
        # (bool, np.bool_(False), "bool"),
    ],
)
def test_io_spec(store, value, encoding_type):
    key = f"key_for_{encoding_type}"
    write_elem(store, key, value, dataset_kwargs={})

    assert encoding_type == _read_attr(store[key].attrs, "encoding-type")

    from_disk = read_elem(store[key])
    assert_equal(value, from_disk)


def test_io_spec_raw(store):
    adata = gen_adata((3, 2))
    adata.raw = adata

    write_elem(store, "adata", adata)

    assert "raw" == _read_attr(store["adata/raw"].attrs, "encoding-type")

    from_disk = read_elem(store["adata"])
    assert_equal(from_disk.raw, adata.raw)


def test_write_to_root(store):
    adata = gen_adata((3, 2))

    write_elem(store, "/", adata)
    from_disk = read_elem(store)

    assert "anndata" == _read_attr(store.attrs, "encoding-type")
    assert_equal(from_disk, adata)
