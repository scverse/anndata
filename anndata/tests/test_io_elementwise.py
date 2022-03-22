"""
Tests that each element in an anndata is written correctly
"""
from __future__ import annotations

import re

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import zarr

from anndata._io.utils import AnnDataReadError, AnnDataWriteError
from anndata.compat import _read_attr, H5Group, ZarrGroup
from anndata._io.specs import write_elem, read_elem
from anndata._io.specs.registry import NoSuchIO
from anndata.tests.helpers import assert_equal, gen_adata


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.fixture(scope="function", params=["h5", "zarr"])
def store(request, tmp_path) -> H5Group | ZarrGroup:
    if request.param == "h5":
        file = h5py.File(tmp_path / "test.h5", "w")
        store = file["/"]
    elif request.param == "zarr":
        store = zarr.open(tmp_path / "test.zarr", "w")
    else:
        assert False

    try:
        yield store
    finally:
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


@pytest.mark.parametrize(
    ["mess_with", "modifiers", "pattern"],
    [
        # (lambda s: ???, r"Unknown source type"),
        (
            lambda s: s["obs"].attrs.__setitem__("encoding-type", "floob"),
            (),
            r"Unknown encoding type “floob”",
        ),
        (
            lambda s: s["obs"].attrs.__setitem__("encoding-version", "10000.0"),
            (),
            r"Unknown encoding version 10000\.0",
        ),
        (lambda s: None, {"x"}, r"Unknown modifier set.+\{'x'\}"),
    ],
)
def test_read_io_error(store, mess_with, modifiers, pattern):
    adata = gen_adata((3, 2))

    write_elem(store, "/", adata)
    mess_with(store)
    full_pattern = rf"No such read function registered: {pattern}.*try updating"
    with pytest.raises(
        AnnDataReadError, match=r"while reading key '/(obs)?'"
    ) as exc_info:
        read_elem(store, modifiers=frozenset(modifiers))
    msg = str(exc_info.value.__cause__)
    assert re.match(full_pattern, msg)


@pytest.mark.parametrize(
    ["obj", "modifiers", "pattern"],
    [
        (b"x", (), r"Source type <class 'bytes'> not found for dest.*Group"),
        # (???, (), r"Dest.*not found for.*source.*asdfsdf"),
        (np.ndarray([1]), {"x"}, r"Unknown modifier set.+\{'x'\}"),
    ],
)
def test_write_io_error(store, obj, modifiers, pattern):
    full_pattern = re.compile(rf"No such write function registered: {pattern}")
    with pytest.raises(AnnDataWriteError, match=r"while writing key '/el'") as exc_info:
        write_elem(store, "/el", obj, modifiers=frozenset(modifiers))
    msg = str(exc_info.value.__cause__)
    assert re.match(full_pattern, msg)
