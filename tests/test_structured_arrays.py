from __future__ import annotations

import warnings
from itertools import combinations, product
from typing import TYPE_CHECKING

import numpy as np

import anndata as ad
from anndata import AnnData
from anndata.compat import is_zarr_v2
from anndata.tests.helpers import gen_vstr_recarray

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal


def assert_str_contents_equal(
    a: Iterable[Iterable[object]], b: Iterable[Iterable[object]]
) -> None:
    l_a = [
        [str(el) if not isinstance(el, bytes) else el.decode("utf-8") for el in row]
        for row in a
    ]
    l_b = [
        [str(el) if not isinstance(el, bytes) else el.decode("utf-8") for el in row]
        for row in b
    ]
    assert l_a == l_b


def test_io(
    tmp_path, diskfmt: Literal["zarr", "h5ad"], diskfmt2: Literal["zarr", "h5ad"]
) -> None:
    if not is_zarr_v2():
        from zarr.core.dtype.common import UnstableSpecificationWarning

        warnings.filterwarnings(  # raised by “S10” dtype in the recarray below
            "default", r".*NullTerminatedBytes", UnstableSpecificationWarning
        )

    read1 = lambda pth: getattr(ad, f"read_{diskfmt}")(pth)
    write1 = lambda adata, pth: getattr(adata, f"write_{diskfmt}")(pth)
    read2 = lambda pth: getattr(ad, f"read_{diskfmt2}")(pth)
    write2 = lambda adata, pth: getattr(adata, f"write_{diskfmt2}")(pth)

    filepth1 = tmp_path / f"test1.{diskfmt}"
    filepth2 = tmp_path / f"test2.{diskfmt2}"

    str_recarray = gen_vstr_recarray(3, 5)
    u_recarray = str_recarray.astype([(n, "U10") for n in str_recarray.dtype.fields])
    s_recarray = str_recarray.astype([(n, "S10") for n in str_recarray.dtype.fields])

    initial = AnnData(np.zeros((3, 3)))
    initial.uns = dict(str_rec=str_recarray, u_rec=u_recarray, s_rec=s_recarray)

    write1(initial, filepth1)
    disk_once = read1(filepth1)
    write2(disk_once, filepth2)
    disk_twice = read2(filepth2)

    adatas = [initial, disk_once, disk_twice]
    keys = [
        "str_rec",
        "u_rec",
        # "s_rec"
    ]

    for (ad1, key1), (ad2, key2) in combinations(product(adatas, keys), 2):
        assert_str_contents_equal(ad1.uns[key1], ad2.uns[key2])
