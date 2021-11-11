from functools import singledispatch
from collections.abc import Mapping

import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy import sparse

from anndata import AnnData, Raw
from anndata._core.views import ArrayView
from anndata._core.sparse_dataset import SparseDataset
from anndata._core.aligned_mapping import AlignedMapping
from anndata.testing.exceptions import report_name, format_msg
from anndata.utils import asarray


@report_name
def _assert_equal(a, b):
    """Allows reporting elem name for simple assertion."""
    assert a == b


@singledispatch
def assert_equal(a, b, *, exact=False, elem_name=None):
    _assert_equal(a, b, _elem_name=elem_name)


@assert_equal.register(np.ndarray)
def assert_equal_ndarray(a, b, *, exact=False, elem_name=None):
    b = asarray(b)
    if not exact and is_numeric_dtype(a) and is_numeric_dtype(b):
        assert a.shape == b.shape, format_msg(elem_name)
        assert np.allclose(a, b, equal_nan=True), format_msg(elem_name)
    elif (  # Structured dtype
        not exact
        and hasattr(a, "dtype")
        and hasattr(b, "dtype")
        and len(a.dtype) > 1
        and len(b.dtype) > 0
    ):
        assert_equal(pd.DataFrame(a), pd.DataFrame(b), exact=exact, elem_name=elem_name)
    else:
        assert np.all(a == b), format_msg(elem_name)


@assert_equal.register(ArrayView)
def assert_equal_arrayview(a, b, *, exact=False, elem_name=None):
    assert_equal(asarray(a), asarray(b), exact=exact, elem_name=elem_name)


@assert_equal.register(SparseDataset)
@assert_equal.register(sparse.spmatrix)
def assert_equal_sparse(a, b, *, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact=exact, elem_name=elem_name)


@assert_equal.register(h5py.Dataset)
def assert_equal_h5py_dataset(a, b, *, exact=False, elem_name=None):
    a = asarray(a)
    assert_equal(b, a, exact=exact, elem_name=elem_name)


@assert_equal.register(pd.DataFrame)
def are_equal_dataframe(a, b, *, exact=False, elem_name=None):
    if not isinstance(b, pd.DataFrame):
        assert_equal(b, a, exact=exact, elem_name=elem_name)  # , a.values maybe?

    report_name(pd.testing.assert_frame_equal)(
        a,
        b,
        check_index_type=exact,
        check_exact=exact,
        _elem_name=elem_name,
        check_frame_type=False,
    )


@assert_equal.register(Mapping)
def assert_equal_mapping(a, b, *, exact=False, elem_name=None):
    assert set(a.keys()) == set(b.keys()), format_msg(elem_name)
    for k in a.keys():
        if elem_name is None:
            elem_name = ""
        assert_equal(a[k], b[k], exact=exact, elem_name=f"{elem_name}/{k}")


@assert_equal.register(AlignedMapping)
def assert_equal_aligned_mapping(a, b, *, exact=False, elem_name=None):
    a_indices = (a.parent.obs_names, a.parent.var_names)
    b_indices = (b.parent.obs_names, b.parent.var_names)
    for axis_idx in a.axes:
        assert_equal(
            a_indices[axis_idx], b_indices[axis_idx], exact=exact, elem_name=axis_idx
        )
    assert a.attrname == b.attrname, format_msg(elem_name)
    assert_equal_mapping(a, b, exact=exact, elem_name=elem_name)


@assert_equal.register(pd.Index)
def assert_equal_index(a, b, *, exact=False, elem_name=None):
    if not exact:
        report_name(pd.testing.assert_index_equal)(
            a, b, check_names=False, check_categorical=False, _elem_name=elem_name
        )
    else:
        report_name(pd.testing.assert_index_equal)(a, b, _elem_name=elem_name)


@assert_equal.register(Raw)
def assert_equal_raw(a, b, *, exact=False, elem_name=None):
    def assert_is_not_none(x):  # can't put an assert in a lambda
        assert x is not None

    report_name(assert_is_not_none)(b, _elem_name=elem_name)
    for attr in ["X", "var", "varm", "obs_names"]:
        assert_equal(
            getattr(a, attr),
            getattr(b, attr),
            exact=exact,
            elem_name=f"{elem_name}/{attr}",
        )


@assert_equal.register(AnnData)
def assert_adata_equal(a: AnnData, b: AnnData, *, exact: bool = False):
    """\
    Check whether two AnnData objects are equivalent,
    raising an AssertionError if they aren’t.

    Params
    ------
    a
    b
    exact
        Whether comparisons should be exact or not. This has a somewhat flexible
        meaning and should probably get refined in the future.
    """
    # There may be issues comparing views, since np.allclose
    # can modify ArrayViews if they contain `nan`s
    assert_equal(a.obs_names, b.obs_names, exact=exact, elem_name="obs_names")
    assert_equal(a.var_names, b.var_names, exact=exact, elem_name="var_names")
    if not exact:
        # Reorder all elements if neccesary
        idx = [slice(None), slice(None)]
        # Since it’s a pain to compare a list of pandas objects
        change_flag = False
        if not np.all(a.obs_names == b.obs_names):
            idx[0] = a.obs_names
            change_flag = True
        if not np.all(a.var_names == b.var_names):
            idx[1] = a.var_names
            change_flag = True
        if change_flag:
            b = b[tuple(idx)].copy()
    for attr in [
        "X",
        "obs",
        "var",
        "obsm",
        "varm",
        "layers",
        "uns",
        "obsp",
        "varp",
        "raw",
    ]:
        assert_equal(
            getattr(a, attr),
            getattr(b, attr),
            exact=exact,
            elem_name=attr,
        )
