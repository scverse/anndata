from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata import AnnData, concat
from anndata._core.merge import _resolve_axis
from anndata.experimental.merge import as_group, concat_on_disk
from anndata.io import read_elem, write_elem
from anndata.tests.helpers import assert_equal, gen_adata
from anndata.utils import asarray

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal


GEN_ADATA_OOC_CONCAT_ARGS = dict(
    obsm_types=(
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
    ),
    varm_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame),
    layers_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame),
)


@pytest.fixture(params=[0, 1])
def axis(request) -> Literal[0, 1]:
    return request.param


@pytest.fixture(params=["array", "sparse", "sparse_array"])
def array_type(request) -> Literal["array", "sparse", "sparse_array"]:
    return request.param


@pytest.fixture(params=["inner", "outer"])
def join_type(request) -> Literal["inner", "outer"]:
    return request.param


@pytest.fixture(params=["zarr", "h5ad"])
def file_format(request) -> Literal["zarr", "h5ad"]:
    return request.param


# 1000 is enough to guarantee that the feature is being used
@pytest.fixture(params=[1_000, 100_000_000])
def max_loaded_elems(request) -> int:
    return request.param


def _adatas_to_paths(adatas, tmp_path, file_format):
    """
    Gets list of adatas, writes them and returns their paths as zarr
    """
    paths = None

    if isinstance(adatas, Mapping):
        paths = {}
        for k, v in adatas.items():
            p = tmp_path / (f"{k}." + file_format)
            write_elem(as_group(p, mode="a"), "", v)
            paths[k] = p
    else:
        paths = []
        for i, a in enumerate(adatas):
            p = tmp_path / (f"{i}." + file_format)
            write_elem(as_group(p, mode="a"), "", a)
            paths += [p]
    return paths


def assert_eq_concat_on_disk(
    adatas,
    tmp_path: Path,
    file_format: Literal["zarr", "h5ad"],
    max_loaded_elems: int | None = None,
    *args,
    **kwargs,
):
    # create one from the concat function
    res1 = concat(adatas, *args, **kwargs)
    # create one from the on disk concat function
    paths = _adatas_to_paths(adatas, tmp_path, file_format)
    out_name = tmp_path / f"out.{file_format}"
    if max_loaded_elems is not None:
        kwargs["max_loaded_elems"] = max_loaded_elems
    concat_on_disk(paths, out_name, *args, **kwargs)
    res2 = read_elem(as_group(out_name, mode="r"))
    assert_equal(res1, res2, exact=False)


def get_array_type(array_type, axis):
    if array_type == "sparse":
        return sparse.csr_matrix if axis == 0 else sparse.csc_matrix
    if array_type == "sparse_array":
        return sparse.csr_array if axis == 0 else sparse.csc_array
    if array_type == "array":
        return asarray
    msg = f"array_type {array_type} not implemented"
    raise NotImplementedError(msg)


@pytest.mark.parametrize("reindex", [True, False], ids=["reindex", "no_reindex"])
def test_anndatas(
    *,
    axis: Literal[0, 1],
    array_type: Literal["array", "sparse", "sparse_array"],
    join_type: Literal["inner", "outer"],
    tmp_path: Path,
    max_loaded_elems: int,
    file_format: Literal["zarr", "h5ad"],
    reindex: bool,
):
    _, off_axis_name = _resolve_axis(1 - axis)
    random_axes = {0, 1} if reindex else {axis}
    sparse_fmt = "csr" if axis == 0 else "csc"
    kw = (
        GEN_ADATA_OOC_CONCAT_ARGS
        if not reindex
        else dict(
            obsm_types=(get_array_type("sparse", 1 - axis), np.ndarray, pd.DataFrame),
            varm_types=(get_array_type("sparse", 1 - axis), np.ndarray, pd.DataFrame),
            layers_types=(get_array_type("sparse", axis), np.ndarray, pd.DataFrame),
        )
    )

    adatas = []
    for i in range(5):
        M, N = (np.random.randint(1, 100) if a in random_axes else 50 for a in (0, 1))
        a = gen_adata(
            (M, N), X_type=get_array_type(array_type, axis), sparse_fmt=sparse_fmt, **kw
        )
        # ensure some names overlap, others do not, for the off-axis so that inner/outer is properly tested
        off_names = getattr(a, f"{off_axis_name}_names").array
        off_names[1::2] = f"{i}-" + off_names[1::2]
        setattr(a, f"{off_axis_name}_names", off_names)
        adatas.append(a)

    assert_eq_concat_on_disk(
        adatas,
        tmp_path,
        file_format,
        max_loaded_elems,
        axis=axis,
        join=join_type,
    )


def test_concat_ordered_categoricals_retained(tmp_path, file_format):
    a = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("aabcd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5)],
        ),
    )
    b = AnnData(
        X=np.ones((5, 1)),
        obs=pd.DataFrame(
            {
                "cat_ordered": pd.Categorical(list("abcdd"), ordered=True),
            },
            index=[f"cell{i:02}" for i in range(5, 10)],
        ),
    )

    adatas = [a, b]
    assert_eq_concat_on_disk(adatas, tmp_path, file_format)


@pytest.fixture
def xxxm_adatas():
    def gen_index(n):
        return [f"cell{i}" for i in range(n)]

    return [
        AnnData(
            X=sparse.csr_matrix((3, 5)),
            obs=pd.DataFrame(index=gen_index(3)),
            obsm={
                "dense": np.arange(6).reshape(3, 2),
                "sparse": sparse.csr_matrix(np.arange(6).reshape(3, 2)),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(3),
                        "b": list("abc"),
                        "c": pd.Categorical(list("aab")),
                    },
                    index=gen_index(3),
                ),
            },
        ),
        AnnData(
            X=sparse.csr_matrix((4, 10)),
            obs=pd.DataFrame(index=gen_index(4)),
            obsm=dict(
                dense=np.arange(12).reshape(4, 3),
                df=pd.DataFrame(dict(a=np.arange(3, 7)), index=gen_index(4)),
            ),
        ),
        AnnData(
            X=sparse.csr_matrix((2, 100)),
            obs=pd.DataFrame(index=gen_index(2)),
            obsm={
                "sparse": np.arange(8).reshape(2, 4),
                "dense": np.arange(4, 8).reshape(2, 2),
                "df": pd.DataFrame(
                    {
                        "a": np.arange(7, 9),
                        "b": list("cd"),
                        "c": pd.Categorical(list("ab")),
                    },
                    index=gen_index(2),
                ),
            },
        ),
    ]


def test_concatenate_xxxm(xxxm_adatas, tmp_path, file_format, join_type):
    if join_type == "outer":
        for i in range(len(xxxm_adatas)):
            xxxm_adatas[i] = xxxm_adatas[i].T
            xxxm_adatas[i].X = sparse.csr_matrix(xxxm_adatas[i].X)
    assert_eq_concat_on_disk(xxxm_adatas, tmp_path, file_format, join=join_type)


def test_output_dir_exists(tmp_path):
    in_pth = tmp_path / "in.h5ad"
    out_pth = tmp_path / "does_not_exist" / "out.h5ad"

    AnnData(X=np.ones((5, 1))).write_h5ad(in_pth)

    with pytest.raises(FileNotFoundError, match=f"{out_pth}"):
        concat_on_disk([in_pth], out_pth)


def test_failure_w_no_args(tmp_path):
    with pytest.raises(ValueError, match=r"No objects to concatenate"):
        concat_on_disk([], tmp_path / "out.h5ad")
