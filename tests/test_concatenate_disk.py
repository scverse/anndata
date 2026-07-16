from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import awkward as ak
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData, concat
from anndata._core import merge
from anndata._core.merge import _resolve_axis
from anndata.experimental.merge import as_group, concat_on_disk
from anndata.io import read_elem, write_elem
from anndata.tests.helpers import assert_equal, check_all_sharded, gen_adata
from anndata.utils import asarray

if TYPE_CHECKING:
    from collections.abc import Collection
    from pathlib import Path
    from typing import Literal

    from anndata._types import Join_T


GEN_ADATA_OOC_CONCAT_ARGS = dict(
    obsm_types=(
        sparse.csr_matrix,
        np.ndarray,
        pd.DataFrame,
    ),
    varm_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame),
    layers_types=(sparse.csr_matrix, np.ndarray, pd.DataFrame),
)


@pytest.fixture(params=list(merge.MERGE_STRATEGIES.keys()))
def merge_strategy(request):
    return request.param


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


# 5 is enough to guarantee that the feature is being used since the
# `test_anndatas` generates a minimum of 5 on at least one of the axes.
# Thus there will be at least 5 elems.
@pytest.fixture(params=[5, 1_000_000])
def max_loaded_elems(request) -> int:
    return request.param


def _adatas_to_paths(
    adatas: Mapping[str, AnnData] | Collection[AnnData],
    tmp_path: Path,
    file_format: str,
) -> dict[str, Path] | list[Path]:
    """Gets list of adatas, writes them and returns their paths as zarr."""
    paths = {}
    for k, v in adatas.items() if isinstance(adatas, Mapping) else enumerate(adatas):
        p = tmp_path / f"{k}.{file_format}"
        with as_group(p, mode="a") as f:
            write_elem(f, "/", v)
        paths[k] = p
    return paths if isinstance(adatas, Mapping) else list(paths.values())


def assert_eq_concat_on_disk(
    adatas: Mapping[str, AnnData] | Collection[AnnData],
    tmp_path: Path,
    file_format: Literal["zarr", "h5ad"],
    max_loaded_elems: int | None = None,
    *args,
    merge_strategy: merge.StrategiesLiteral | None = None,
    **kwargs,
):
    # create one from the concat function
    res1 = concat(adatas, *args, merge=merge_strategy, **kwargs)
    # create one from the on disk concat function
    paths = _adatas_to_paths(adatas, tmp_path, file_format)
    out_name = tmp_path / f"out.{file_format}"
    if max_loaded_elems is not None:
        kwargs["max_loaded_elems"] = max_loaded_elems
    concat_on_disk(paths, out_name, *args, merge=merge_strategy, **kwargs)
    with as_group(out_name, mode="r") as rg:
        res2 = read_elem(rg)
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
@pytest.mark.filterwarnings("ignore:Misaligned chunks detected")
def test_anndatas(
    *,
    axis: Literal[0, 1],
    array_type: Literal["array", "sparse", "sparse_array"],
    join_type: Join_T,
    tmp_path: Path,
    max_loaded_elems: int,
    file_format: Literal["zarr", "h5ad"],
    reindex: bool,
    merge_strategy: merge.StrategiesLiteral,
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
    for i in range(3):
        M, N = (np.random.randint(5, 10) if a in random_axes else 50 for a in (0, 1))
        a = gen_adata(
            (M, N),
            X_type=get_array_type(array_type, axis),
            sparse_fmt=sparse_fmt,
            obs_dtypes=[pd.CategoricalDtype(ordered=False)],
            var_dtypes=[pd.CategoricalDtype(ordered=False)],
            **kw,
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
        merge_strategy=merge_strategy,
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
                "sparse": sparse.csr_matrix(np.arange(8).reshape(2, 4)),
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


def test_concat_on_disk_outer_mapping_missing_keys(tmp_path, file_format):
    """`.obsm` keys absent from some objects are retained (and filled) under an outer join.

    Regression test for https://github.com/scverse/anndata/issues/2394. ``dense`` is present
    in every object but with differing widths (outer reindexing), while ``sparse`` and the
    dataframe ``df`` are missing from the middle object (union-of-keys + fill) -- the
    concatenation-axis mapping path that :func:`test_concatenate_xxxm` sidesteps by
    transposing. Values are numeric so the filled result is writable across supported
    ``h5py``/``zarr`` versions.
    """

    def make(i, n, dense_w, sparse_w, *, with_df):
        idx = [f"{i}-cell{j}" for j in range(n)]
        obsm = {"dense": np.arange(n * dense_w).reshape(n, dense_w).astype(float)}
        if sparse_w is not None:
            obsm["sparse"] = sparse.csr_matrix(
                np.arange(n * sparse_w).reshape(n, sparse_w).astype(float)
            )
        if with_df:
            obsm["df"] = pd.DataFrame({"s": np.arange(n, dtype=float)}, index=idx)
        return AnnData(
            sparse.csr_matrix((n, 5)),
            obs=pd.DataFrame(index=idx),
            obsm=obsm,
        )

    adatas = [
        make(0, 3, 2, 2, with_df=True),
        make(1, 4, 3, None, with_df=False),
        make(2, 2, 2, 4, with_df=True),
    ]
    assert_eq_concat_on_disk(adatas, tmp_path, file_format, join="outer")


def test_concat_on_disk_outer_layers_missing_keys(tmp_path, file_format):
    """`.layers` keys absent from some objects are retained (and filled) under an outer join.

    Regression test for https://github.com/scverse/anndata/issues/2394; ``var`` differs
    between objects so the mapping is reindexed along the alternate axis as well.
    """
    a = AnnData(
        np.ones((3, 3)),
        obs=pd.DataFrame(index=[f"a{i}" for i in range(3)]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
        layers={"counts": np.arange(9).reshape(3, 3).astype(float)},
    )
    b = AnnData(
        np.ones((2, 3)),
        obs=pd.DataFrame(index=[f"b{i}" for i in range(2)]),
        var=pd.DataFrame(index=["g2", "g3", "g4"]),
    )
    assert_eq_concat_on_disk([a, b], tmp_path, file_format, join="outer")


def test_concatenate_zarr_stays_sharded_v3(xxxm_adatas, tmp_path):
    import zarr

    assert_eq_concat_on_disk(xxxm_adatas, tmp_path, file_format="zarr")
    g = zarr.open(tmp_path)
    assert g.metadata.zarr_format == 3

    check_all_sharded(g)


def test_singleton(xxxm_adatas, tmp_path, file_format):
    assert_eq_concat_on_disk(xxxm_adatas[:1], tmp_path, file_format=file_format)


def test_output_dir_exists(tmp_path):
    in_pth = tmp_path / "in.h5ad"
    out_pth = tmp_path / "does_not_exist" / "out.h5ad"

    AnnData(X=np.ones((5, 1))).write_h5ad(in_pth)

    with pytest.raises(FileNotFoundError, match=str(out_pth)):
        concat_on_disk([in_pth], out_pth)


def test_no_open_h5_file_handles_after_error(tmp_path):
    in_pth = tmp_path / "in.h5ad"
    in_pth2 = tmp_path / "in2.h5ad"
    out_pth = tmp_path / "out.h5ad"

    adata = AnnData(
        X=np.ones((2, 1)),
        obsm={
            "awk": ak.Array([
                [{"a": 1, "b": "foo"}],
                [{"a": 2, "b": "bar"}],
            ])
        },
    )
    adata.write_h5ad(in_pth)
    adata.write_h5ad(in_pth2)

    # Intentionally write an unsupported array type, which could leave dangling file handles:
    # https://github.com/scverse/anndata/issues/2198
    try:
        concat_on_disk([in_pth, in_pth2], out_pth)
    except NotImplementedError:
        for path in [in_pth, in_pth2, out_pth]:
            # should not error out because there are no file handles open
            f = h5py.File(path, mode="w")
            f.close()


def test_write_using_groups(tmp_path, file_format):
    in_pth = tmp_path / f"in.{file_format}"
    in_pth2 = tmp_path / f"in2.{file_format}"
    out_pth = tmp_path / f"out.{file_format}"

    adata = AnnData(X=np.ones((2, 1)))
    getattr(adata, f"write_{file_format}")(in_pth)
    getattr(adata, f"write_{file_format}")(in_pth2)

    with (
        as_group(in_pth, mode="r") as f1,
        as_group(in_pth2, mode="r") as f2,
        as_group(out_pth, mode="w") as fout,
    ):
        concat_on_disk([f1, f2], fout)
    adata_out = getattr(ad, f"read_{file_format}")(out_pth)
    assert_equal(adata_out, concat([adata, adata]))


def test_failure_w_no_args(tmp_path):
    with pytest.raises(ValueError, match=r"No objects to concatenate"):
        concat_on_disk([], tmp_path / "out.h5ad")
