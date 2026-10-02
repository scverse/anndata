from __future__ import annotations

import narwhals as nw
import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata._core._dataframe_backend import (
    DataFrameLike,
    _from_backend,
    _ingest_as_pandas,
    _ingest_axis_frame,
    axis_index,
    column_backed_dim,
    copy_frame,
    frame_annotation_columns,
    frame_column,
    frame_equal,
    native_backend,
    relabel_axis,
    set_axis_index,
    subset_frame,
    to_backend,
    true_axis_index,
)
from anndata._core.raw import Raw
from anndata._core.xarray import Dataset2D
from anndata.compat import XDataset, XVariable
from anndata.tests.helpers import assert_equal, gen_typed_df

pytest.importorskip("xarray")


@pytest.fixture
def df():
    return gen_typed_df(10)


@pytest.fixture
def dataset2d(df):
    return Dataset2D(XDataset.from_dataframe(df))


@pytest.fixture
def lazy_dataset2d():
    """A `Dataset2D` mixing a chunked dask column with an in-memory categorical one."""
    da = pytest.importorskip("dask.array")
    n = 50
    return Dataset2D(
        XDataset(
            {
                "val": ("idx", da.from_array(np.arange(n, dtype="float64"), chunks=20)),
                "grp": ("idx", pd.array(["a", "b"] * (n // 2), dtype="category")),
            },
            coords={"idx": [f"c{i}" for i in range(n)]},
        )
    )


@pytest.fixture
def named_obs():
    return pd.DataFrame(
        {"cell_type": pd.Categorical(["A", "B", "A"]), "n_genes": [10, 20, 30]},
        index=pd.Index(["AAAC", "AAAG", "AAAT"], name="obs_names"),
    )


@pytest.fixture
def dask_obs(named_obs):
    dd = pytest.importorskip("dask.dataframe")
    return dd.from_pandas(named_obs, npartitions=2)


def make_native_frame(backend, data):
    module = pytest.importorskip(backend)
    return module.DataFrame(data) if backend == "polars" else module.table(data)


def test_narwhals_from_native_roundtrip(df, dataset2d):
    frame = nw.from_native(dataset2d)
    assert isinstance(frame, nw.DataFrame)
    assert set(frame.columns) == set(df.columns)
    assert_equal(
        frame.to_native().sort_index(axis=1), dataset2d.to_memory().sort_index(axis=1)
    )


def test_narwhals_op_parity_with_pandas(df, dataset2d):
    num = df.select_dtypes("number").columns.tolist()
    predicate = nw.col(num[0]) >= df[num[0]].median()
    from_ds = nw.from_native(dataset2d).filter(predicate).select(num).to_native()
    from_df = (
        nw.from_native(df, eager_only=True).filter(predicate).select(num).to_native()
    )
    assert_equal(from_ds.reset_index(drop=True), from_df.reset_index(drop=True))


def test_narwhals_index_preserved(dataset2d):
    idx = nw.maybe_get_index(nw.from_native(dataset2d))
    assert_equal(idx, dataset2d.to_memory().index)


def test_narwhals_categorical_preserved(named_obs):
    ds = Dataset2D(XDataset.from_dataframe(named_obs))
    out = nw.from_native(ds).to_native()
    assert isinstance(out["cell_type"].dtype, pd.CategoricalDtype)
    assert out["cell_type"].tolist() == ["A", "B", "A"]


def test_narwhals_lazy_dataset2d_is_a_lazyframe(lazy_dataset2d):
    frame = nw.from_native(lazy_dataset2d)
    assert isinstance(frame, nw.LazyFrame)
    assert frame.implementation is nw.Implementation.DASK
    assert set(frame.columns) == {"val", "grp"}

    collected = frame.collect().to_native()
    assert collected["val"].tolist() == list(range(lazy_dataset2d.shape[0]))
    assert isinstance(collected["grp"].dtype, pd.CategoricalDtype)


def test_narwhals_lazy_dataset2d_defers_the_read(lazy_dataset2d):
    """Wrapping and building a query must not touch the data behind the dask columns."""
    frame = nw.from_native(lazy_dataset2d).filter(nw.col("val") >= 45)
    assert isinstance(frame, nw.LazyFrame)
    assert frame.collect()["val"].to_list() == [45.0, 46.0, 47.0, 48.0, 49.0]


def test_dataset2d_to_dask_dataframe_partitions_follow_chunks(lazy_dataset2d):
    ddf = lazy_dataset2d.to_dask_dataframe()
    assert ddf.npartitions == 3  # 50 rows in chunks of 20
    assert isinstance(ddf.dtypes["grp"], pd.CategoricalDtype)

    computed = ddf.compute()
    assert computed.index.tolist() == [f"c{i}" for i in range(50)]
    assert computed["val"].tolist() == list(range(50))


def test_dataset2d_in_memory_is_not_lazy(dataset2d):
    assert not dataset2d.is_lazy
    assert isinstance(nw.from_native(dataset2d), nw.DataFrame)
    assert dataset2d.to_dask_dataframe().npartitions == 1


@pytest.mark.parametrize("backend", ["pandas", "pyarrow", "polars"])
def test_to_backend_from_dataset2d(dataset2d, backend):
    pytest.importorskip(backend)

    native = to_backend(dataset2d, backend)
    assert type(native).__module__.split(".")[0] == backend
    got = nw.from_native(native, eager_only=True).to_arrow()
    src = nw.from_native(dataset2d).to_arrow()
    assert set(dataset2d.columns) <= set(got.column_names)
    for name in dataset2d.columns:
        assert got.column(name).to_pylist() == src.column(name).to_pylist()


def test_to_backend_named_index_becomes_names_column(named_obs):
    ds = Dataset2D(XDataset.from_dataframe(named_obs))
    for backend in ("polars", "pyarrow"):
        pytest.importorskip(backend)
        got = nw.from_native(to_backend(ds, backend), eager_only=True).to_arrow()
        assert got.column("obs_names").to_pylist() == ["AAAC", "AAAG", "AAAT"]
        assert got.column("n_genes").to_pylist() == [10, 20, 30]
    pandas_frame = to_backend(ds, "pandas")
    assert pandas_frame.index.name == "obs_names"
    assert "obs_names" not in pandas_frame.columns


def test_obs_var_as_in_memory(named_obs):
    pytest.importorskip("polars")
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=named_obs)

    obs_pl = adata.obs_as("polars")
    assert type(obs_pl).__module__.split(".")[0] == "polars"
    assert nw.from_native(obs_pl, eager_only=True)["cell_type"].to_list() == [
        "A",
        "B",
        "A",
    ]

    obs_pd = adata.obs_as("pandas")
    assert isinstance(obs_pd, pd.DataFrame)
    assert obs_pd.index.name == "obs_names"

    var_pl = adata.var_as("polars")
    assert type(var_pl).__module__.split(".")[0] == "polars"


def test_obs_var_as_read_lazy(tmp_path, diskfmt, named_obs):
    pytest.importorskip("polars")
    from anndata.acc import A

    var = pd.DataFrame(
        {"gene_type": pd.Categorical(["tf", "rp", "tf"])},
        index=pd.Index(["G1", "G2", "G3"], name="var_names"),
    )
    p = tmp_path / f"a.{diskfmt}"
    getattr(
        ad.AnnData(X=np.zeros((3, 3), "f4"), obs=named_obs, var=var),
        f"write_{diskfmt}",
    )(p)
    lazy = ad.experimental.read_lazy(p)
    assert isinstance(lazy.obs, Dataset2D)
    assert isinstance(lazy.var, Dataset2D)
    obs_pl = lazy.obs_as("polars")
    var_pl = lazy.var_as("polars")
    assert type(obs_pl).__module__.split(".")[0] == "polars"
    assert "obs_names" in obs_pl.columns
    assert var_pl["var_names"].to_list() == ["G1", "G2", "G3"]
    assert var_pl["gene_type"].to_list() == ["tf", "rp", "tf"]
    assert lazy[A.obs.index].tolist() == ["AAAC", "AAAG", "AAAT"]
    assert lazy[A.var["gene_type"]].values.tolist() == ["tf", "rp", "tf"]


def test_to_backend_unnamed_index_named(named_obs):
    pytest.importorskip("polars")
    obs = pd.DataFrame({"g": [1, 2, 3]}, index=["c0", "c1", "c2"])
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=obs)
    assert adata.obs.index.name is None
    cols = adata.obs_as("polars").columns
    assert "obs_names" in cols
    assert not any(c.startswith("__index_level_") for c in cols)
    assert "var_names" in adata.var_as("polars").columns


def test_to_backend_renames_arbitrarily_named_index():
    pl = pytest.importorskip("polars")
    obs = pd.DataFrame({"kind": ["a", "b"]}, index=pd.Index(["c0", "c1"], name="cells"))
    out = to_backend(obs, "polars", dim="obs")
    assert isinstance(out, pl.DataFrame)
    assert out.columns == ["kind", "obs_names"]
    assert out["obs_names"].to_list() == ["c0", "c1"]


def test_to_backend_rejects_reserved_names_column_collision():
    pytest.importorskip("polars")
    obs = pd.DataFrame(
        {"obs_names": ["annotation-a", "annotation-b"]},
        index=pd.Index(["c0", "c1"], name="cells"),
    )
    with pytest.raises(ValueError, match="reserved column 'obs_names'"):
        to_backend(obs, "polars", dim="obs")


def test_obs_as_narwhals_wraps_storage_as_it_is(named_obs):
    """Every storage form answers the same narwhals query; only lazy storage stays lazy."""
    pytest.importorskip("polars")
    dd = pytest.importorskip("dask.dataframe")
    X = np.zeros((3, 2), "f4")
    stored = {
        "pandas": named_obs,
        "polars": make_native_frame(
            "polars",
            {"obs_names": ["AAAC", "AAAG", "AAAT"], "n_genes": [10, 20, 30]},
        ),
        "dask": dd.from_pandas(named_obs, npartitions=2),
    }
    for name, obs in stored.items():
        frame = ad.AnnData(X, obs=obs).obs_as("narwhals")
        assert isinstance(frame, nw.LazyFrame if name == "dask" else nw.DataFrame)
        query = frame.filter(nw.col("n_genes") >= 20).select("n_genes")
        eager = query.collect() if isinstance(query, nw.LazyFrame) else query
        assert list(eager.to_native()["n_genes"]) == [20, 30]


def test_obs_as_narwhals_keeps_pandas_index(named_obs):
    frame = ad.AnnData(np.zeros((3, 2), "f4"), obs=named_obs).obs_as("narwhals")
    assert list(nw.maybe_get_index(frame)) == ["AAAC", "AAAG", "AAAT"]


@pytest.mark.parametrize("backend", ["polrs", "duckdb"])
def test_to_backend_rejects_unsupported_backend(dataset2d, backend):
    with pytest.raises(ValueError, match=rf"Unsupported DataFrame backend '{backend}'"):
        to_backend(dataset2d, backend)


@pytest.mark.parametrize("lazy", [False, True], ids=["eager", "lazy"])
def test_narwhals_frame_ingests_as_its_native_frame(lazy):
    """A narwhals frame is unwrapped on the way in, never stored as the wrapper."""
    pl = pytest.importorskip("polars")
    frame = pl.DataFrame({"obs_names": ["c0", "c1", "c2"], "g": [1, 2, 3]})
    adata = ad.AnnData(
        np.zeros((3, 2), "f4"),
        obs=nw.from_native(frame.lazy() if lazy else frame, eager_only=not lazy),
    )
    assert isinstance(adata.obs, pl.DataFrame)
    assert adata.obs_names.tolist() == ["c0", "c1", "c2"]


def test_narwhals_pandas_frame_ingests_by_the_pandas_path(named_obs):
    """Unwrapping dispatches again, so a pandas-backed frame keeps its index."""
    adata = ad.AnnData(
        np.zeros((3, 2), "f4"), obs=nw.from_native(named_obs, eager_only=True)
    )
    expected = ad.AnnData(np.zeros((3, 2), "f4"), obs=named_obs)
    assert isinstance(adata.obs, pd.DataFrame)
    assert list(adata.obs.columns) == list(expected.obs.columns)
    assert adata.obs.index.equals(expected.obs.index)


def test_narwhals_frame_in_obsm_is_unwrapped(named_obs):
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(np.zeros((3, 2), "f4"), obs=named_obs)
    adata.obsm["k"] = nw.from_native(
        pl.DataFrame({"obs_names": list(named_obs.index), "x": [1, 2, 3]}),
        eager_only=True,
    )
    assert isinstance(adata.obsm["k"], pl.DataFrame)


def test_narwhals_frame_in_uns_is_written(tmp_path, diskfmt, named_obs):
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(np.zeros((3, 2), "f4"), obs=named_obs)
    adata.uns["frame"] = nw.from_native(pl.DataFrame({"x": [1, 2, 3]}), eager_only=True)
    path = tmp_path / f"a.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)
    back = getattr(ad.io, f"read_{diskfmt}")(path)
    assert back.uns["frame"]["x"].tolist() == [1, 2, 3]


def test_narwhals_plugin_entry_point_registered(dataset2d):
    from importlib.metadata import entry_points

    from anndata._core import xarray as xarray_module

    eps = {e.name: e for e in entry_points(group="narwhals.plugins")}
    assert "anndata" in eps
    mod = eps["anndata"].load()
    assert mod is xarray_module
    assert mod.NATIVE_PACKAGE == "anndata"
    assert mod.is_native(dataset2d)
    assert not mod.is_native(pd.DataFrame({"a": [1]}))
    assert callable(mod.__narwhals_namespace__)


@pytest.mark.parametrize("backend", ["cudf", "modin"])
def test_pandas_like_backends_route_via_arrow(backend):
    impl = nw.Implementation.from_backend(backend)
    assert impl is not nw.Implementation.UNKNOWN
    assert impl is not nw.Implementation.PANDAS


def test_to_backend_empty_obs():
    pytest.importorskip("polars")

    zero_col = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame(index=pd.Index(["c0", "c1", "c2"], name="obs_names"))
        )
    )
    out = to_backend(zero_col, "polars", dim="obs")
    assert out.shape == (3, 1)
    assert out.columns == ["obs_names"]

    zero_row = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame({"g": pd.Series([], dtype="int64")}).rename_axis("obs_names")
        )
    )
    out0 = to_backend(zero_row, "polars", dim="obs")
    assert out0.shape[0] == 0
    assert set(out0.columns) >= {"g"}


def test_dataframe_like_isinstance(named_obs, dataset2d):
    assert isinstance(named_obs, DataFrameLike)
    assert isinstance(dataset2d, DataFrameLike)
    pl = pytest.importorskip("polars")
    pa = pytest.importorskip("pyarrow")
    assert isinstance(pl.DataFrame({"a": [1]}), DataFrameLike)
    assert isinstance(pa.table({"a": [1]}), DataFrameLike)


def test_from_backend_passthrough(named_obs, dataset2d):
    assert _from_backend(named_obs) is named_obs
    assert _from_backend(dataset2d) is dataset2d


def test_dataframe_x_index_sets_default_axis_names():
    x = pd.DataFrame(np.ones((2, 1)), index=["a", "b"], columns=["g"])
    adata = ad.AnnData(x)
    assert adata.obs_names.tolist() == ["a", "b"]


def test_from_backend_foreign_restores_index(named_obs):
    pytest.importorskip("polars")

    pl_frame = to_backend(named_obs, "polars", dim="obs")
    back = _from_backend(pl_frame, dim="obs")
    assert isinstance(back, DataFrameLike)
    assert back.index.name == "obs_names"
    assert list(back.index) == ["AAAC", "AAAG", "AAAT"]
    assert "obs_names" not in back.columns
    assert back["cell_type"].tolist() == ["A", "B", "A"]


def test_indexed_frame_helpers(named_obs):
    frame = named_obs.copy()
    assert axis_index(frame, dim="obs") is frame.index
    assert frame_annotation_columns(frame, dim="obs") == [
        "cell_type",
        "n_genes",
    ]
    assert column_backed_dim(frame) is None
    assert relabel_axis(frame, source="obs", target="var") is frame

    renamed = set_axis_index(frame, pd.Index(["x", "y", "z"]), dim="obs")
    assert renamed is frame
    assert frame.index.tolist() == ["x", "y", "z"]
    pd.testing.assert_frame_equal(
        subset_frame(frame, (np.array([2, 0]), slice(None))), frame.iloc[[2, 0], :]
    )
    assert frame_equal(frame, frame.copy())
    assert not frame_equal(frame, frame.iloc[:2])
    assert copy_frame(frame) is not frame


@pytest.mark.parametrize("backend", ["polars", "pyarrow"])
def test_indexless_frame_helpers(backend):
    frame = make_native_frame(backend, {"obs_names": ["a", "b"], "value": [1, 2]})
    assert axis_index(frame, dim="obs").tolist() == ["a", "b"]
    assert frame_annotation_columns(frame, dim="obs") == ["value"]
    assert column_backed_dim(frame) == "obs"
    assert relabel_axis(frame, source="obs", target="obs") is frame

    renamed = set_axis_index(frame, pd.Index(["x", "y"]), dim="obs")
    assert axis_index(renamed, dim="obs").tolist() == ["x", "y"]
    subset = subset_frame(frame, (np.array([1]), slice(None)))
    assert axis_index(subset, dim="obs").tolist() == ["b"]
    assert frame_equal(frame, copy_frame(frame))
    assert not frame_equal(frame, make_native_frame(backend, {"value": [1, 2]}))
    assert not frame_equal(frame, object())
    converted = to_backend(
        frame,
        "pyarrow" if backend == "polars" else "pandas",
        dim="obs",
    )
    assert axis_index(converted, dim="obs").tolist() == ["a", "b"]

    relabeled = relabel_axis(frame, source="obs", target="var")
    assert column_backed_dim(relabeled) == "var"
    missing = make_native_frame(backend, {"value": [1]})
    with pytest.raises(ValueError, match="missing required column"):
        axis_index(missing, dim="obs")
    with pytest.raises(ValueError, match="missing required column"):
        relabel_axis(missing, source="obs", target="var")
    collision = make_native_frame(backend, {"obs_names": ["a"], "var_names": ["g"]})
    with pytest.raises(ValueError, match="reserved column 'var_names'"):
        relabel_axis(collision, source="obs", target="var")


def test_copy_frame_isolates_polars_mutation():
    pl = pytest.importorskip("polars")
    frame = pl.DataFrame({"obs_names": ["a", "b"], "value": [1, 2]})

    copied = copy_frame(frame)
    copied.insert_column(2, pl.Series("group", ["x", "y"]))

    assert copied.columns == ["obs_names", "value", "group"]
    assert frame.columns == ["obs_names", "value"]


def test_copy_frame_preserves_dataset2d(dataset2d):
    copied = copy_frame(dataset2d)

    assert isinstance(copied, Dataset2D)
    assert copied is not dataset2d


def test_optional_frame_conversion(named_obs):
    pl = pytest.importorskip("polars")
    assert _ingest_as_pandas(np.ones((2, 2))) is None
    assert _ingest_axis_frame(None, dim="obs") is None
    assert _ingest_axis_frame(named_obs, dim="obs") is named_obs
    with pytest.raises(ValueError, match="Cannot convert"):
        ad.AnnData(X=np.zeros((2, 1)), obs=np.ones((2, 2)))

    native = pl.DataFrame({"value": [1, 2]})
    ensured = _ingest_axis_frame(native, dim="obs")
    assert isinstance(ensured, pl.DataFrame)
    assert ensured["obs_names"].to_list() == ["0", "1"]

    converted = _ingest_as_pandas(
        pl.DataFrame({"obs_names": ["a", "b"], "value": [1, 2]}),
        dim="obs",
    )
    assert isinstance(converted, pd.DataFrame)
    assert converted.index.tolist() == ["a", "b"]
    assert _from_backend(native, dim="obs").index.tolist() == [0, 1]


@pytest.mark.parametrize("key", ["obs", "metadata"])
def test_write_elem_foreign_backend(tmp_path, diskfmt, key):
    from anndata.io import read_elem, write_elem

    if diskfmt == "zarr":
        import zarr

        g = zarr.open_group(str(tmp_path / "t.zarr"), mode="w")
    else:
        import h5py

        g = h5py.File(tmp_path / "t.h5ad", "w")
    frame = make_native_frame("polars", {"obs_names": ["c0", "c1"], "value": [1, 2]})
    write_elem(g, key, frame)
    back = read_elem(g[key])
    assert isinstance(back, pd.DataFrame)
    assert back["value"].tolist() == [1, 2]
    if key == "obs":
        assert back.index.tolist() == ["c0", "c1"]
        assert back.index.name == "obs_names"
    else:
        assert back.index.tolist() == [0, 1]
        assert back["obs_names"].tolist() == ["c0", "c1"]


def test_native_polars_without_names_column_gets_default_names():
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=pl.DataFrame({"kind": ["a", "b"]}))
    assert isinstance(adata.obs, pl.DataFrame)
    assert adata.obs_names.tolist() == ["0", "1"]
    assert adata.obs["obs_names"].to_list() == ["0", "1"]


def test_native_polars_numeric_names_are_stringified():
    pl = pytest.importorskip("polars")
    with pytest.warns(
        ad.ImplicitModificationWarning, match="Transforming to str index"
    ):
        adata = ad.AnnData(X=np.zeros((2, 1)), obs=pl.DataFrame({"obs_names": [1, 2]}))
    assert adata.obs_names.tolist() == ["1", "2"]
    assert adata.obs["obs_names"].to_list() == ["1", "2"]


def test_native_polars_names_make_unique():
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(
        X=np.zeros((3, 3)),
        obs=pl.DataFrame({"obs_names": ["a", "a", "b"]}),
        var=pl.DataFrame({"var_names": ["g", "g", "h"]}),
    )
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    assert adata.obs_names.tolist() == ["a", "a-1", "b"]
    assert adata.var_names.tolist() == ["g", "g-1", "h"]


@pytest.mark.parametrize("backend", ["polars", "pyarrow"])
def test_native_frame_concat_preserves_backend_and_names(backend):
    left = ad.AnnData(
        X=np.ones((1, 2)),
        obs=make_native_frame(backend, {"obs_names": ["a"], "group": ["x"]}),
        var=make_native_frame(backend, {"var_names": ["g0", "g1"], "kind": ["A", "B"]}),
    )
    right = ad.AnnData(
        X=np.ones((1, 2)),
        obs=make_native_frame(backend, {"obs_names": ["b"], "score": [2]}),
        var=make_native_frame(backend, {"var_names": ["g0", "g1"], "kind": ["A", "B"]}),
    )
    left.obsm["meta"] = make_native_frame(backend, {"obs_names": ["a"], "value": [1]})
    right.obsm["meta"] = make_native_frame(backend, {"obs_names": ["b"], "value": [2]})

    result = ad.concat(
        [left, right], join="outer", merge="same", label="source", keys=["l", "r"]
    )

    assert type(result.obs).__module__.split(".")[0] == backend
    assert type(result.var).__module__.split(".")[0] == backend
    assert type(result.obsm["meta"]).__module__.split(".")[0] == backend
    assert result.obs_names.tolist() == ["a", "b"]
    assert result.var_names.tolist() == ["g0", "g1"]
    obs = nw.from_native(result.obs)
    assert obs["group"].to_list() == ["x", None]
    assert obs["score"].to_list() == [None, 2]
    assert obs["source"].to_list() == ["l", "r"]
    meta = nw.from_native(result.obsm["meta"])
    assert meta["obs_names"].to_list() == ["a", "b"]
    assert meta["value"].to_list() == [1, 2]


@pytest.mark.parametrize("backend", ["polars", "pyarrow"])
def test_native_frame_memory_transpose_writeability_and_io(backend, tmp_path, diskfmt):
    from anndata.acc import A

    adata = ad.AnnData(
        X=np.arange(6).reshape(3, 2),
        obs=make_native_frame(
            backend,
            {"obs_names": ["a", "b", "c"], "group": ["x", "y", "x"]},
        ),
        var=make_native_frame(backend, {"var_names": ["g0", "g1"], "kind": ["A", "B"]}),
    )
    adata.obsm["obs_meta"] = make_native_frame(
        backend,
        {"obs_names": ["a", "b", "c"], "value": [1, 2, 3]},
    )
    adata.varm["var_meta"] = make_native_frame(
        backend, {"var_names": ["g0", "g1"], "value": [4, 5]}
    )

    view = adata[["c", "a"], ["g1"]]
    copied = adata.copy()
    in_memory = adata.to_memory(copy=True)
    transposed = adata.T

    assert type(view.obs).__module__.split(".")[0] == backend
    assert view.obs_names.tolist() == ["c", "a"]
    assert view.var_names.tolist() == ["g1"]
    assert type(copied.obs).__module__.split(".")[0] == backend
    assert type(in_memory.obs).__module__.split(".")[0] == backend
    assert in_memory.obs_names.tolist() == ["a", "b", "c"]
    assert adata[A.obs["group"]].tolist() == ["x", "y", "x"]
    assert type(transposed.obs).__module__.split(".")[0] == backend
    assert type(transposed.var).__module__.split(".")[0] == backend
    assert transposed.obs_names.tolist() == ["g0", "g1"]
    assert transposed.var_names.tolist() == ["a", "b", "c"]
    assert nw.from_native(transposed.obsm["var_meta"])["obs_names"].to_list() == [
        "g0",
        "g1",
    ]
    assert nw.from_native(transposed.varm["obs_meta"])["var_names"].to_list() == [
        "a",
        "b",
        "c",
    ]
    assert adata.obsm.to_df()["obs_meta1"].tolist() == [1, 2, 3]
    adata.raw = adata
    assert type(adata.raw.var).__module__.split(".")[0] == backend
    raw = Raw(adata)
    assert type(raw.var).__module__.split(".")[0] == backend
    assert "kind" in str(raw)
    with pytest.warns(FutureWarning):
        assert adata.obs_keys() == ["group"]
    with pytest.warns(FutureWarning):
        assert adata.var_keys() == ["kind"]
    adata.obs_names = ["x", "y", "z"]
    assert nw.from_native(adata.obs)["obs_names"].to_list() == ["x", "y", "z"]
    assert not adata.unwriteable()

    path = tmp_path / f"native-{backend}.{diskfmt}"
    getattr(adata, f"write_{diskfmt}")(path)
    roundtripped = getattr(ad, f"read_{diskfmt}")(path)
    assert roundtripped.obs_names.tolist() == ["x", "y", "z"]
    assert roundtripped.var_names.tolist() == ["g0", "g1"]
    assert roundtripped.obs["group"].tolist() == ["x", "y", "x"]


def test_native_frame_rename_categories_has_explicit_boundary():
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pl.DataFrame({"obs_names": ["a", "b"], "group": ["x", "y"]}),
        var=pl.DataFrame({"var_names": ["g0", "g1"], "kind": ["x", "y"]}),
    )
    for key in ("group", "kind"):
        with pytest.raises(NotImplementedError, match="pandas-backed"):
            adata.rename_categories(key, ["left", "right"])


def test_dataset2d_obsm_validation_and_index_resync():
    adata = ad.AnnData(X=np.zeros((3, 1)))
    sidecar = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame({"score": [1, 2, 3]}, index=adata.obs_names)
        )
    )
    adata.obsm["sidecar"] = sidecar

    new_names = pd.Index(["a", "b", "c"], name="cells")
    adata.obs_names = new_names
    pd.testing.assert_index_equal(adata.obsm["sidecar"].index, new_names)

    mismatched = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame({"score": [1, 2, 3]}, index=["x", "y", "z"])
        )
    )
    with pytest.raises(ValueError, match="index does not match"):
        adata.obsm["mismatched"] = mismatched


def test_concat_normalizes_mixed_dataframe_like_before_equality():
    var = pd.DataFrame({"kind": ["a", "b"]}, index=["g0", "g1"])
    lazy_var = Dataset2D(XDataset.from_dataframe(var))
    left = ad.AnnData(X=np.ones((1, 2)), var=lazy_var)
    right = ad.AnnData(X=np.ones((1, 2)), var=var.copy())

    result = ad.concat([left, right], merge="same")
    pd.testing.assert_frame_equal(result.var, var)


def test_dask_obs_is_stored_as_a_lazy_dataset2d(dask_obs):
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dask_obs)

    assert isinstance(adata.obs, Dataset2D)
    assert adata.obs.is_lazy
    assert adata.obs_names.tolist() == ["AAAC", "AAAG", "AAAT"]
    assert list(adata.obs.columns) == ["cell_type", "n_genes"]
    # Categoricals cannot stay dask arrays without losing their categories, so they are read;
    # everything else is still deferred.
    assert isinstance(adata.obs["cell_type"].dtype, pd.CategoricalDtype)
    assert isinstance(nw.from_native(adata.obs), nw.LazyFrame)


def test_dask_obs_preserves_every_dtype_family():
    """Each dtype needs a different representation, since xarray can hold only some eagerly."""
    dd = pytest.importorskip("dask.dataframe")
    pdf = pd.DataFrame(
        {
            "text": ["x", "y", "x"],
            "nullable_text": pd.array(["a", None, "c"], dtype="string"),
            "kind": pd.Categorical(["p", "q", "p"]),
            "count": pd.array([1, 2, None], dtype="Int64"),
            "score": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["c0", "c1", "c2"], name="obs_names"),
    )
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dd.from_pandas(pdf, npartitions=2))

    # Only the plain float column can stay deferred; the rest xarray must hold eagerly.
    assert adata.obs.is_lazy

    in_memory = adata.obs.to_memory()
    assert in_memory["text"].tolist() == ["x", "y", "x"]
    assert in_memory["nullable_text"].tolist() == ["a", pd.NA, "c"]
    assert in_memory["kind"].tolist() == ["p", "q", "p"]
    assert in_memory["count"].tolist() == [1, 2, pd.NA]
    assert in_memory["score"].tolist() == [1.0, 2.0, 3.0]
    assert isinstance(in_memory["kind"].dtype, pd.CategoricalDtype)
    assert in_memory["count"].dtype == "Int64"
    # Strings round trip as a nullable string dtype rather than as `object`.
    assert isinstance(in_memory["text"].dtype, pd.StringDtype)

    # Subsetting has to work column by column, which is where a raw Arrow-backed
    # string array would fail.
    assert adata[["c2", "c0"]].obs.to_memory()["text"].tolist() == ["x", "x"]
    assert adata.obs_as("dask").compute()["text"].tolist() == ["x", "y", "x"]


def test_dask_obs_takes_names_from_the_names_column():
    dd = pytest.importorskip("dask.dataframe")
    frame = dd.from_pandas(
        pd.DataFrame({"obs_names": ["a", "b", "c"], "value": [1, 2, 3]}), npartitions=2
    )
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=frame)

    assert adata.obs_names.tolist() == ["a", "b", "c"]
    assert list(adata.obs.columns) == ["value"]


def test_dask_obs_positional_index_is_stringified():
    dd = pytest.importorskip("dask.dataframe")
    frame = dd.from_pandas(pd.DataFrame({"value": [1, 2, 3]}), npartitions=1)
    with pytest.warns(
        ad.ImplicitModificationWarning, match="Transforming to str index"
    ):
        adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=frame)

    assert adata.obs_names.tolist() == ["0", "1", "2"]


def test_dask_obs_supports_views_and_copies(dask_obs):
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dask_obs)

    view = adata[["AAAT", "AAAC"]]
    assert isinstance(view.obs, Dataset2D)
    assert view.obs_names.tolist() == ["AAAT", "AAAC"]
    assert view.obs.to_memory()["n_genes"].tolist() == [30, 10]

    copied = adata.copy()
    assert isinstance(copied.obs, Dataset2D)
    assert copied.obs_names.tolist() == adata.obs_names.tolist()

    assert isinstance(adata.to_memory().obs, pd.DataFrame)


def test_dask_obs_shape_mismatch_is_reported(dask_obs):
    with pytest.raises(ValueError, match="must have as many rows"):
        ad.AnnData(X=np.zeros((2, 2), "f4"), obs=dask_obs)


def test_dask_obsm_is_ingested(dask_obs, named_obs):
    dd = pytest.importorskip("dask.dataframe")
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=named_obs)
    adata.obsm["meta"] = dd.from_pandas(
        pd.DataFrame({"score": [1, 2, 3]}, index=named_obs.index), npartitions=2
    )

    assert isinstance(adata.obsm["meta"], Dataset2D)
    assert adata.obsm["meta"].to_memory()["score"].tolist() == [1, 2, 3]


def test_polars_lazyframe_is_collected_on_ingest():
    pl = pytest.importorskip("polars")
    frame = pl.LazyFrame({"obs_names": ["a", "b"], "value": [1, 2]})
    adata = ad.AnnData(X=np.zeros((2, 2), "f4"), obs=frame)

    assert isinstance(adata.obs, pl.DataFrame)
    assert adata.obs_names.tolist() == ["a", "b"]


def test_obs_as_dask_stays_lazy(dask_obs):
    dd = pytest.importorskip("dask.dataframe")
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dask_obs)

    out = adata.obs_as("dask")
    assert isinstance(out, dd.DataFrame)
    computed = out.compute()
    assert computed.index.tolist() == ["AAAC", "AAAG", "AAAT"]
    assert computed["n_genes"].tolist() == [10, 20, 30]


def test_to_backend_dask_from_eager_frames(named_obs):
    dd = pytest.importorskip("dask.dataframe")

    from_pandas = to_backend(named_obs, "dask", dim="obs")
    assert isinstance(from_pandas, dd.DataFrame)
    assert from_pandas.compute().index.tolist() == ["AAAC", "AAAG", "AAAT"]

    pl = pytest.importorskip("polars")
    from_polars = to_backend(
        pl.DataFrame({"obs_names": ["a", "b"], "value": [1, 2]}), "dask", dim="obs"
    )
    assert from_polars.compute().index.tolist() == ["a", "b"]


def test_adapt_rejects_an_uningested_lazy_frame():
    dd = pytest.importorskip("dask.dataframe")
    frame = dd.from_pandas(
        pd.DataFrame({"obs_names": ["a"], "value": [1]}), npartitions=1
    )

    with pytest.raises(TypeError, match="defers evaluation"):
        axis_index(frame, dim="obs")


def test_dask_obs_survives_concat(dask_obs):
    dd = pytest.importorskip("dask.dataframe")
    other = dd.from_pandas(
        pd.DataFrame(
            {"cell_type": pd.Categorical(["B"]), "n_genes": [40]},
            index=pd.Index(["TTTT"], name="obs_names"),
        ),
        npartitions=1,
    )
    left = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dask_obs)
    right = ad.AnnData(X=np.zeros((1, 2), "f4"), obs=other)

    result = ad.concat([left, right])

    assert isinstance(result.obs, Dataset2D)
    assert result.obs_names.tolist() == ["AAAC", "AAAG", "AAAT", "TTTT"]
    assert result.obs.to_memory()["n_genes"].tolist() == [10, 20, 30, 40]


def test_dask_obs_reports_itself_unwriteable(dask_obs, tmp_path, diskfmt):
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=dask_obs)

    assert adata.unwriteable()
    with pytest.raises(NotImplementedError, match="Dataset2D not supported"):
        getattr(adata, f"write_{diskfmt}")(tmp_path / f"a.{diskfmt}")


@pytest.mark.parametrize("backend", ["polars", "pyarrow"])
def test_frame_column_gives_numpy_for_index_less_frames(backend):
    frame = make_native_frame(backend, {"obs_names": ["a", "b"], "value": [1, 2]})

    column = frame_column(frame, "value")
    assert isinstance(column, np.ndarray)
    assert column.tolist() == [1, 2]


def test_frame_column_gives_the_pandas_array_for_indexed_frames(named_obs):
    column = frame_column(named_obs, "cell_type")

    assert isinstance(column, pd.Categorical)
    assert column.tolist() == ["A", "B", "A"]


def test_frame_column_gives_the_variable_for_a_dataset2d(dataset2d, df):
    column = frame_column(dataset2d, df.columns[0])

    assert isinstance(column, XVariable)
    assert column.values.tolist() == df[df.columns[0]].tolist()


def test_true_axis_index_matches_axis_index_when_not_deferred(named_obs):
    pd.testing.assert_index_equal(
        true_axis_index(named_obs, dim="obs"), axis_index(named_obs, dim="obs")
    )


def test_true_axis_index_prefers_out_of_memory_names():
    """A backed frame indexes by a stand-in, but reports the names it defers to."""
    ds = Dataset2D(
        XDataset(
            {"obs_names": ("idx", ["AAAC", "AAAG"]), "n": ("idx", [1, 2])},
            coords={"idx": [0, 1]},
        )
    )
    ds.true_index_dim = "obs_names"

    assert axis_index(ds, dim="obs").tolist() == [0, 1]
    assert true_axis_index(ds, dim="obs").tolist() == ["AAAC", "AAAG"]


@pytest.mark.parametrize("backend", ["polars", "pyarrow"])
def test_native_backend_names_only_index_less_representations(
    backend, named_obs, dataset2d
):
    frame = make_native_frame(backend, {"obs_names": ["a"], "value": [1]})

    assert native_backend(frame) is nw.Implementation.from_backend(backend)
    assert native_backend(named_obs) is None
    assert native_backend(dataset2d) is None
