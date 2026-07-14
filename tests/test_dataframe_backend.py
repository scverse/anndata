from __future__ import annotations

import narwhals as nw
import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata._core._dataframe_backend import (
    DataFrameLike,
    axis_index,
    column_backed_axis_name,
    copy_frame,
    frame_annotation_columns,
    frame_equal,
    from_backend,
    relabel_axis_identity,
    set_axis_index,
    subset_frame,
    to_backend,
    try_ensure_axis_frame,
    try_from_backend,
)
from anndata._core.raw import Raw
from anndata._core.xarray import Dataset2D
from anndata.compat import XDataset
from anndata.tests.helpers import assert_equal, gen_typed_df

pytest.importorskip("xarray")


@pytest.fixture
def df():
    return gen_typed_df(10)


@pytest.fixture
def dataset2d(df):
    return Dataset2D(XDataset.from_dataframe(df))


@pytest.fixture
def named_obs():
    return pd.DataFrame(
        {"cell_type": pd.Categorical(["A", "B", "A"]), "n_genes": [10, 20, 30]},
        index=pd.Index(["AAAC", "AAAG", "AAAT"], name="obs_names"),
    )


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


def test_narwhals_realises_backed_dataset2d():
    da = pytest.importorskip("dask.array")
    n = 50
    ds = Dataset2D(
        XDataset(
            {
                "val": ("idx", da.from_array(np.arange(n, dtype="float64"), chunks=20)),
                "grp": ("idx", pd.array(["a", "b"] * (n // 2), dtype="category")),
            },
            coords={"idx": [f"c{i}" for i in range(n)]},
        )
    )
    frame = nw.from_native(ds)
    assert isinstance(frame, nw.DataFrame)
    native = frame.to_native()
    assert type(native).__module__.split(".")[0] == "pandas"
    assert "dask" not in type(native["val"].values).__module__
    assert native["val"].tolist() == list(range(n))
    assert isinstance(native["grp"].dtype, pd.CategoricalDtype)


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


def test_to_backend_named_index_identity(named_obs):
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


def test_to_backend_canonicalizes_arbitrarily_named_index():
    pl = pytest.importorskip("polars")
    obs = pd.DataFrame({"kind": ["a", "b"]}, index=pd.Index(["c0", "c1"], name="cells"))
    out = to_backend(obs, "polars", index_name="obs_names")
    assert isinstance(out, pl.DataFrame)
    assert out.columns == ["kind", "obs_names"]
    assert out["obs_names"].to_list() == ["c0", "c1"]


def test_to_backend_rejects_reserved_identity_column_collision():
    pytest.importorskip("polars")
    obs = pd.DataFrame(
        {"obs_names": ["annotation-a", "annotation-b"]},
        index=pd.Index(["c0", "c1"], name="cells"),
    )
    with pytest.raises(ValueError, match="reserved column 'obs_names'"):
        to_backend(obs, "polars", index_name="obs_names")


@pytest.mark.parametrize("backend", ["polrs", "duckdb"])
def test_to_backend_rejects_unsupported_backend(dataset2d, backend):
    with pytest.raises(ValueError, match=rf"Unsupported DataFrame backend '{backend}'"):
        to_backend(dataset2d, backend)


def test_narwhals_plugin_entry_point_registered():
    from importlib.metadata import entry_points

    eps = {e.name: e for e in entry_points(group="narwhals.plugins")}
    assert "anndata" in eps
    mod = eps["anndata"].load()
    assert mod.NATIVE_PACKAGE == "anndata"
    assert callable(mod.is_native)
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
    out = to_backend(zero_col, "polars", index_name="obs_names")
    assert out.shape == (3, 1)
    assert out.columns == ["obs_names"]

    zero_row = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame({"g": pd.Series([], dtype="int64")}).rename_axis("obs_names")
        )
    )
    out0 = to_backend(zero_row, "polars", index_name="obs_names")
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
    assert from_backend(named_obs) is named_obs
    assert from_backend(dataset2d) is dataset2d


def test_dataframe_x_index_sets_default_axis_names():
    x = pd.DataFrame(np.ones((2, 1)), index=["a", "b"], columns=["g"])
    adata = ad.AnnData(x)
    assert adata.obs_names.tolist() == ["a", "b"]


def test_from_backend_foreign_restores_index(named_obs):
    pytest.importorskip("polars")

    pl_frame = to_backend(named_obs, "polars", index_name="obs_names")
    back = from_backend(pl_frame, index_name="obs_names")
    assert isinstance(back, DataFrameLike)
    assert back.index.name == "obs_names"
    assert list(back.index) == ["AAAC", "AAAG", "AAAT"]
    assert "obs_names" not in back.columns
    assert back["cell_type"].tolist() == ["A", "B", "A"]


def test_indexed_frame_helpers(named_obs):
    frame = named_obs.copy()
    assert axis_index(frame, index_name="obs_names") is frame.index
    assert frame_annotation_columns(frame, index_name="obs_names") == [
        "cell_type",
        "n_genes",
    ]
    assert column_backed_axis_name(frame) is None
    assert (
        relabel_axis_identity(frame, source_name="obs_names", target_name="var_names")
        is frame
    )

    renamed = set_axis_index(frame, pd.Index(["x", "y", "z"]), index_name="obs_names")
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
    assert axis_index(frame, index_name="obs_names").tolist() == ["a", "b"]
    assert frame_annotation_columns(frame, index_name="obs_names") == ["value"]
    assert column_backed_axis_name(frame) == "obs_names"
    assert (
        relabel_axis_identity(frame, source_name="obs_names", target_name="obs_names")
        is frame
    )

    renamed = set_axis_index(frame, pd.Index(["x", "y"]), index_name="obs_names")
    assert axis_index(renamed, index_name="obs_names").tolist() == ["x", "y"]
    subset = subset_frame(frame, (np.array([1]), slice(None)))
    assert axis_index(subset, index_name="obs_names").tolist() == ["b"]
    assert frame_equal(frame, copy_frame(frame))
    assert not frame_equal(frame, make_native_frame(backend, {"value": [1, 2]}))
    assert not frame_equal(frame, object())
    converted = to_backend(frame, "pyarrow" if backend == "polars" else "polars")
    assert axis_index(converted, index_name="obs_names").tolist() == ["a", "b"]

    relabeled = relabel_axis_identity(
        frame, source_name="obs_names", target_name="var_names"
    )
    assert column_backed_axis_name(relabeled) == "var_names"
    missing = make_native_frame(backend, {"value": [1]})
    with pytest.raises(ValueError, match="missing required identity column"):
        axis_index(missing, index_name="obs_names")
    with pytest.raises(ValueError, match="missing required identity column"):
        relabel_axis_identity(missing, source_name="obs_names", target_name="var_names")
    collision = make_native_frame(backend, {"obs_names": ["a"], "var_names": ["g"]})
    with pytest.raises(ValueError, match="reserved column 'var_names'"):
        relabel_axis_identity(
            collision, source_name="obs_names", target_name="var_names"
        )


def test_optional_frame_conversion(named_obs):
    pl = pytest.importorskip("polars")
    assert try_from_backend(np.ones((2, 2))) is None
    assert try_ensure_axis_frame(None, index_name="obs_names") is None
    assert try_ensure_axis_frame(named_obs, index_name="obs_names") is named_obs
    with pytest.raises(ValueError, match="Cannot convert"):
        ad.AnnData(X=np.zeros((2, 1)), obs=np.ones((2, 2)))

    native = pl.DataFrame({"value": [1, 2]})
    ensured = try_ensure_axis_frame(native, index_name="obs_names")
    assert isinstance(ensured, pl.DataFrame)
    assert ensured["obs_names"].to_list() == ["0", "1"]

    converted = try_from_backend(
        pl.DataFrame({"obs_names": ["a", "b"], "value": [1, 2]}),
        index_name="obs_names",
    )
    assert isinstance(converted, pd.DataFrame)
    assert converted.index.tolist() == ["a", "b"]
    assert from_backend(native, index_name="obs_names").index.tolist() == [0, 1]


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


def test_native_polars_without_identity_column_gets_default_names():
    pl = pytest.importorskip("polars")
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=pl.DataFrame({"kind": ["a", "b"]}))
    assert isinstance(adata.obs, pl.DataFrame)
    assert adata.obs_names.tolist() == ["0", "1"]
    assert adata.obs["obs_names"].to_list() == ["0", "1"]


def test_native_polars_numeric_identity_is_stringified():
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
def test_native_frame_concat_preserves_backend_and_identity(backend):
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
