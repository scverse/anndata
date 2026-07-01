"""DataFrame-backend layer: narwhals plugin + obs/var backend conversion.

Covers:
- the narwhals plugin (``nw.from_native(Dataset2D)``)
- the ``DataFrameLike`` contract + ``from_backend`` (ingest)
- ``to_backend`` / ``AnnData.obs_as`` / ``AnnData.var_as`` (outgest to any eager backend)
"""

from __future__ import annotations

import narwhals as nw
import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata._core._dataframe_backend import DataFrameLike, from_backend, to_backend
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
    """A small obs with a categorical column and a *named* string index (obs_names)."""
    return pd.DataFrame(
        {"cell_type": pd.Categorical(["A", "B", "A"]), "n_genes": [10, 20, 30]},
        index=pd.Index(["AAAC", "AAAG", "AAAT"], name="obs_names"),
    )


def test_narwhals_from_native_roundtrip(df, dataset2d):
    """``nw.from_native`` accepts a Dataset2D, yielding an eager frame matching to_memory()."""
    frame = nw.from_native(dataset2d)
    assert isinstance(frame, nw.DataFrame)
    assert set(frame.columns) == set(df.columns)
    assert_equal(
        frame.to_native().sort_index(axis=1), dataset2d.to_memory().sort_index(axis=1)
    )


def test_narwhals_op_parity_with_pandas(df, dataset2d):
    """A narwhals op on a Dataset2D matches the same op on the source pandas frame (numeric
    columns, which ``to_memory`` doesn't recast — so this is a real transform, not self-comparison)."""
    num = df.select_dtypes("number").columns.tolist()
    predicate = nw.col(num[0]) >= df[num[0]].median()
    from_ds = nw.from_native(dataset2d).filter(predicate).select(num).to_native()
    from_df = (
        nw.from_native(df, eager_only=True).filter(predicate).select(num).to_native()
    )
    assert_equal(from_ds.reset_index(drop=True), from_df.reset_index(drop=True))


def test_narwhals_index_preserved(dataset2d):
    """The row index survives wrapping, recoverable via maybe_get_index (narwhals has no index)."""
    idx = nw.maybe_get_index(nw.from_native(dataset2d))
    assert_equal(idx, dataset2d.to_memory().index)


def test_narwhals_categorical_preserved(named_obs):
    """Categorical columns (the load-bearing dtype) round-trip through the plugin."""
    ds = Dataset2D(XDataset.from_dataframe(named_obs))
    out = nw.from_native(ds).to_native()
    assert isinstance(out["cell_type"].dtype, pd.CategoricalDtype)
    assert out["cell_type"].tolist() == ["A", "B", "A"]


def test_narwhals_realises_backed_dataset2d():
    """``from_native`` realises a backed/lazy (dask) Dataset2D to eager pandas — no dask survives."""
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
    assert isinstance(frame, nw.DataFrame)  # eager, not a LazyFrame
    native = frame.to_native()
    assert type(native).__module__.split(".")[0] == "pandas"
    assert "dask" not in type(native["val"].values).__module__  # realised
    assert native["val"].tolist() == list(range(n))
    assert isinstance(native["grp"].dtype, pd.CategoricalDtype)


@pytest.mark.parametrize("backend", ["pandas", "pyarrow", "polars"])
def test_to_backend_from_dataset2d(dataset2d, backend):
    """to_backend fans a Dataset2D out to any eager backend (values compared at the Arrow level)."""
    pytest.importorskip(backend)

    native = to_backend(dataset2d, backend)
    assert type(native).__module__.split(".")[0] == backend
    got = nw.from_native(native, eager_only=True).to_arrow()
    src = nw.from_native(dataset2d).to_arrow()
    assert set(dataset2d.columns) <= set(got.column_names)
    for name in dataset2d.columns:
        assert got.column(name).to_pylist() == src.column(name).to_pylist()


def test_to_backend_named_index_identity(named_obs):
    """A named index (obs_names) rides into index-less backends as a column; pandas keeps it as the index."""
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
    """AnnData.obs_as / var_as on an in-memory (pandas-backed) AnnData."""
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


def test_obs_as_read_lazy(tmp_path, diskfmt, named_obs):
    """obs_as works the same when obs is a Dataset2D (read lazily) as when it is pandas."""
    pytest.importorskip("polars")
    p = tmp_path / f"a.{diskfmt}"
    getattr(ad.AnnData(X=np.zeros((3, 2), "f4"), obs=named_obs), f"write_{diskfmt}")(p)
    lazy = ad.experimental.read_lazy(p)
    assert isinstance(lazy.obs, Dataset2D)
    obs_pl = lazy.obs_as("polars")
    assert type(obs_pl).__module__.split(".")[0] == "polars"
    assert "obs_names" in obs_pl.columns


def test_var_as_read_lazy(tmp_path, diskfmt):
    """var_as on a backed Dataset2D (axis 1) carries var_names + categorical values."""
    pytest.importorskip("polars")
    var = pd.DataFrame(
        {"gene_type": pd.Categorical(["tf", "rp", "tf"])},
        index=pd.Index(["G1", "G2", "G3"], name="var_names"),
    )
    p = tmp_path / f"a.{diskfmt}"
    getattr(ad.AnnData(X=np.zeros((2, 3), "f4"), var=var), f"write_{diskfmt}")(p)
    lazy = ad.experimental.read_lazy(p)
    assert isinstance(lazy.var, Dataset2D)
    var_pl = lazy.var_as("polars")
    assert var_pl["var_names"].to_list() == ["G1", "G2", "G3"]
    assert var_pl["gene_type"].to_list() == ["tf", "rp", "tf"]


def test_to_backend_unnamed_index_named(named_obs):
    """A default (unnamed) index surfaces as obs_names/var_names, not ``__index_level_0__``."""
    pytest.importorskip("polars")
    obs = pd.DataFrame({"g": [1, 2, 3]}, index=["c0", "c1", "c2"])  # unnamed index
    adata = ad.AnnData(X=np.zeros((3, 2), "f4"), obs=obs)
    assert adata.obs.index.name is None
    cols = adata.obs_as("polars").columns
    assert "obs_names" in cols
    assert not any(c.startswith("__index_level_") for c in cols)
    assert "var_names" in adata.var_as("polars").columns


def test_to_backend_unknown_backend_raises(dataset2d):
    """An unrecognised backend name is a clear ValueError, not a cryptic AssertionError."""
    with pytest.raises(ValueError, match=r"Unsupported DataFrame backend 'polrs'"):
        to_backend(dataset2d, "polrs")


def test_to_backend_lazy_backend_raises(dataset2d):
    """A *lazy* backend (duckdb/dask/...) is rejected with our clear error — ``to_backend``
    is eager-only — instead of failing deep inside narwhals' ``from_arrow``."""
    with pytest.raises(ValueError, match=r"Unsupported DataFrame backend 'duckdb'"):
        to_backend(dataset2d, "duckdb")


def test_narwhals_plugin_entry_point_registered():
    """The narwhals.plugins entry point is declared and exposes the plugin surface."""
    from importlib.metadata import entry_points

    eps = {e.name: e for e in entry_points(group="narwhals.plugins")}
    assert "anndata" in eps
    mod = eps["anndata"].load()
    assert mod.NATIVE_PACKAGE == "anndata"
    assert callable(mod.is_native)
    assert callable(mod.__narwhals_namespace__)


@pytest.mark.parametrize("backend", ["cudf", "modin"])
def test_pandas_like_backends_route_via_arrow(backend):
    """cuDF/modin are recognised and (being non-pandas) route through Arrow — without needing
    the backend installed."""
    impl = nw.Implementation.from_backend(backend)
    assert impl is not nw.Implementation.UNKNOWN
    assert (
        impl is not nw.Implementation.PANDAS
    )  # takes the Arrow (identity-as-column) branch


def test_to_backend_empty_obs():
    """Empty obs shapes (0 columns; 0 rows) convert cleanly to a backend."""
    pytest.importorskip("polars")

    # 0 columns, with a named index
    zero_col = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame(index=pd.Index(["c0", "c1", "c2"], name="obs_names"))
        )
    )
    out = to_backend(zero_col, "polars", index_name="obs_names")
    assert out.shape == (3, 1)
    assert out.columns == ["obs_names"]

    # 0 rows, with columns
    zero_row = Dataset2D(
        XDataset.from_dataframe(
            pd.DataFrame({"g": pd.Series([], dtype="int64")}).rename_axis("obs_names")
        )
    )
    out0 = to_backend(zero_row, "polars", index_name="obs_names")
    assert out0.shape[0] == 0
    assert set(out0.columns) >= {"g"}


def test_dataframe_like_isinstance(named_obs, dataset2d):
    """pandas + Dataset2D conform to DataFrameLike; index-less backends do not."""
    assert isinstance(named_obs, DataFrameLike)  # pandas
    assert isinstance(dataset2d, DataFrameLike)  # Dataset2D
    pl = pytest.importorskip("polars")
    pa = pytest.importorskip("pyarrow")
    assert not isinstance(pl.DataFrame({"a": [1]}), DataFrameLike)
    assert not isinstance(pa.table({"a": [1]}), DataFrameLike)


def test_from_backend_passthrough(named_obs, dataset2d):
    """Already-conforming frames (pandas, Dataset2D) ingest unchanged (same object)."""
    assert from_backend(named_obs) is named_obs
    assert from_backend(dataset2d) is dataset2d


def test_from_backend_foreign_restores_index(named_obs):
    """A foreign backend ingests to a DataFrameLike (pandas), restoring obs_names to the index;
    round-trips with to_backend."""
    pytest.importorskip("polars")

    pl_frame = to_backend(
        named_obs, "polars", index_name="obs_names"
    )  # obs_names is a column
    back = from_backend(pl_frame, index_name="obs_names")
    assert isinstance(back, DataFrameLike)
    assert back.index.name == "obs_names"
    assert list(back.index) == ["AAAC", "AAAG", "AAAT"]
    assert "obs_names" not in back.columns
    assert back["cell_type"].tolist() == ["A", "B", "A"]


def test_write_elem_foreign_backend(tmp_path, diskfmt):
    """write_elem normalizes a foreign frame (polars) to pandas before writing — the
    dispatch-miss path."""
    pytest.importorskip("polars")
    import polars as pl

    from anndata.io import read_elem, write_elem

    if diskfmt == "zarr":
        import zarr

        g = zarr.open_group(str(tmp_path / "t.zarr"), mode="w")
    else:
        import h5py

        g = h5py.File(tmp_path / "t.h5ad", "w")
    write_elem(g, "obs", pl.DataFrame({"obs_names": ["c0", "c1"], "ct": ["T", "B"]}))
    back = read_elem(g["obs"])
    assert isinstance(back, pd.DataFrame)
    assert back["ct"].tolist() == ["T", "B"]
    assert "obs_names" in back.columns  # carried as a column (no AnnData axis context)


def test_assign_foreign_backend_to_obs():
    """Ingest: assigning a foreign frame to obs stores pandas with obs_names as the index."""
    pytest.importorskip("polars")
    import polars as pl

    adata = ad.AnnData(np.zeros((3, 2), "f4"))
    adata.obs = pl.DataFrame({
        "obs_names": ["c0", "c1", "c2"],
        "ct": ["T", "B", "T"],
        "n": [1, 2, 3],
    })
    assert isinstance(adata.obs, pd.DataFrame)
    assert list(adata.obs_names) == ["c0", "c1", "c2"]
    assert adata.obs["ct"].tolist() == ["T", "B", "T"]
    assert "obs_names" not in adata.obs.columns
