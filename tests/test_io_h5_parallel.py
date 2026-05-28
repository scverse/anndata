"""Tests for the auto-detected parallel HDF5 chunk-decompression read path.

The parallel path (``anndata._io.h5_parallel``) must be bit-identical to the
serial path on success and silently fall back to serial when ineligible.
"""

from __future__ import annotations

import warnings
from math import ceil
from typing import TYPE_CHECKING, NamedTuple

import h5py
import hdf5plugin
import numcodecs
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata._io import h5_parallel as h5p

if TYPE_CHECKING:
    from pathlib import Path


# -------------------------------------------------------------------------------
# Shared fixtures / helpers
# -------------------------------------------------------------------------------


# gzip is built into libhdf5; zstd and blosc come from hdf5plugin filters.
CODEC_KWARGS: dict[str, dict] = {
    "gzip": dict(compression="gzip"),
    "zstd": dict(hdf5plugin.Zstd(clevel=3)),
    "blosc": dict(
        hdf5plugin.Blosc(cname="lz4", clevel=5, shuffle=hdf5plugin.Blosc.SHUFFLE)
    ),
}


def _csr(
    n_obs: int = 2_000, n_var: int = 1_000, density: float = 0.1
) -> sparse.csr_matrix:
    return sparse.random(
        n_obs,
        n_var,
        density=density,
        format="csr",
        dtype=np.float32,
        random_state=np.random.default_rng(0),
    )


def _rechunk_x_with_codec(
    path: Path, codec_kwargs: dict, *, n_subchunks: int = 4
) -> None:
    # Only X's 1-D arrays get the filter (tiny obs/var datasets don't support
    # every plugin); >= 2 chunks keeps them eligible at min_mb=0.
    with h5py.File(path, "r+") as f:
        for sub in ("data", "indices"):
            arr = f[f"X/{sub}"][...]
            del f[f"X/{sub}"]
            chunk = max(1, ceil(arr.shape[0] / n_subchunks))
            f["X"].create_dataset(sub, data=arr, chunks=(chunk,), **codec_kwargs)


def _spy_parallel(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    # Lets a test assert the parallel path ran (else, after a silent fallback,
    # a serial-vs-serial comparison would pass vacuously).
    calls: list[int] = []
    real = h5p._do_parallel_read

    def spy(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(h5p, "_do_parallel_read", spy)
    return calls


def _assert_csr_byte_equal(a: ad.AnnData, b: ad.AnnData) -> None:
    Xa, Xb = a.X.tocsr(), b.X.tocsr()
    np.testing.assert_array_equal(Xa.data, Xb.data)
    np.testing.assert_array_equal(Xa.indices, Xb.indices)
    np.testing.assert_array_equal(Xa.indptr, Xb.indptr)


# -------------------------------------------------------------------------------
# Orchestration / public entry — parallel read is bit-identical to serial
# -------------------------------------------------------------------------------


@pytest.mark.parametrize("codec", list(CODEC_KWARGS), ids=list(CODEC_KWARGS))
def test_parallel_byte_identical_to_serial(
    tmp_path: Path, codec: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    p = tmp_path / f"{codec}.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS[codec])

    with ad.settings.override(parallel_h5_read=False):
        serial = ad.read_h5ad(p)

    calls = _spy_parallel(monkeypatch)
    with ad.settings.override(parallel_h5_read_min_mb=0):
        parallel = ad.read_h5ad(p)

    assert calls, "parallel path did not engage; comparison would be vacuous"
    _assert_csr_byte_equal(serial, parallel)


def test_dense_2d_default_threshold_engages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # 4200 x 4200 f32 ~= 70.6 MB > 64 MiB default: engages with default settings.
    X = np.random.default_rng(0).standard_normal((4200, 4200), dtype=np.float32)
    p = tmp_path / "dense.h5ad"
    ad.AnnData(X=X).write_h5ad(p, compression="gzip")
    with h5py.File(p, "r") as f:
        assert isinstance(f["X"], h5py.Dataset), "expected dense X"
        assert f["X"].chunks is not None

    with ad.settings.override(parallel_h5_read=False):
        serial = ad.read_h5ad(p)

    calls = _spy_parallel(monkeypatch)
    parallel = ad.read_h5ad(p)

    assert calls, "parallel path did not engage at the default threshold"
    np.testing.assert_array_equal(serial.X, parallel.X)


def test_csc_gzip_byte_identical(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    X = sparse.random(
        2_000,
        1_000,
        density=0.1,
        format="csc",
        dtype=np.float32,
        random_state=np.random.default_rng(0),
    )
    p = tmp_path / "csc.h5ad"
    ad.AnnData(X=X).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS["gzip"])

    with ad.settings.override(parallel_h5_read=False):
        serial = ad.read_h5ad(p)

    calls = _spy_parallel(monkeypatch)
    with ad.settings.override(parallel_h5_read_min_mb=0):
        parallel = ad.read_h5ad(p)

    assert calls
    Xa, Xb = serial.X.tocsc(), parallel.X.tocsc()
    np.testing.assert_array_equal(Xa.data, Xb.data)
    np.testing.assert_array_equal(Xa.indices, Xb.indices)
    np.testing.assert_array_equal(Xa.indptr, Xb.indptr)


def test_opt_out_disables_parallel(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    p = tmp_path / "x.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS["gzip"])

    calls = _spy_parallel(monkeypatch)
    with ad.settings.override(parallel_h5_read=False, parallel_h5_read_min_mb=0):
        ad.read_h5ad(p)
    assert calls == []


def test_unexpected_error_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    p = tmp_path / "x.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS["gzip"])

    def boom(*args, **kwargs):
        msg = "synthetic"
        raise RuntimeError(msg)

    monkeypatch.setattr(h5p, "_do_parallel_read", boom)
    # warn_once installs a process-wide ignore filter; reset so pytest.warns sees it.
    with warnings.catch_warnings():
        warnings.resetwarnings()
        with (
            pytest.warns(RuntimeWarning, match="parallel HDF5 read failed"),
            ad.settings.override(parallel_h5_read_min_mb=0),
        ):
            adata = ad.read_h5ad(p)
    with ad.settings.override(parallel_h5_read=False):
        ref = ad.read_h5ad(p)
    _assert_csr_byte_equal(ref, adata)


def test_resize_between_reads_succeeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A worker-count change between reads rebuilds the shared pool; both reads
    # must engage and agree (pre-refactor this could silently fall back).
    p = tmp_path / "x.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS["gzip"])

    calls = _spy_parallel(monkeypatch)
    with ad.settings.override(parallel_h5_read_min_mb=0, parallel_h5_read_workers=2):
        first = ad.read_h5ad(p)
    with ad.settings.override(parallel_h5_read_min_mb=0, parallel_h5_read_workers=4):
        second = ad.read_h5ad(p)
    assert calls, "parallel path did not engage"
    _assert_csr_byte_equal(first, second)


# -------------------------------------------------------------------------------
# Eligibility gating — ineligible datasets fall back to serial
# -------------------------------------------------------------------------------


def test_small_dataset_falls_back(tmp_path: Path) -> None:
    X = sparse.random(
        1_000,
        500,
        density=0.05,
        format="csr",
        dtype=np.float32,
        random_state=np.random.default_rng(0),
    )
    p = tmp_path / "small.h5ad"
    ad.AnnData(X=X).write_h5ad(p, compression="gzip")
    with ad.settings.override(parallel_h5_read=False):
        serial = ad.read_h5ad(p)
    parallel = ad.read_h5ad(p)  # default 64 MiB threshold; too small to engage
    _assert_csr_byte_equal(serial, parallel)
    with h5py.File(p, "r") as f:
        assert h5p._eligible_for_multithreading(f["X/data"], min_mb=64) is None


def test_unchunked_dataset_falls_back(tmp_path: Path) -> None:
    p = tmp_path / "contig.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("d", data=np.arange(2_000, dtype=np.float32))
        assert f["d"].chunks is None
        assert h5p._eligible_for_multithreading(f["d"], min_mb=0) is None


def test_unsupported_codec_falls_back(tmp_path: Path) -> None:
    # lzf is HDF5-builtin but has no numcodecs equivalent; stands in for any
    # future filter without a matching decoder.
    p = tmp_path / "lzf.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset(
            "d",
            data=np.arange(4_000, dtype=np.float32),
            chunks=(1_000,),
            compression="lzf",
        )
        assert h5p._eligible_for_multithreading(f["d"], min_mb=0) is None


def test_multifilter_pipeline_falls_back(tmp_path: Path) -> None:
    # Two-filter pipeline (gzip + shuffle): we accept exactly one filter.
    p = tmp_path / "gzip_shuffle.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset(
            "d",
            data=np.arange(4_000, dtype=np.float32),
            chunks=(1_000,),
            compression="gzip",
            shuffle=True,
        )
        assert h5p._eligible_for_multithreading(f["d"], min_mb=0) is None


def test_standalone_lz4_falls_back(tmp_path: Path) -> None:
    # HDF5 LZ4's multi-block framing isn't parseable by numcodecs.LZ4, so the
    # gate must reject it and the serial path must read it.
    X = _csr()
    p = tmp_path / "lz4.h5ad"
    ad.AnnData(X=X).write_h5ad(p)
    _rechunk_x_with_codec(p, dict(hdf5plugin.LZ4()))
    with h5py.File(p, "r") as f:
        assert h5p._eligible_for_multithreading(f["X/data"], min_mb=0) is None
    Xp = ad.read_h5ad(p).X.tocsr()
    np.testing.assert_array_equal(Xp.data, X.data)
    np.testing.assert_array_equal(Xp.indices, X.indices)


def test_partial_chunk_allocation_falls_back() -> None:
    # Unallocated chunks would leave garbage in the pre-allocated output. Mocked
    # since h5py's DatasetID is a C type.
    from types import SimpleNamespace
    from unittest.mock import MagicMock

    plist = MagicMock()
    plist.get_nfilters.return_value = 1
    plist.get_filter.return_value = (h5p._FILTER_DEFLATE, 0, (), b"")
    dset_id = MagicMock()
    dset_id.get_num_chunks.return_value = 9  # one missing: shape+chunks expect 10
    dset_id.get_create_plist.return_value = plist
    fake = SimpleNamespace(
        chunks=(100_000,),
        shape=(1_000_000,),
        dtype=np.dtype("float32"),
        id=dset_id,
    )
    assert h5p._eligible_for_multithreading(fake, min_mb=0) is None
    dset_id.get_num_chunks.return_value = 10  # all allocated -> eligible
    assert h5p._eligible_for_multithreading(fake, min_mb=0) is not None


# -------------------------------------------------------------------------------
# Step 1 — locate chunks on disk (coalesce extents, split, disjointness guard)
# -------------------------------------------------------------------------------


class _StoreInfo(NamedTuple):
    # Matches what h5py's chunk_iter delivers, so tests don't depend on its
    # internal StoreInfo class.
    chunk_offset: tuple[int, ...]
    filter_mask: int
    byte_offset: int
    size: int


def _meta(offset, byte_off, size, mask=0) -> _StoreInfo:
    return _StoreInfo(offset, mask, byte_off, size)


def test_coalesce_adjacent_chunks_contiguous() -> None:
    metas = [_meta((i * 100,), i * 1000, 1000) for i in range(10)]
    extents = h5p._coalesce_adjacent_chunks(metas)
    assert len(extents) == 1
    assert extents[0].file_start == 0
    assert extents[0].file_end == 10_000
    assert extents[0].length == 10_000


def test_coalesce_adjacent_chunks_unsorted_input() -> None:
    metas = [
        _meta((0,), 9000, 1000),
        _meta((100,), 8000, 1000),
        _meta((200,), 7000, 1000),
    ]
    extents = h5p._coalesce_adjacent_chunks(metas)
    assert len(extents) == 1
    assert extents[0].file_start == 7000
    assert extents[0].file_end == 10_000


def test_coalesce_adjacent_chunks_large_gap_splits() -> None:
    # A gap above _EXTENT_ALLOWED_CHUNK_GAPS_BYTES splits the extent so we don't
    # pread megabytes of free space between fragmented chunks.
    metas = [
        _meta((0,), 0, 1000),
        _meta((100,), h5p._EXTENT_ALLOWED_CHUNK_GAPS_BYTES + 10_000, 1000),
    ]
    extents = h5p._coalesce_adjacent_chunks(metas)
    assert len(extents) == 2


def test_split_extents_for_workers_no_cross_extent_subextents() -> None:
    members_a = [_meta((i * 100,), i * 1000, 1000) for i in range(5)]
    members_b = [_meta(((i + 5) * 100,), 10_000_000 + i * 1000, 1000) for i in range(5)]
    extents = [
        h5p._Extent(0, 5000, members_a),
        h5p._Extent(10_000_000, 10_005_000, members_b),
    ]
    sub_extents = h5p._split_extents_for_workers(extents, n_workers=4, n_chunks=10)
    for sub in sub_extents:
        in_a = sub.file_start >= 0 and sub.file_end <= 5000
        in_b = sub.file_start >= 10_000_000 and sub.file_end <= 10_005_000
        assert in_a ^ in_b, f"sub-extent {sub} crosses an extent boundary"


def test_assert_disjoint_chunks_happy_partition_passes() -> None:
    members = [_meta((i * 100,), i * 1000, 1000) for i in range(10)]
    sub_extents = h5p._split_extents_for_workers(
        [h5p._Extent.from_members(members)], n_workers=4, n_chunks=10
    )
    h5p._assert_disjoint_chunks(sub_extents, n_chunks=10)


def test_assert_disjoint_chunks_overlap_raises() -> None:
    # Two sub-extents referencing the same chunk would race on the output.
    overlapping = [
        h5p._Extent.from_members([_meta((100,), 1000, 1000)]),
        h5p._Extent.from_members([_meta((100,), 2000, 1000)]),
    ]
    with pytest.raises(RuntimeError, match="parallel read partition is invalid"):
        h5p._assert_disjoint_chunks(overlapping, n_chunks=2)


def test_assert_disjoint_chunks_missing_chunk_raises() -> None:
    # Fewer chunks than n_chunks would leave part of the output uninitialised.
    truncated = [h5p._Extent.from_members([_meta((100,), 1000, 1000)])]
    with pytest.raises(RuntimeError, match="parallel read partition is invalid"):
        h5p._assert_disjoint_chunks(truncated, n_chunks=2)


def test_disjoint_guard_triggers_serial_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    p = tmp_path / "x.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(p)
    _rechunk_x_with_codec(p, CODEC_KWARGS["gzip"])

    def force_overlap(extents, n_workers, n_chunks):
        # Duplicate the first chunk so the disjointness guard fails.
        first = extents[0].members[0]
        return [h5p._Extent.from_members([first, first])]

    monkeypatch.setattr(h5p, "_split_extents_for_workers", force_overlap)

    with warnings.catch_warnings():
        warnings.resetwarnings()
        with (
            pytest.warns(RuntimeWarning, match="parallel HDF5 read failed"),
            ad.settings.override(parallel_h5_read_min_mb=0),
        ):
            adata = ad.read_h5ad(p)
    with ad.settings.override(parallel_h5_read=False):
        ref = ad.read_h5ad(p)
    _assert_csr_byte_equal(ref, adata)


# -------------------------------------------------------------------------------
# Step 2 — bulk pread chunk bytes into memory (disk-backed driver required)
# -------------------------------------------------------------------------------


def test_core_driver_falls_back(tmp_path: Path) -> None:
    # The in-memory 'core' driver has no on-disk path to pread() against.
    p = tmp_path / "for_core.h5"
    data = np.random.default_rng(0).standard_normal(4_000, dtype=np.float32)
    with h5py.File(p, "w") as f:
        f.create_dataset("d", data=data, chunks=(1_000,), compression="gzip")
    with (
        h5py.File(p, "r", driver="core") as f,
        pytest.raises(h5p._UnsupportedDriver),
    ):
        h5p._open_for_bulk_read(f["d"])


# -------------------------------------------------------------------------------
# Step 3 — decompress + place into output (chunk -> array slice math)
# -------------------------------------------------------------------------------


def test_filter_mask_chunk_bypassed(tmp_path: Path) -> None:
    # A chunk HDF5 stored raw (filter_mask bit set) must skip the decoder.
    path = tmp_path / "mixed.h5"
    chunks = (1000,)
    shape = (4000,)
    dtype = np.dtype("<f8")
    ref = np.random.default_rng(0).standard_normal(shape).astype(dtype)

    with h5py.File(path, "w") as f:
        ds = f.create_dataset(
            "x", shape=shape, dtype=dtype, chunks=chunks, compression="gzip"
        )
        for i in range(3):
            ds[i * 1000 : (i + 1) * 1000] = ref[i * 1000 : (i + 1) * 1000]
        # Fourth chunk: bypass the gzip pipeline and tell HDF5 we did.
        ds.id.write_direct_chunk((3000,), ref[3000:4000].tobytes(), filter_mask=1)

    with h5py.File(path, "r") as f:
        dset = f["x"]
        metas = h5p._gather_chunk_metadata(dset)
        assert sorted(m.filter_mask for m in metas) == [0, 0, 0, 1]
        out = h5p._do_parallel_read(
            dset, numcodecs.Zlib().decode, n_chunks=len(metas), n_workers=2
        )
        serial = dset[...]
    np.testing.assert_array_equal(serial, out)
    np.testing.assert_array_equal(ref, out)


def test_chunk_assignment_round_trip_reconstructs_dense_array() -> None:
    # (250, 300) / (100, 100) chunks gives a clipped row, column, and corner;
    # edge chunks are zero-padded like HDF5 fill, then reassembled.
    rng = np.random.default_rng(0)
    array_shape = (250, 300)
    chunk_shape = (100, 100)
    arr = rng.standard_normal(array_shape, dtype=np.float32)

    out = np.empty(array_shape, dtype=arr.dtype)
    for i in range(0, array_shape[0], chunk_shape[0]):
        for j in range(0, array_shape[1], chunk_shape[1]):
            chunk = np.zeros(chunk_shape, dtype=arr.dtype)
            valid = arr[i : i + chunk_shape[0], j : j + chunk_shape[1]]
            chunk[: valid.shape[0], : valid.shape[1]] = valid
            array_region, chunk_region = h5p._chunk_assignment(
                (i, j), chunk_shape, array_shape
            )
            out[array_region] = chunk[chunk_region]

    np.testing.assert_array_equal(out, arr)


def test_chunk_assignment_excludes_chunk_padding() -> None:
    array_shape = (150,)
    chunk_shape = (100,)

    interior_chunk = np.arange(100, 200, dtype=np.int32)
    # Edge chunk overhangs by 50: first 50 are real, the rest is fill padding.
    edge_chunk = np.full(chunk_shape, -999, dtype=np.int32)
    edge_chunk[:50] = np.arange(200, 250, dtype=np.int32)

    out = np.full(array_shape, -1, dtype=np.int32)

    array_region, chunk_region = h5p._chunk_assignment((0,), chunk_shape, array_shape)
    out[array_region] = interior_chunk[chunk_region]

    array_region, chunk_region = h5p._chunk_assignment((100,), chunk_shape, array_shape)
    out[array_region] = edge_chunk[chunk_region]

    np.testing.assert_array_equal(out[:100], np.arange(100, 200))
    np.testing.assert_array_equal(out[100:150], np.arange(200, 250))
    assert -999 not in out  # padding never leaked into the destination


# -------------------------------------------------------------------------------
# Codec dispatch — blosc thread configuration is lazy
# -------------------------------------------------------------------------------


def test_blosc_threading_configured_only_when_blosc_used(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # set_nthreads(1) must fire the first time we decode blosc, not on import
    # and not for non-blosc reads.
    gzip_path = tmp_path / "gzip.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(gzip_path)
    _rechunk_x_with_codec(gzip_path, CODEC_KWARGS["gzip"])
    blosc_path = tmp_path / "blosc.h5ad"
    ad.AnnData(X=_csr()).write_h5ad(blosc_path)
    _rechunk_x_with_codec(blosc_path, CODEC_KWARGS["blosc"])

    calls: list[int] = []

    def fake_set_nthreads(n: int) -> None:
        calls.append(n)

    monkeypatch.setattr(numcodecs.blosc, "set_nthreads", fake_set_nthreads)
    # Reset the lazy gate so the spy can observe it.
    monkeypatch.setattr(h5p, "_blosc_threads_configured", False)

    with ad.settings.override(parallel_h5_read_min_mb=0):
        ad.read_h5ad(gzip_path)
        assert calls == []  # gzip read must not touch blosc thread state

        ad.read_h5ad(blosc_path)
        assert calls == [1]  # first blosc read configures it

        ad.read_h5ad(blosc_path)
        assert calls == [1]  # idempotent
