"""Parallel HDF5 chunk decompression for whole-array reads.

By default, HDF5 reads from a single process will use only a single
thread (due to HDF5's design for multi-processing, and h5py's PHIL).
This leads to decompression of all data serializing onto one core.

When eligible (chunked, supported codec, large enough), this module
uses multi-threading to read HDF5 data into memory, in 3 steps:

1. **Locate chunks on disk** — get every chunk's on-disk byte
   offset/size in a single C-side iteration.
2. **Bulk-pread chunk bytes into memory** — spin up a Python threads
   and evenly distribute memory to be read (os.pread()) to the threads.
   This bypasses h5py's process-global lock (PHIL).
3. **Decompress + place** — each worker walks its chunks, decompresses
   them with the relevant GIL-releasing Python decoder, and writes each
   decoded chunk into its destination slice of the output array.

Note that numcodecs (Zarr's codec package is reused here; any codec from
that package can be added here with a few lines of code.
   
Falls back silently to the caller's serial path if any issues occur.
"""

from __future__ import annotations

import math
import os
import threading
from concurrent.futures import wait
from itertools import batched, pairwise
from typing import TYPE_CHECKING, NamedTuple

import numcodecs
import numcodecs.blosc
import numpy as np

from anndata._settings import settings
from anndata.utils import warn_once

from ._pool import PoolManager

if TYPE_CHECKING:
    from collections.abc import Callable

    import h5py

_pool_manager = PoolManager(thread_name_prefix="anndata-h5-parallel")

# HDF5 filter IDs
# https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_z.html
# https://github.com/HDFGroup/hdf5_plugins/blob/master/docs/RegisteredFilterPlugins.md
_FILTER_DEFLATE = 1
_FILTER_BLOSC = 32001
_FILTER_ZSTD = 32015

# HDF5 generally has several-KB gaps between some of it's on-disk chunks.
# These primarily come from the HDF5 index, e.g. in the default case of
# it's B-tree, index nodes are written between occasional chunks.
# So for efficiency of reading memory, we read these extra bytes
# (which we won't use), because it's fastest.
# Note that gaps can also occur from unusual HDF5 file writing patterns, so
# in the case of gaps larger than 256KB, we simply do separate reads excluding
# those gaps to avoid reading excessive unused memory.
_EXTENT_ALLOWED_CHUNK_GAPS_BYTES = 256 * 1024

# file drivers we can pread() against
#   - 'sec2' is the default POSIX driver on Unix.
#   - 'stdio' goes through fopen() but still backs a single regular file.
_DISK_BACKED_DRIVERS: frozenset[str] = frozenset({"sec2", "stdio"})


class _UnsupportedDriver(RuntimeError):
    """Raised when the h5py file uses a non-disk driver (e.g., 'core',
    'mpio', 'ros3', 'fileobj') for which our bulk ``pread()`` approach
    isn't applicable. Caller falls back to serial."""


class _ChunkIterUnavailable(RuntimeError):
    """Raised when ``Dataset.id.chunk_iter`` is not present (older h5py
    < 3.8 / older HDF5 builds). Recommended to update for multi-threading
    (significant performance improvements)."""


# -------------------------------------------------------------------------------
# Codec dispatch (HDF5 filter id -> GIL-releasing decoder)
# -------------------------------------------------------------------------------


_FILTER_TO_DECODER: dict[int, Callable[[bytes], bytes]] = {
    _FILTER_DEFLATE: numcodecs.Zlib().decode, # Note that HDF5 'gzip' is actually zlib
    _FILTER_BLOSC: numcodecs.Blosc().decode,  # inner codec self-described in payload
    _FILTER_ZSTD: numcodecs.Zstd().decode,
}

_blosc_threads_lock = threading.Lock()
_blosc_threads_configured = False


def _ensure_blosc_single_threaded() -> None:
    global _blosc_threads_configured
    if _blosc_threads_configured:
        return
    with _blosc_threads_lock:
        if _blosc_threads_configured:
            return
        numcodecs.blosc.set_nthreads(1)
        _blosc_threads_configured = True


# -------------------------------------------------------------------------------
# Eligibility gating
# -------------------------------------------------------------------------------

def _supported_filter(dset: h5py.Dataset) -> Callable[[bytes], bytes] | None:
    # AnnData does not generally deal with HDF5 files with multiple filter
    # steps (e.g. shuffle + compressor). Single-filter pipelines only.
    plist = dset.id.get_create_plist()
    if plist.get_nfilters() != 1:
        return None
    filter_id = plist.get_filter(0)[0]
    decoder = _FILTER_TO_DECODER.get(filter_id)
    if decoder is None:
        return None
    if filter_id == _FILTER_BLOSC:
        _ensure_blosc_single_threaded()
    return decoder


def _eligible_for_multithreading(
    dset: h5py.Dataset, *, min_mb: float
) -> tuple[Callable[[bytes], bytes], int] | None:
    # Nothing to parallelize if no chunks
    if dset.chunks is None:
        return None
    n_chunks = dset.id.get_num_chunks()
    if n_chunks < 2:
        return None

    # Reject datasets with initialized but unfilled chunks
    expected = math.prod(
        math.ceil(s / c)
        for s, c in zip(dset.shape, dset.chunks, strict=True)
    )
    if n_chunks != expected:
        return None

    nbytes = int(np.prod(dset.shape)) * dset.dtype.itemsize
    if nbytes < min_mb * (1 << 20):
        return None

    # Codec must be supported
    decoder = _supported_filter(dset)
    if decoder is None:
        return None
    return decoder, n_chunks


# -------------------------------------------------------------------------------
# Step 1 — locate chunks on disk
# -------------------------------------------------------------------------------


class _Extent(NamedTuple):
    """A contiguous byte range on disk plus the chunks it covers.

    Produced first by adjacency-coalescing all chunks (one extent per
    contiguous on-disk region), then refined by splitting each extent
    across workers (one ``pread()`` per resulting sub-extent).
    ``members`` are ``StoreInfo`` namedtuples from h5py's ``chunk_iter``
    callback (fields: ``chunk_offset``, ``filter_mask``,
    ``byte_offset``, ``size``), and must be sorted by ``byte_offset``
    and non-overlapping — both hold for HDF5 dataset chunks by
    construction.
    """

    file_start: int  # absolute byte offset in the file
    file_end: int  # exclusive
    members: list  # list[StoreInfo]

    @property
    def length(self) -> int:
        return self.file_end - self.file_start

    @classmethod
    def from_members(cls, members) -> _Extent:
        members = list(members)  # batched() yields tuples; normalize
        return cls(
            file_start=members[0].byte_offset,
            file_end=members[-1].byte_offset + members[-1].size,
            members=members,
        )


def _gather_chunk_metadata(dset: h5py.Dataset) -> list:
    """Enumerate every chunk's metadata in a single C-side iteration."""
    if not hasattr(dset.id, "chunk_iter"):
        msg = (
            "h5py.Dataset.id.chunk_iter not available; need h5py >= 3.8"
            "upgrading is strongly recommended for better performance."
        )
        raise _ChunkIterUnavailable(msg)
    metas: list = []
    dset.id.chunk_iter(metas.append)
    return metas


def _coalesce_adjacent_chunks(metas: list) -> list[_Extent]:
    """\
    This determines the continuous sections of memory to be read (extents).

    This will generally return one extent unless the HDF5 file has gaps
    between chunks greater than ``_EXTENT_ALLOWED_CHUNK_GAPS_BYTES``. See
    above comment on this parameter.
    """
    if not metas:
        return []
    sorted_metas = sorted(metas, key=lambda m: m.byte_offset)
    splits = [
        i
        for i, (prev, cur) in enumerate(pairwise(sorted_metas), start=1)
        if cur.byte_offset - (prev.byte_offset + prev.size)
        > _EXTENT_ALLOWED_CHUNK_GAPS_BYTES
    ]
    boundaries = [0, *splits, len(sorted_metas)]
    return [_Extent.from_members(sorted_metas[a:b]) for a, b in pairwise(boundaries)]


def _split_extents_for_workers(
    extents: list[_Extent], n_workers: int, n_chunks: int
) -> list[_Extent]:
    step = max(1, math.ceil(n_chunks / n_workers))
    return [
        _Extent.from_members(sub)
        for extent in extents
        for sub in batched(extent.members, step)
    ]


def _assert_disjoint_chunks(sub_extents: list[_Extent], n_chunks: int) -> None:
    """Refuse to dispatch unless every chunk appears in exactly one sub-extent;
    otherwise workers would race on the output array."""
    seen: set[tuple[int, ...]] = {
        m.chunk_offset for ext in sub_extents for m in ext.members
    }
    if len(seen) != n_chunks:
        msg = (
            f"parallel read partition is invalid: {len(seen)} unique chunks, "
            f"expected {n_chunks}."
        )
        raise RuntimeError(msg)


# -------------------------------------------------------------------------------
# Step 2 — bulk pread chunk bytes into memory
# -------------------------------------------------------------------------------


def _open_for_bulk_read(dset: h5py.Dataset) -> int:
    """Open a raw read-only file descriptor at the path backing ``dset``.

    This is separate from h5py's fd because h5py routes
    I/O through libhdf5's VFD layer; here we want a read that bypasses
    h5py's PHIL. Multiple read-only fds on the same file
    are universally supported on POSIX and Windows.

    Falls back via ``_UnsupportedDriver`` if the file isn't on a normal
    disk-backed driver (e.g., in-memory, MPI-IO, S3-backed).
    """
    driver = dset.file.driver
    if driver not in _DISK_BACKED_DRIVERS:
        msg = f"file driver {driver!r} is not supported by the bulk-read path"
        raise _UnsupportedDriver(msg)
    fd = os.open(dset.file.filename, os.O_RDONLY)

    # On linux, speed up cold cache (start populating page cache)
    if hasattr(os, "posix_fadvise"):
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_WILLNEED)
        except OSError:  # pragma: no cover - some FSes refuse fadvise
            pass
    return fd


# TODO: before completion of PR, adjust for proper Windows support
def _pread(fd: int, n: int, offset: int) -> bytes:
    """Cross-platform positional read. ``os.pread`` is Unix-only;
    Windows would need lseek+read."""
    if hasattr(os, "pread"):
        return os.pread(fd, n, offset)
    os.lseek(fd, offset, os.SEEK_SET)  # pragma: no cover - Unix in CI
    return os.read(fd, n)  # pragma: no cover - Unix in CI


# -------------------------------------------------------------------------------
# Step 3 — decompress + place into output
# -------------------------------------------------------------------------------


def _chunk_assignment(
    chunk_offset: tuple[int, ...],
    chunk_shape: tuple[int, ...],
    array_shape: tuple[int, ...],
) -> tuple[tuple[slice, ...], tuple[slice, ...]]:
    """Chunked HDF5 is designed such that N-dimensional arrays are split into N-dimensional
    chunks, all of the same size (e.g. tiles).

    Thus, chunks around the edges typically contain extra, unused bytes.

    Here we identify how much of the chunk to keep, and the position of the N-dimensional
    output array where it should reside. This results in the chunks being properly
    re-assembled into the original array.
    """
    array_region: list[slice] = []
    chunk_region: list[slice] = []
    for chunk_start, chunk__len, array_axis_len in zip(
        chunk_offset, chunk_shape, array_shape, strict=True
    ):
        chunk_boundary = min(chunk__len, array_axis_len - chunk_start)
        array_region.append(slice(chunk_start, chunk_start + chunk_boundary))
        chunk_region.append(slice(0, chunk_boundary))
    return tuple(array_region), tuple(chunk_region)


def _process_extent(
    extent: _Extent,
    fd: int,
    decoder: Callable[[bytes], bytes],
    dtype: np.dtype,
    chunk_shape: tuple[int, ...],
    array_shape: tuple[int, ...],
    chunk_elems: int,
    out: np.ndarray,
) -> None:
    """Read one extent with a single ``pread`` and assign each of its
    chunks into the corresponding region of ``out``.

    Workers run concurrently against the same ``out``; disjointness of
    their writes is enforced upstream by :func:`_assert_disjoint_chunks`.
    """
    extent_bytes = _pread(fd, extent.length, extent.file_start)
    if len(extent_bytes) != extent.length:
        msg = (
            f"short pread: got {len(extent_bytes)} bytes, "
            f"expected {extent.length} at offset {extent.file_start}"
        )
        raise OSError(msg)

    extent_view = memoryview(extent_bytes)
    for chunk_info in extent.members:
        rel_offset = chunk_info.byte_offset - extent.file_start
        raw_chunk = extent_view[rel_offset : rel_offset + chunk_info.size]

        # Decode chunk. Skips if HDF5 did not code this chunk.
        decoded = raw_chunk if chunk_info.filter_mask else decoder(raw_chunk)

        # Restore it to it's original shape (C row-major order, same as numpy)
        chunk_array = np.frombuffer(
            decoded, dtype=dtype, count=chunk_elems
        ).reshape(chunk_shape)

        # Put it in the correct position of 'out' (noted in the chunk's metadata)
        array_region, chunk_region = _chunk_assignment(
            chunk_info.chunk_offset, chunk_shape, array_shape
        )
        out[array_region] = chunk_array[chunk_region]


# -------------------------------------------------------------------------------
# Orchestration / public entry
# -------------------------------------------------------------------------------


def _resolve_num_workers(workers: int | None) -> int:
    if workers is not None:
        return workers
    return min(os.cpu_count() or 1, 16)


def _warn_fallback(exc: Exception) -> None:
    msg = (
        "anndata: parallel HDF5 read failed; falling back to serial. "
        f"Reason: {type(exc).__name__}: {exc}"
    )
    warn_once(msg, RuntimeWarning)


def _do_parallel_read(
    dset: h5py.Dataset,
    decoder: Callable[[bytes], bytes],
    n_chunks: int,
    n_workers: int,
) -> np.ndarray:
    metas = _gather_chunk_metadata(dset)
    extents = _coalesce_adjacent_chunks(metas)
    pread_extents = _split_extents_for_workers(extents, n_workers, n_chunks)
    _assert_disjoint_chunks(pread_extents, n_chunks)

    # worker threads will share the same fd; ``os.pread`` is thread-safe.
    fd = _open_for_bulk_read(dset)

    array_shape = dset.shape
    chunk_shape = dset.chunks
    dtype = dset.dtype
    chunk_elems = int(np.prod(chunk_shape))
    out = np.empty(array_shape, dtype=dtype)

    try:
        futs = _pool_manager.distribute_tasks_to_threads(
            n_workers,
            _process_extent,
            [
                (extent, fd, decoder, dtype, chunk_shape, array_shape, chunk_elems, out)
                for extent in pread_extents
            ],
        )

        wait(futs)
        for f in futs:
            f.result()
    finally:
        os.close(fd)
    return out


def parallel_read_full(dset: h5py.Dataset) -> np.ndarray | None:
    """Read a whole ``h5py.Dataset`` into memory using bulk-pread + parallel
    decompression if eligible, else return ``None``.

    Bit-identical to ``dset[...]`` on success. Returns ``None`` (so callers
    can fall back to their serial path) when:

    * ``settings.parallel_h5_read`` is disabled,
    * the dataset is ineligible (see :func:`_eligible_for_multithreading`),
    * the file driver isn't disk-backed (we need a real fd for ``pread``),
    * ``chunk_iter`` isn't available (h5py < 3.8 / older HDF5 build),
    * or any unexpected error occurs during the parallel read.
    """
    if not settings.parallel_h5_read:
        return None
    try:
        info = _eligible_for_multithreading(dset, min_mb=settings.parallel_h5_read_min_mb)
        if info is None:
            return None
        decoder, n_chunks = info
        workers = _resolve_num_workers(settings.parallel_h5_read_workers)
        return _do_parallel_read(dset, decoder, n_chunks, workers)
    except Exception as exc:  # noqa: BLE001
        _warn_fallback(exc)
        return None
