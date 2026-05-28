<!--
Rename this file to <PR_NUMBER>.feat.md when the PR is opened (towncrier
issue_format expects a numeric id). For example: 2463.feat.md
-->
Multi-threaded HDF5 reads in ``read_h5ad``. Chunked datasets
compressed with **gzip / deflate**, **zstd**, or **blosc** (any inner
codec — lz4, zstd, blosclz, snappy) are now read into memory in
parallel: chunk locations are enumerated via
``Dataset.id.chunk_iter``, adjacent chunks are coalesced into on-disk
extents, and worker threads each ``pread()`` a sub-extent (bypassing
h5py's process-global lock) and decompress their chunks into disjoint
regions of the output array using GIL-releasing decoders from
:mod:`numcodecs`.

Bit-identical to the serial path; silently falls back to serial for
ineligible datasets (unsupported or multi-filter pipeline, unchunked,
partially written, smaller than ``parallel_h5_read_min_mb``, non-disk
file driver, missing ``chunk_iter`` on h5py < 3.8) or on any
unexpected error.

Controlled by ``settings.parallel_h5_read`` (default ``True``),
``settings.parallel_h5_read_min_mb`` (default ``64``), and
``settings.parallel_h5_read_workers`` (default ``min(cpu_count, 16)``).
