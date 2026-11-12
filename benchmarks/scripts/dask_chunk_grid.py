from __future__ import annotations

import argparse
import csv
import itertools
import json
import platform
import shutil
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import h5py
import numpy as np
import zarr

import anndata as ad
from anndata.experimental import read_elem_lazy

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

ChunkSpec = tuple[int | None, ...]
StoreType = Literal["h5ad", "zarr"]


DEFAULT_ON_DISK_CHUNKS: tuple[ChunkSpec, ...] = ((256, 1024), (1024, 1024))
DEFAULT_DASK_CHUNKS: tuple[ChunkSpec | None, ...] = (None, (1024, -1), (4096, -1))
DEFAULT_WORKLOADS = (
    "sum_axis0",
    "sum_axis1",
    "normalize_log1p_slice",
    "scanpy_normalize_log1p",
)


@dataclass(frozen=True)
class StoreConfig:
    store_type: StoreType
    on_disk_chunks: ChunkSpec
    zarr_format: int
    zarr_shards: ChunkSpec | None


@dataclass(frozen=True)
class DaskConfig:
    chunks: ChunkSpec | None
    n_workers: int
    threads_per_worker: int
    processes: bool


@dataclass(frozen=True)
class WorkloadConfig:
    shape: tuple[int, int]
    slice_obs: int
    slice_vars: int
    target_sum: float
    memory_limit: str
    dataset_nbytes: int
    run_started_at: str


@dataclass(frozen=True)
class BenchmarkResult:
    run_started_at: str
    python_version: str
    platform: str
    anndata_version: str
    numpy_version: str
    h5py_version: str
    zarr_version: str
    dask_version: str
    distributed_version: str
    scanpy_version: str
    store_type: str
    zarr_format: int | str
    zarr_shards: str
    shape: str
    dataset_nbytes: int
    store_nbytes: int
    on_disk_chunks: str
    dask_chunks_arg: str
    dask_chunksize: str
    workers: int
    threads_per_worker: int
    processes: bool
    workload: str
    repeat: int
    task_count: int
    elapsed_s: float
    result_shape: str
    result_nbytes: int
    driver_max_rss_mb: float | None
    worker_rss_before_mb: float | None
    worker_rss_after_mb: float | None


def parse_chunk_spec(value: str) -> ChunkSpec | None:
    normalized = value.strip().lower()
    if normalized in {"default", "none", "null"}:
        return None
    if normalized == "full":
        return (-1, -1)

    parts = normalized.replace("x", ",").split(",")
    if len(parts) != 2:
        msg = f"Expected two comma-separated dimensions, got {value!r}"
        raise argparse.ArgumentTypeError(msg)

    parsed: list[int | None] = []
    for part in parts:
        part = part.strip()
        if part in {"default", "none", "null"}:
            parsed.append(None)
            continue
        if part == "full":
            parsed.append(-1)
            continue
        try:
            dim = int(part)
        except ValueError as exc:
            msg = f"Chunk dimension {part!r} is not an integer"
            raise argparse.ArgumentTypeError(msg) from exc
        if dim == 0 or dim < -1:
            msg = f"Chunk dimension must be positive, -1, or None, got {dim}"
            raise argparse.ArgumentTypeError(msg)
        parsed.append(dim)
    return tuple(parsed)


def parse_shape(value: str) -> tuple[int, int]:
    parsed = parse_chunk_spec(value)
    if parsed is None or any(dim is None or dim <= 0 for dim in parsed):
        msg = f"Shape must contain two positive integers, got {value!r}"
        raise argparse.ArgumentTypeError(msg)
    return (int(parsed[0]), int(parsed[1]))


def materialize_chunks(chunks: ChunkSpec, shape: tuple[int, int]) -> tuple[int, int]:
    return tuple(
        axis_size if chunk in {-1, None} else min(int(chunk), axis_size)
        for chunk, axis_size in zip(chunks, shape, strict=True)
    )


def format_chunks(chunks: Sequence[int | None] | None) -> str:
    if chunks is None:
        return "default"
    return "x".join("none" if chunk is None else str(chunk) for chunk in chunks)


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def generate_counts(shape: tuple[int, int], seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.poisson(lam=1.0, size=shape).astype(np.float32, copy=False)


def path_nbytes(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(child.stat().st_size for child in path.rglob("*") if child.is_file())


def package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "not-installed"


def runtime_metadata() -> dict[str, str]:
    return {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "anndata_version": package_version("anndata"),
        "numpy_version": package_version("numpy"),
        "h5py_version": package_version("h5py"),
        "zarr_version": package_version("zarr"),
        "dask_version": package_version("dask"),
        "distributed_version": package_version("distributed"),
        "scanpy_version": package_version("scanpy"),
    }


def store_path(base_dir: Path, config: StoreConfig) -> Path:
    chunk_label = format_chunks(config.on_disk_chunks)
    if config.store_type == "h5ad":
        return base_dir / f"X_h5ad_chunks-{chunk_label}.h5ad"

    shard_label = (
        f"_shards-{format_chunks(config.zarr_shards)}"
        if config.zarr_shards is not None
        else ""
    )
    return (
        base_dir
        / f"X_zarr-v{config.zarr_format}_chunks-{chunk_label}{shard_label}.zarr"
    )


def write_store(
    path: Path,
    config: StoreConfig,
    shape: tuple[int, int],
    seed: int,
    *,
    force: bool,
) -> None:
    if path.exists():
        if not force:
            return
        remove_path(path)

    data = generate_counts(shape, seed)
    on_disk_chunks = materialize_chunks(config.on_disk_chunks, shape)
    if config.store_type == "h5ad":
        with h5py.File(path, "w") as f:
            ad.io.write_elem(f, "X", data, dataset_kwargs={"chunks": on_disk_chunks})
        return

    dataset_kwargs: dict[str, object] = {"chunks": on_disk_chunks}
    if config.zarr_shards is not None:
        dataset_kwargs["shards"] = materialize_chunks(config.zarr_shards, shape)

    root = zarr.open_group(path, mode="w", zarr_format=config.zarr_format)
    with ad.settings.override(auto_shard_zarr_v3=False):
        ad.io.write_elem(root, "X", data, dataset_kwargs=dataset_kwargs)
    zarr.consolidate_metadata(root.store)


@contextmanager
def open_x(path: Path, store_type: StoreType) -> Iterator[object]:
    if store_type == "h5ad":
        with h5py.File(path, "r") as f:
            yield f["X"]
        return

    root = zarr.open_group(path, mode="r")
    yield root["X"]


def max_rss_mb() -> float | None:
    try:
        import resource
    except ImportError:
        return None

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() == "Darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def _process_rss_bytes() -> int | None:
    try:
        import psutil
    except ImportError:
        return None
    return psutil.Process().memory_info().rss


def worker_rss_mb(client) -> float | None:
    values = client.run(_process_rss_bytes).values()
    rss_values = [value for value in values if value is not None]
    if not rss_values:
        return None
    return sum(rss_values) / (1024 * 1024)


@contextmanager
def dask_client(config: DaskConfig, memory_limit: str) -> Iterator[object]:
    try:
        from dask.distributed import Client, LocalCluster
    except ImportError as exc:
        msg = "dask.distributed is required; run with `uv run --group test-min python`."
        raise RuntimeError(msg) from exc

    cluster = LocalCluster(
        n_workers=config.n_workers,
        threads_per_worker=config.threads_per_worker,
        processes=config.processes,
        dashboard_address=None,
        memory_limit=memory_limit,
    )
    client = Client(cluster)
    try:
        yield client
    finally:
        client.close()
        cluster.close()


def build_workload(
    array,
    workload: str,
    *,
    slice_obs: int,
    slice_vars: int,
    target_sum: float,
):
    import dask.array as da

    if workload == "sum_axis0":
        return array.sum(axis=0)
    if workload == "sum_axis1":
        return array.sum(axis=1)
    if workload == "subset_mean":
        return array[:slice_obs, :slice_vars].mean(axis=0)
    if workload == "normalize_log1p_slice":
        subset = array[:slice_obs, :]
        row_sums = subset.sum(axis=1)
        safe_row_sums = da.where(row_sums == 0, 1, row_sums)
        normalized = subset / safe_row_sums[:, None] * target_sum
        return da.log1p(normalized[:, :slice_vars]).mean(axis=0)
    if workload == "scanpy_normalize_log1p":
        try:
            import scanpy as sc
        except ImportError as exc:
            msg = (
                "scanpy_normalize_log1p requires scanpy; run with "
                "`uv run --group test-min python`."
            )
            raise RuntimeError(msg) from exc

        adata = ad.AnnData(X=array)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        return adata.X[:, :slice_vars].mean(axis=0)

    msg = f"Unknown workload {workload!r}"
    raise ValueError(msg)


def task_count(collection) -> int:
    return len(collection.__dask_graph__())


def result_nbytes(result: object) -> int:
    return int(getattr(result, "nbytes", 0))


def result_shape(result: object) -> str:
    shape = getattr(result, "shape", ())
    return format_chunks(tuple(int(dim) for dim in shape))


def run_one(
    path: Path,
    store_config: StoreConfig,
    dask_config: DaskConfig,
    *,
    workload: str,
    repeat: int,
    workload_config: WorkloadConfig,
) -> BenchmarkResult:
    with (
        dask_client(dask_config, memory_limit=workload_config.memory_limit) as client,
        open_x(path, store_config.store_type) as elem,
    ):
        metadata = runtime_metadata()
        array = read_elem_lazy(elem, chunks=dask_config.chunks)
        slice_obs = min(workload_config.slice_obs, array.shape[0])
        slice_vars = min(workload_config.slice_vars, array.shape[1])
        workload_array = build_workload(
            array,
            workload,
            slice_obs=slice_obs,
            slice_vars=slice_vars,
            target_sum=workload_config.target_sum,
        )

        before_worker_rss = worker_rss_mb(client)
        started = time.perf_counter()
        computed = workload_array.compute()
        elapsed_s = time.perf_counter() - started
        after_worker_rss = worker_rss_mb(client)

        return BenchmarkResult(
            run_started_at=workload_config.run_started_at,
            **metadata,
            store_type=store_config.store_type,
            zarr_format=store_config.zarr_format
            if store_config.store_type == "zarr"
            else "n/a",
            zarr_shards=format_chunks(store_config.zarr_shards),
            shape=format_chunks(workload_config.shape),
            dataset_nbytes=workload_config.dataset_nbytes,
            store_nbytes=path_nbytes(path),
            on_disk_chunks=format_chunks(store_config.on_disk_chunks),
            dask_chunks_arg=format_chunks(dask_config.chunks),
            dask_chunksize=format_chunks(getattr(array, "chunksize", None)),
            workers=dask_config.n_workers,
            threads_per_worker=dask_config.threads_per_worker,
            processes=dask_config.processes,
            workload=workload,
            repeat=repeat,
            task_count=task_count(workload_array),
            elapsed_s=elapsed_s,
            result_shape=result_shape(computed),
            result_nbytes=result_nbytes(computed),
            driver_max_rss_mb=max_rss_mb(),
            worker_rss_before_mb=before_worker_rss,
            worker_rss_after_mb=after_worker_rss,
        )


def iter_store_configs(args: argparse.Namespace) -> Iterator[StoreConfig]:
    for store_type, on_disk_chunks in itertools.product(
        args.store_types, args.on_disk_chunks
    ):
        if store_type == "h5ad":
            yield StoreConfig(
                store_type="h5ad",
                on_disk_chunks=on_disk_chunks,
                zarr_format=args.zarr_format,
                zarr_shards=None,
            )
            continue

        for zarr_shards in args.zarr_shards:
            if args.zarr_format == 2 and zarr_shards is not None:
                msg = "Zarr v2 does not support shards; omit --zarr-shards or use --zarr-format 3"
                raise ValueError(msg)
            yield StoreConfig(
                store_type="zarr",
                on_disk_chunks=on_disk_chunks,
                zarr_format=args.zarr_format,
                zarr_shards=zarr_shards,
            )


def iter_dask_configs(args: argparse.Namespace) -> Iterator[DaskConfig]:
    for chunks, n_workers, threads_per_worker in itertools.product(
        args.dask_chunks, args.workers, args.threads_per_worker
    ):
        yield DaskConfig(
            chunks=chunks,
            n_workers=n_workers,
            threads_per_worker=threads_per_worker,
            processes=args.processes,
        )


def write_results(
    output: Path, rows: Iterable[BenchmarkResult], *, append: bool
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    fieldnames = list(BenchmarkResult.__dataclass_fields__)
    write_header = not append or not output.exists()

    with output.open(mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
            f.flush()


def describe_plan(args: argparse.Namespace) -> dict[str, object]:
    normalize_args(args)
    store_configs = list(iter_store_configs(args))
    dask_configs = list(iter_dask_configs(args))
    return {
        "shape": format_chunks(args.shape),
        "dataset_nbytes": int(
            np.prod(args.shape, dtype=np.int64) * np.dtype("float32").itemsize
        ),
        "store_configs": [asdict(config) for config in store_configs],
        "dask_configs": [asdict(config) for config in dask_configs],
        "workloads": args.workloads,
        "repeats": args.repeats,
        "result_rows": len(store_configs)
        * len(dask_configs)
        * len(args.workloads)
        * args.repeats,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run an exploratory grid over on-disk chunking, Dask chunks, and "
            "Dask worker settings for AnnData lazy array reads."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("benchmarks/results/dask_chunk_grid_stores"),
        help="Directory for generated HDF5/Zarr stores.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/results/dask_chunk_grid.csv"),
        help="CSV path for benchmark rows.",
    )
    parser.add_argument(
        "--shape",
        type=parse_shape,
        default=(12_000, 3_000),
        help="Dense X shape as OBS,VARS.",
    )
    parser.add_argument(
        "--store-types",
        nargs="+",
        choices=("h5ad", "zarr"),
        default=["h5ad", "zarr"],
        help="Storage backends to benchmark.",
    )
    parser.add_argument(
        "--on-disk-chunks",
        action="append",
        type=parse_chunk_spec,
        default=None,
        help="On-disk chunks as OBS,VARS. Repeat for a grid.",
    )
    parser.add_argument(
        "--dask-chunks",
        action="append",
        type=parse_chunk_spec,
        default=None,
        help="Virtual Dask chunks as OBS,VARS, or 'default'. Repeat for a grid.",
    )
    parser.add_argument(
        "--zarr-format",
        type=int,
        choices=(2, 3),
        default=3,
        help="Zarr format for generated Zarr stores.",
    )
    parser.add_argument(
        "--zarr-shards",
        action="append",
        type=parse_chunk_spec,
        default=None,
        help="Zarr v3 shard shape as OBS,VARS. Repeat for a grid.",
    )
    parser.add_argument(
        "--workers",
        action="append",
        type=int,
        default=None,
        help="Dask worker counts to test. Repeat for a grid.",
    )
    parser.add_argument(
        "--threads-per-worker",
        action="append",
        type=int,
        default=None,
        help="Dask threads per worker to test. Repeat for a grid.",
    )
    parser.add_argument(
        "--no-processes",
        action="store_false",
        dest="processes",
        help="Use threaded workers instead of process workers.",
    )
    parser.set_defaults(processes=True)
    parser.add_argument(
        "--workloads",
        nargs="+",
        choices=(
            "sum_axis0",
            "sum_axis1",
            "subset_mean",
            "normalize_log1p_slice",
            "scanpy_normalize_log1p",
        ),
        default=list(DEFAULT_WORKLOADS),
        help="Dask workloads to run.",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--slice-obs", type=int, default=4096)
    parser.add_argument("--slice-vars", type=int, default=1024)
    parser.add_argument("--target-sum", type=float, default=10_000.0)
    parser.add_argument("--memory-limit", default="auto")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--append", action="store_true")
    parser.add_argument("--force", action="store_true", help="Regenerate stores.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the expanded grid without generating data or running Dask.",
    )
    return parser


def normalize_args(args: argparse.Namespace) -> None:
    if args.on_disk_chunks is None:
        args.on_disk_chunks = list(DEFAULT_ON_DISK_CHUNKS)
    if args.dask_chunks is None:
        args.dask_chunks = list(DEFAULT_DASK_CHUNKS)
    if args.zarr_shards is None:
        args.zarr_shards = [None]
    if args.workers is None:
        args.workers = [1, 4]
    if args.threads_per_worker is None:
        args.threads_per_worker = [1]


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    normalize_args(args)

    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if any(worker < 1 for worker in args.workers):
        parser.error("--workers values must be at least 1")
    if any(threads < 1 for threads in args.threads_per_worker):
        parser.error("--threads-per-worker values must be at least 1")

    if args.dry_run:
        print(json.dumps(describe_plan(args), indent=2, default=format_chunks))
        return 0

    args.work_dir.mkdir(parents=True, exist_ok=True)
    store_configs = list(iter_store_configs(args))
    dask_configs = list(iter_dask_configs(args))
    workload_config = WorkloadConfig(
        shape=args.shape,
        slice_obs=args.slice_obs,
        slice_vars=args.slice_vars,
        target_sum=args.target_sum,
        memory_limit=args.memory_limit,
        dataset_nbytes=int(
            np.prod(args.shape, dtype=np.int64) * np.dtype("float32").itemsize
        ),
        run_started_at=datetime.now(UTC).isoformat(timespec="seconds"),
    )

    def rows() -> Iterator[BenchmarkResult]:
        for store_config in store_configs:
            path = store_path(args.work_dir, store_config)
            write_store(path, store_config, args.shape, args.seed, force=args.force)
            for dask_config, workload, repeat in itertools.product(
                dask_configs, args.workloads, range(args.repeats)
            ):
                yield run_one(
                    path,
                    store_config,
                    dask_config,
                    workload=workload,
                    repeat=repeat,
                    workload_config=workload_config,
                )

    write_results(args.output, rows(), append=args.append)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
