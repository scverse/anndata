from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pytest


def load_script_module():
    script = Path(__file__).parents[1] / "benchmarks" / "scripts" / "dask_chunk_grid.py"
    spec = importlib.util.spec_from_file_location("dask_chunk_grid", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_parse_chunk_spec_accepts_default_and_axis_chunks():
    module = load_script_module()

    assert module.parse_chunk_spec("default") is None
    assert module.parse_chunk_spec("1024,-1") == (1024, -1)
    assert module.parse_chunk_spec("512xnone") == (512, None)


def test_parse_chunk_spec_rejects_wrong_rank():
    module = load_script_module()

    with pytest.raises(argparse.ArgumentTypeError):
        module.parse_chunk_spec("1024")


def test_materialize_chunks_clamps_to_shape_and_expands_full_axis():
    module = load_script_module()

    assert module.materialize_chunks((2048, -1), (1000, 300)) == (1000, 300)
    assert module.materialize_chunks((None, 128), (1000, 300)) == (1000, 128)


def test_path_nbytes_counts_files_and_directories(tmp_path):
    module = load_script_module()
    file_path = tmp_path / "one.bin"
    dir_path = tmp_path / "nested"
    file_path.write_bytes(b"1234")
    dir_path.mkdir()
    (dir_path / "two.bin").write_bytes(b"12")
    (dir_path / "three.bin").write_bytes(b"123")

    assert module.path_nbytes(file_path) == 4
    assert module.path_nbytes(dir_path) == 5


def test_describe_plan_counts_expanded_grid():
    module = load_script_module()
    parser = module.build_parser()
    args = parser.parse_args([
        "--shape",
        "100,20",
        "--store-types",
        "h5ad",
        "--on-disk-chunks",
        "25,20",
        "--dask-chunks",
        "50,-1",
        "--workers",
        "2",
        "--workloads",
        "sum_axis0",
        "--repeats",
        "2",
    ])

    plan = module.describe_plan(args)

    assert plan["result_rows"] == 2
    assert plan["shape"] == "100x20"
    assert plan["dataset_nbytes"] == 8000
