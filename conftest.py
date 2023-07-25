# This file exists just to allow ignoring warnings without test collection failing on CI
# TODO: Fix that
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--only-gpu",
        action="store_true",
        default=False,
        help="Only run GPU tests.",
    )


def pytest_runtest_setup(item):
    gpu = len([mark for mark in item.iter_markers(name="gpu")]) > 0

    if not gpu:  # Skip non gpu tests on gpu CI
        if item.config.getoption("--only-gpu"):
            pytest.skip("Only running GPU tests.")
    else:  # Skip GPU tests if cupy not installed
        from importlib.util import find_spec

        if not find_spec("cupy"):
            pytest.skip("Cupy not installed.")
