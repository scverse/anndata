# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a pytest plugin/config that applies to doctests as well
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.

from pathlib import Path
from collections.abc import Iterable
import pytest


doctest_skip_marker = pytest.mark.usefixtures("doctest_env")


@pytest.fixture
def doctest_env(cache: pytest.Cache, tmp_path: Path) -> None:
    from scanpy import settings

    old_wd = tmp_path.chdir()
    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    yield
    old_wd.chdir()
    settings.datasetdir = old_dd


def pytest_itemcollected(item):
    """Define behavior of pytest.mark.gpu and doctests."""
    from importlib.util import find_spec

    is_gpu = len([mark for mark in item.iter_markers(name="gpu")]) > 0
    if is_gpu:
        item.add_marker(
            pytest.mark.skipif(not find_spec("cupy"), reason="Cupy not installed.")
        )

    if isinstance(item, pytest.DoctestItem):
        item.add_marker(doctest_skip_marker)
