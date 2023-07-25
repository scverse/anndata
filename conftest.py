# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a pytest plugin/config that applies to doctests as well
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.

from pathlib import Path
from collections.abc import Iterable
import pytest


@pytest.fixture
def doctest_env(cache: pytest.Cache, tmp_path: Path) -> None:
    from scanpy import settings

    old_wd = tmp_path.chdir()
    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    yield
    old_wd.chdir()
    settings.datasetdir = old_dd


def pytest_collection_modifyitems(items: Iterable[pytest.Item]) -> None:
    skip_marker = pytest.mark.usefixtures("doctest_env")

    for item in items:
        if isinstance(item, pytest.DoctestItem):
            item.add_marker(skip_marker)
