# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a central
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.

from pathlib import Path
from collections.abc import Iterable
import pytest


@pytest.fixture
def chdir_to_tmp(tmp_path: Path) -> None:
    old = tmp_path.chdir()
    yield
    old.chdir()


def pytest_collection_modifyitems(items: Iterable[pytest.Item]) -> None:
    skip_marker = pytest.mark.usefixtures("chdir_to_tmp")

    for item in items:
        if isinstance(item, pytest.DoctestItem):
            item.add_marker(skip_marker)
