# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a pytest plugin/config that applies to doctests as well
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.
from __future__ import annotations

import re
import warnings
from typing import TYPE_CHECKING

import pytest

from anndata.compat import chdir
from anndata.utils import import_name

if TYPE_CHECKING:
    from pathlib import Path

doctest_marker = pytest.mark.usefixtures("doctest_env")


@pytest.fixture
def doctest_env(
    request: pytest.FixtureRequest, cache: pytest.Cache, tmp_path: Path
) -> None:
    from scanpy import settings

    # request.node.parent is either a DoctestModule or a DoctestTextFile.
    # Only DoctestModule has a .obj attribute (the imported module).
    if request.node.parent.obj:
        func = import_name(request.node.name)
        if warning_detail := getattr(func, "__deprecated", None):
            cat, msg, _ = warning_detail  # type: tuple[type[Warning], str, bool]
            warnings.filterwarnings("ignore", category=cat, message=re.escape(msg))

    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
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
        item.add_marker(doctest_marker)
