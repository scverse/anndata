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
    from collections.abc import Generator, Iterable
    from pathlib import Path


doctest_marker = pytest.mark.usefixtures("doctest_env")


@pytest.fixture
def doctest_env(
    request: pytest.FixtureRequest, cache: pytest.Cache, tmp_path: Path
) -> Generator[None, None, None]:
    from scanpy import settings

    assert isinstance(request.node.parent, pytest.Module)
    # request.node.parent is either a DoctestModule or a DoctestTextFile.
    # Only DoctestModule has a .obj attribute (the imported module).
    if request.node.parent.obj:
        func = import_name(request.node.name)
        warning_detail: tuple[type[Warning], str, bool] | None
        if warning_detail := getattr(func, "__deprecated", None):
            cat, msg, _ = warning_detail
            warnings.filterwarnings("ignore", category=cat, message=re.escape(msg))

    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
    settings.datasetdir = old_dd


def pytest_itemcollected(item: pytest.Item) -> None:
    """Define behavior of pytest.mark.gpu and doctests."""
    from importlib.util import find_spec

    is_gpu = len([mark for mark in item.iter_markers(name="gpu")]) > 0
    if is_gpu:
        item.add_marker(
            pytest.mark.skipif(not find_spec("cupy"), reason="Cupy not installed.")
        )

    if isinstance(item, pytest.DoctestItem):
        item.add_marker(doctest_marker)


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--strict-warnings",
        action="store_true",
        default=False,
        help="Turn most warnings into errors",
    )


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: Iterable[pytest.Item]
):
    if not config.getoption("--strict-warnings"):
        return

    filters_from_config: list[str] = []
    if fs := config.getini("filterwarnings"):
        assert isinstance(fs, list)
        filters_from_config = [str(f) for f in fs]

    warning_filters = [
        r"error",
        *filters_from_config,
        r"default::anndata._warnings.ImplicitModificationWarning",
        r"default:Transforming to str index:UserWarning",
        r"default:(Observation|Variable) names are not unique. To make them unique:UserWarning",
        r"default::scipy.sparse.SparseEfficiencyWarning",
        r"default::dask.array.core.PerformanceWarning",
    ]
    warning_marks = [pytest.mark.filterwarnings(f) for f in warning_filters]

    for item in items:
        # reverse and then prepend means essentially `marks[0:0] = warning_marks`
        # this ensures that markers that are applied later override these
        for mark in reversed(warning_marks):
            item.add_marker(mark, append=False)
