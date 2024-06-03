# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a pytest plugin/config that applies to doctests as well
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.
from __future__ import annotations

import re
import warnings
from typing import TYPE_CHECKING, cast

import pytest

from anndata.compat import chdir
from anndata.utils import import_name

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable
    from pathlib import Path


# TODO: Should be done in pyproject.toml eventually
# See https://github.com/pytest-dev/pytest-cov/issues/437
def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "filterwarnings", "ignore::anndata._warnings.OldFormatWarning"
    )
    config.addinivalue_line(
        "filterwarnings", "ignore::anndata._warnings.ExperimentalFeatureWarning"
    )


@pytest.fixture(autouse=True)
def _suppress_env_for_doctests(request: pytest.FixtureRequest) -> None:
    if isinstance(request.node, pytest.DoctestItem):
        request.getfixturevalue("_doctest_env")


@pytest.fixture()
def _doctest_env(
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
    """Define behavior of pytest.mark.gpu."""
    from importlib.util import find_spec

    is_gpu = len([mark for mark in item.iter_markers(name="gpu")]) > 0
    if is_gpu:
        item.add_marker(
            pytest.mark.skipif(not find_spec("cupy"), reason="Cupy not installed.")
        )


def pytest_addoption(parser: pytest.Parser) -> None:
    """Hook to register custom CLI options and config values"""
    parser.addoption(
        "--strict-warnings",
        action="store_true",
        default=False,
        help="Turn warnings into errors that are not overridden by `filterwarnings` or `filterwarnings_when_strict`.",
    )

    parser.addini(
        "filterwarnings_when_strict",
        "Filters to apply after `-Werror` when --strict-warnings is active",
        type="linelist",
        default=[],
    )


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: Iterable[pytest.Item]
):
    if not config.getoption("--strict-warnings"):
        return

    warning_filters = [
        "error",
        *_config_get_strlist(config, "filterwarnings"),
        *_config_get_strlist(config, "filterwarnings_when_strict"),
    ]
    warning_marks = [pytest.mark.filterwarnings(f) for f in warning_filters]

    # Add warning filters defined in the config to all tests items.
    # Test items might already have @pytest.mark.filterwarnings applied,
    # so we prepend ours to ensure that an item’s explicit filters override these.
    # Reversing then individually prepending ensures that the order is preserved.
    for item in items:
        for mark in reversed(warning_marks):
            item.add_marker(mark, append=False)


def _config_get_strlist(config: pytest.Config, name: str) -> list[str]:
    if strs := config.getini(name):
        assert isinstance(strs, list)
        assert all(isinstance(item, str) for item in strs)
        return cast(list[str], strs)
    return []
