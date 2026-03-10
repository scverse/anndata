"""Private anndata pytest plugin.

This file exists
1. to allow ignoring warnings without test collection failing on CI
2. as a pytest plugin/config that applies to doctests as well

It lives outside of the anndata package in order to avoid importing anndata too early.
"""

from __future__ import annotations

import re
import warnings
from importlib.metadata import version
from importlib.util import find_spec
from typing import TYPE_CHECKING, cast

import pandas as pd
import pytest
from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence
    from pathlib import Path

    from ._doctest import WarningFilter

# Use a marker present in the environment so VS Code’s tests behave identical
IS_PRE = Version(version("zarr")).is_prerelease

# Hack, but I didn’t feel like adding rST syntax to define warning filters
# TODO: remove filters (here and elsewhere) once https://github.com/scverse/scanpy/issues/3879 is fixed
_RST_FILTERS: Sequence[WarningFilter] = (
    ("ignore", r"Moving element.*uns.*to.*obsp", FutureWarning, "", 0),
)


def setup_env() -> None:
    import anndata

    anndata.settings.reset(anndata.settings._registered_options.keys())

    if IS_PRE:
        # https://pandas.pydata.org/docs/whatsnew/v2.3.0.html#upcoming-changes-in-pandas-3-0
        pd.options.future.infer_string = True


@pytest.fixture(scope="session", autouse=True)
def _anndata_session_env(request: pytest.FixtureRequest) -> None:
    setup_env()


@pytest.fixture(autouse=True)
def _anndata_test_env(request: pytest.FixtureRequest) -> None:
    if isinstance(request.node, pytest.DoctestItem):
        request.getfixturevalue("_doctest_env")

    setup_env()


@pytest.fixture
def _doctest_env(
    request: pytest.FixtureRequest, cache: pytest.Cache, tmp_path: Path
) -> Generator[None, None, None]:
    from contextlib import chdir

    from scanpy import settings

    from anndata.utils import import_name

    assert isinstance(request.node, pytest.DoctestItem)
    assert isinstance(request.node.parent, pytest.Module)
    # request.node.parent is either a DoctestModule or a DoctestTextFile.
    # Only DoctestModule has a .obj attribute (the imported module).
    if request.node.parent.obj:
        func = import_name(request.node.name)
        if msg := cast("str | None", getattr(func, "__deprecated__", None)):
            warnings.filterwarnings(
                "ignore", category=FutureWarning, message=re.escape(msg)
            )
        if (
            mod := cast("str | None", getattr(func, "_doctest_needs", None))
        ) is not None and not find_spec(mod):
            request.applymarker(pytest.skip(reason=f"doctest needs {mod} to run"))
        for filter in cast(
            "Sequence[WarningFilter]", getattr(func, "_doctest_warning_filter", ())
        ):
            warnings.filterwarnings(*filter)
    elif request.node.name.endswith(".rst"):
        for filter in _RST_FILTERS:
            warnings.filterwarnings(*filter)

    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
    settings.datasetdir = old_dd


if find_spec("jax"):
    import jax

    jax.config.update("jax_enable_x64", True)  # noqa: FBT003


def pytest_itemcollected(item: pytest.Item) -> None:
    """Define behavior of pytest.mark.{gpu,array_api}."""
    for mark, package in [("gpu", "cupy"), ("array_api", "jax")]:
        is_marked = len(list(item.iter_markers(name=mark))) > 0
        if is_marked:
            item.add_marker(
                pytest.mark.skipif(
                    not find_spec(package), reason=f"{package} not installed."
                )
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
    for item in items:
        if "zarr" in item.name:
            item.add_marker("zarr_io")

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
        return cast("list[str]", strs)
    return []
