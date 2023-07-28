# This file exists just to allow ignoring warnings without test collection failing on CI
# TODO: Fix that
import pytest


def pytest_itemcollected(item):
    # Defniing behavior of pytest.mark.gpu
    from importlib.util import find_spec

    gpu = len([mark for mark in item.iter_markers(name="gpu")]) > 0

    if gpu:
        item.add_marker(
            pytest.mark.skipif(not find_spec("cupy"), reason="Cupy not installed.")
        )
