# This file exists
# 1. to allow ignoring warnings without test collection failing on CI
# 2. as a central
# TODO: Fix that, e.g. with the `pytest -p anndata.testing._pytest` pattern.

import pytest


@pytest.fixture(autouse=True)
def doctests_add_tmp_path(doctest_namespace, cache, tmp_path):
    doctest_namespace.update(
        tmp_path=tmp_path,
        data_dir=cache.mkdir("scanpy-data"),
    )
