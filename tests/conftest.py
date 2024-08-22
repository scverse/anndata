from __future__ import annotations

import pytest

from anndata.tests.helpers import subset_func  # noqa: F401


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"
