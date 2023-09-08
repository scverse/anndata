from __future__ import annotations

import warnings

import pytest

import anndata
from anndata.tests.helpers import subset_func  # noqa: F401

# TODO: Should be done in pyproject.toml, see anndata/conftest.py
warnings.filterwarnings("ignore", category=anndata.OldFormatWarning)


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"
