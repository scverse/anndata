import pytest


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"
