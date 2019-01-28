from pathlib import Path

import pytest


@pytest.fixture
def backing_h5ad():
    p = Path('./test.h5ad')
    if p.is_file():
        p.unlink()
    return p
