import sys
from pathlib import Path

import pytest


if sys.version_info < (3, 6):
    from _pytest.tmpdir import _mk_tmp  # noqa

    # On 3.5 this is pathlib2, which won’t work.
    @pytest.fixture
    def tmp_path(request, tmp_path_factory):
        return Path(str(_mk_tmp(request, tmp_path_factory)))


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / 'test.h5ad'
