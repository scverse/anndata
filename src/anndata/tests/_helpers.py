from __future__ import annotations

import os

import pytest

xfail_if_numpy2_loompy = pytest.mark.xfail(
    os.environ.get("DEPENDENCIES_VERSION", "latest") == "pre-release",
    reason="loompy still uses `np.string_`, removed in `numpy==2.0`",
    raises=AttributeError,
)
