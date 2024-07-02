from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

xfail_if_numpy2_loompy = pytest.mark.xfail(
    Version(np.__version__) >= Version("2a0"),
    reason="loompy still uses `np.string_`, removed in `numpy==2.0`",
    raises=AttributeError,
)
