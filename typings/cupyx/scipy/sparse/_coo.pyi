# SPDX-License-Identifier: MPL-2.0
from typing import Literal

import cupy.cuda
import scipy.sparse as sps

from ._base import spmatrix

class coo_matrix(spmatrix):
    format: Literal["coo"] = "coo"
    def get(self, stream: cupy.cuda.Stream | None = None) -> sps.coo_matrix: ...
