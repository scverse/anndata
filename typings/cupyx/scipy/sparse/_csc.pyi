# SPDX-License-Identifier: MPL-2.0
from typing import Literal

import cupy.cuda
import scipy.sparse as sps

from ._compressed import _compressed_sparse_matrix

class csc_matrix(_compressed_sparse_matrix):
    format: Literal["csc"] = "csc"
    def copy(self) -> csc_matrix: ...
    def get(self, stream: cupy.cuda.Stream | None = None) -> sps.csc_matrix: ...
