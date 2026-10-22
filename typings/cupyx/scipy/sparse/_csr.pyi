# SPDX-License-Identifier: MPL-2.0
from typing import Literal

import cupy.cuda
import scipy.sparse as sps

from ._compressed import _compressed_sparse_matrix

class csr_matrix(_compressed_sparse_matrix):
    format: Literal["csr"] = "csr"
    def copy(self) -> csr_matrix: ...
    def get(self, stream: cupy.cuda.Stream | None = None) -> sps.csr_matrix: ...
