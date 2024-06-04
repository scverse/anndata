from __future__ import annotations

from .utils import gen_adata


class GarbargeCollectionSuite:
    def peakmem_garbage_collection(self, *_):
        for _ in range(10):
            adata = gen_adata(10000, 10000, "X-csc")  # noqa: F841
