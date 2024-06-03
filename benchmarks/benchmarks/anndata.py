from __future__ import annotations

from .utils import gen_adata


class GarbargeCollectionSuite:
    def track_peakmem_garbage_collection(self, *_):
        for _ in range(50):
            adata = gen_adata(10000, 10000, "X-csc")  # noqa: F841
