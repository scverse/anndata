from __future__ import annotations

from .utils import gen_adata, get_peak_mem, sedate


def run_garbage_collection():
    for _ in range(10):
        adata = gen_adata(10000, 10000, "X-csc")  # noqa: F841


class GarbargeCollectionSuite:
    def track_peakmem_write_compressed(self, *_):
        return get_peak_mem((sedate(run_garbage_collection), (), {}))
