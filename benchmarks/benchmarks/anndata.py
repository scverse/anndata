from __future__ import annotations

import tracemalloc

import numpy as np

from .utils import gen_adata


class GarbargeCollectionSuite:
    runs = 10

    # https://github.com/pythonprofilers/memory_profiler/issues/402 and other backend does not pick this up
    def track_peakmem_write_compressed(self, *_):
        def display_top(snapshot, key_type="lineno"):
            snapshot = snapshot.filter_traces(
                (
                    tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
                    tracemalloc.Filter(False, "<unknown>"),
                )
            )
            top_stats = snapshot.statistics(key_type)
            total = sum(stat.size for stat in top_stats)
            return total

        total = np.zeros(self.runs)
        tracemalloc.start()
        for i in range(self.runs):
            data = gen_adata(10000, 10000, "X-csc")  # noqa: F841
            snapshot = tracemalloc.take_snapshot()
            total[i] = display_top(snapshot)
        tracemalloc.stop()
        return max(total)
