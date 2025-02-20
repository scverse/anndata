from __future__ import annotations

import tracemalloc

import numpy as np

from .utils import gen_adata


class GarbargeCollectionSuite:
    runs = 10

    # custom because `memory_profiler` is a line-by-line profiler (also: https://github.com/pythonprofilers/memory_profiler/issues/402)
    def track_peakmem_garbage_collection(self, *_):
        def display_top(snapshot, key_type="lineno"):
            snapshot = snapshot.filter_traces(
                (
                    tracemalloc.Filter(
                        inclusive=False,
                        filename_pattern="<frozen importlib._bootstrap>",
                    ),
                    tracemalloc.Filter(
                        inclusive=False,
                        filename_pattern="<unknown>",
                    ),
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
