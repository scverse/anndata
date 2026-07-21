from __future__ import annotations

import numpy as np

from anndata._core.index import _process_index_for_h5py


class ProcessIndexForH5py:
    param_names = ("size", "scenario")
    params = (
        [10_000, 100_000, 1_000_000],
        ["sorted_unique", "unsorted_unique", "duplicate_heavy"],
    )

    def setup(self, size, scenario):
        rng = np.random.default_rng(0)
        if scenario == "sorted_unique":
            self.idx = np.arange(size, dtype=np.int64)
        elif scenario == "unsorted_unique":
            self.idx = rng.permutation(size).astype(np.int64)
        else:
            self.idx = rng.integers(0, max(1, size // 10), size=size, dtype=np.int64)

    def time_process_index_for_h5py(self, size, scenario):
        _process_index_for_h5py(self.idx)
