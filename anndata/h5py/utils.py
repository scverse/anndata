import numpy as np


def _chunked_rows(X, chunk_size=1000, read=False):
    start = 0
    type = X.dtype
    n = X.shape[0]
    m = X.shape[1]
    data_array = np.empty((chunk_size, m), type) if read else None
    for _ in range(int(n // chunk_size)):
        end = start + chunk_size
        if read:
            X.read_direct(data_array, source_sel=np.s_[start:end, :])
        else:
            data_array = X[start:end]
        yield (data_array, start, end)
        start = end
    if start < n:
        if read:
            data_array = np.empty((n-start, m), type)
            X.read_direct(data_array, source_sel=np.s_[start:n, :])
        else:
            data_array = X[start:n]
        yield (data_array, start, n)
