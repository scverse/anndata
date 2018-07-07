def _chunked_rows(X, chunk_size=1000):
    start = 0
    n = X.shape[0]
    for _ in range(int(n // chunk_size)):
        end = start + chunk_size
        yield (X[start:end], start, end)
        start = end
    if start < n:
        yield (X[start:n], start, n)
