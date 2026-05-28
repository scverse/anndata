"""Generic, process-global thread-pool management for IO-side parallelism."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from concurrent.futures import Future
    from typing import TypeVar

    R = TypeVar("R")


class PoolManager:
    def __init__(self, *, thread_name_prefix: str) -> None:
        self._lock = threading.Lock()
        self._pool: ThreadPoolExecutor | None = None
        self._n_workers: int = 0
        self._thread_name_prefix = thread_name_prefix

    def distribute_tasks_to_threads(
        self,
        n_workers: int,
        fn: Callable[..., R],
        tasks: Sequence[tuple],
    ) -> list[Future[R]]:
        """Submit one ``fn(*t)`` call per ``t`` of ``tasks`` to a
        pool of ``n_workers`` threads, returning the resulting Futures in
        input order.
        """
        with self._lock:
            self._ensure_pool_locked(n_workers)
            assert self._pool is not None  # for type checkers
            return [self._pool.submit(fn, *task) for task in tasks]

    def reset(self) -> None:
        with self._lock:
            if self._pool is not None:
                self._pool.shutdown(wait=True)
                self._pool = None
                self._n_workers = 0

    def _ensure_pool_locked(self, n_workers: int) -> None:
        if self._pool is not None and self._n_workers == n_workers:
            return
        if self._pool is not None:
            self._pool.shutdown(wait=False, cancel_futures=False)
        self._pool = ThreadPoolExecutor(
            max_workers=n_workers, thread_name_prefix=self._thread_name_prefix
        )
        self._n_workers = n_workers
