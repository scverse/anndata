"""Tests for the IO-side thread-pool manager.

:class:`anndata._io._pool.PoolManager` is exercised end-to-end through
the HDF5 parallel-read path (see :mod:`test_io_h5_parallel`); these
tests pin its standalone semantics:

* lazy creation on first ``distribute_tasks_to_threads``,
* reuse when the requested worker count is unchanged,
* atomic resize-and-submit when the worker count changes,
* ``reset()`` returns the manager to its pre-use state,
* the configured thread-name prefix is applied to worker threads.
"""

from __future__ import annotations

import threading
from concurrent.futures import wait
from typing import TYPE_CHECKING

import pytest

from anndata._io._pool import PoolManager

if TYPE_CHECKING:
    from collections.abc import Generator


@pytest.fixture
def mgr() -> Generator[PoolManager, None, None]:
    """A fresh manager per test — prevents cross-test pool leakage."""
    pm = PoolManager(thread_name_prefix="anndata-test-pool")
    yield pm
    pm.reset()


def _square(x: int) -> int:
    return x * x


def test_returns_one_future_per_task(mgr: PoolManager) -> None:
    """One Future is returned per input task, in input order, each
    resolving to ``fn(*task)``."""
    futs = mgr.distribute_tasks_to_threads(4, _square, [(i,) for i in range(8)])
    assert [f.result() for f in futs] == [i * i for i in range(8)]


def test_empty_tasks_returns_empty_list(mgr: PoolManager) -> None:
    """An empty tasks iterable produces no futures and does not raise.

    Importantly, the pool *is* created on first call (even with no
    work) so that the manager's invariant (`_pool is not None` after
    any ``distribute_tasks_to_threads``) holds uniformly — exercised
    by the reuse test below.
    """
    futs = mgr.distribute_tasks_to_threads(4, _square, [])
    assert futs == []
    assert mgr._pool is not None
    assert mgr._n_workers == 4


def test_reuses_pool_when_workers_unchanged(mgr: PoolManager) -> None:
    """Repeat calls with the same worker count must return the same
    underlying executor — this is the persistent-pool win the manager
    exists for."""
    mgr.distribute_tasks_to_threads(4, _square, [(1,)])
    pool_a = mgr._pool
    mgr.distribute_tasks_to_threads(4, _square, [(2,)])
    pool_b = mgr._pool
    assert pool_a is pool_b


def test_resizes_when_workers_change(mgr: PoolManager) -> None:
    """Different worker count triggers a rebuild — the old pool is
    discarded and a new one of the requested size takes its place."""
    mgr.distribute_tasks_to_threads(2, _square, [(1,)])
    pool_a = mgr._pool
    assert mgr._n_workers == 2
    mgr.distribute_tasks_to_threads(5, _square, [(2,)])
    pool_b = mgr._pool
    assert pool_b is not pool_a
    assert mgr._n_workers == 5


def test_resize_still_completes_in_flight_work(mgr: PoolManager) -> None:
    """After a resize, futures submitted to the *old* pool must still
    complete. ``shutdown(wait=False)`` only blocks new submissions —
    work already enqueued runs to completion."""
    # Submit a batch that takes a moment, so we can race a resize past it.
    sentinel = threading.Event()

    def slow(x: int) -> int:
        sentinel.wait(timeout=5.0)
        return x + 1

    old_futs = mgr.distribute_tasks_to_threads(2, slow, [(10,), (20,)])
    # Resize while the old tasks are parked on the event.
    mgr.distribute_tasks_to_threads(4, _square, [(3,)])
    sentinel.set()
    # Both groups of tasks must complete.
    wait(old_futs)
    assert sorted(f.result() for f in old_futs) == [11, 21]


def test_reset_tears_down_pool(mgr: PoolManager) -> None:
    """``reset()`` waits for in-flight tasks and nulls the cached pool;
    subsequent ``distribute_tasks_to_threads`` lazily creates a fresh
    one."""
    mgr.distribute_tasks_to_threads(2, _square, [(1,)])
    assert mgr._pool is not None
    mgr.reset()
    assert mgr._pool is None
    assert mgr._n_workers == 0

    futs = mgr.distribute_tasks_to_threads(3, _square, [(4,)])
    assert mgr._pool is not None
    assert mgr._n_workers == 3
    assert futs[0].result() == 16


def test_reset_on_unused_manager_is_noop(mgr: PoolManager) -> None:
    """``reset()`` before any submission must not raise — supports a
    uniform test-teardown idiom regardless of whether the test
    exercised the manager."""
    mgr.reset()
    assert mgr._pool is None


def test_thread_name_prefix_is_applied() -> None:
    """The constructor's prefix shows up on worker thread names so
    profilers and crash dumps can attribute the threads back to anndata."""
    mgr = PoolManager(thread_name_prefix="anndata-test-prefix")
    try:
        name_box: list[str] = []
        mgr.distribute_tasks_to_threads(
            1, lambda: name_box.append(threading.current_thread().name), [()]
        )[0].result()
        assert name_box[0].startswith("anndata-test-prefix")
    finally:
        mgr.reset()


def test_concurrent_resize_and_submit_is_race_free() -> None:
    """The central guarantee: two threads alternating
    ``distribute_tasks_to_threads`` with *different* worker counts
    must never trip a ``RuntimeError: cannot schedule new futures
    after shutdown``.

    Pre-refactor, ``_get_pool`` released the lock before the caller
    could submit; a concurrent resize could shut down the pool out
    from under an in-flight submission.
    ``distribute_tasks_to_threads`` holds the lock through submission,
    so the race cannot occur.
    """
    mgr = PoolManager(thread_name_prefix="anndata-test-race")
    try:
        errors: list[BaseException] = []
        # Many alternating calls, each toggling between two sizes.
        # If the race exists, on a contended scheduler at least one
        # submission would land on a freshly-shutdown pool.
        N = 200

        def loop(num_workers: int) -> None:
            for _ in range(N):
                try:
                    futs = mgr.distribute_tasks_to_threads(
                        num_workers, _square, [(1,), (2,), (3,)]
                    )
                    wait(futs)
                    for f in futs:
                        f.result()
                except BaseException as e:  # noqa: BLE001
                    errors.append(e)

        t1 = threading.Thread(target=loop, args=(2,))
        t2 = threading.Thread(target=loop, args=(8,))
        t1.start()
        t2.start()
        t1.join()
        t2.join()
        assert not errors, f"unexpected errors during race test: {errors[:3]}"
    finally:
        mgr.reset()
