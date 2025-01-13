from __future__ import annotations

import logging
import os

from .compat import old_positionals

_previous_memory_usage = None

anndata_logger = logging.getLogger("anndata")
# Don’t pass log messages on to logging.root and its handler
anndata_logger.propagate = False
anndata_logger.addHandler(logging.StreamHandler())  # Logs go to stderr
anndata_logger.handlers[-1].setFormatter(logging.Formatter("%(message)s"))
anndata_logger.handlers[-1].setLevel("INFO")


def get_logger(name: str) -> logging.Logger:
    """\
    Creates a child logger that delegates to anndata_logger
    instead to logging.root
    """
    return anndata_logger.manager.getLogger(name)


def get_memory_usage() -> tuple[float, float]:
    import psutil

    process = psutil.Process(os.getpid())
    try:
        meminfo = process.memory_info()
    except AttributeError:
        meminfo = process.get_memory_info()
    mem = meminfo[0] / 2**30  # output in GB
    mem_diff = mem
    global _previous_memory_usage  # noqa: PLW0603
    if _previous_memory_usage is not None:
        mem_diff = mem - _previous_memory_usage
    _previous_memory_usage = mem
    return mem, mem_diff


@old_positionals("newline")
def format_memory_usage(
    mem_usage: tuple[float, float], msg: str = "", *, newline: bool = False
):
    nl = "\n" if newline else ""
    more = " \n... " if msg != "" else ""
    mem, diff = mem_usage
    return (
        f"{nl}{msg}{more}Memory usage: current {mem:.2f} GB, difference {diff:+.2f} GB"
    )


@old_positionals("newline")
def print_memory_usage(msg: str = "", *, newline: bool = False):
    print(format_memory_usage(get_memory_usage(), msg, newline))
