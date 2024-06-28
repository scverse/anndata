from __future__ import annotations

import logging
import os

_previous_memory_usage = None

anndata_logger = logging.getLogger("anndata")
# Don’t pass log messages on to logging.root and its handler
anndata_logger.propagate = False
anndata_logger.addHandler(logging.StreamHandler())  # Logs go to stderr
anndata_logger.handlers[-1].setFormatter(logging.Formatter("%(message)s"))
anndata_logger.handlers[-1].setLevel("INFO")


def get_logger(name):
    """\
    Creates a child logger that delegates to anndata_logger
    instead to logging.root
    """
    return anndata_logger.manager.getLogger(name)


def get_memory_usage():
    import psutil

    process = psutil.Process(os.getpid())
    try:
        meminfo = process.memory_info()
    except AttributeError:
        meminfo = process.get_memory_info()
    mem = meminfo[0] / 2**30  # output in GB
    mem_diff = mem
    global _previous_memory_usage
    if _previous_memory_usage is not None:
        mem_diff = mem - _previous_memory_usage
    _previous_memory_usage = mem
    return mem, mem_diff


def format_memory_usage(mem_usage, msg="", newline=False):
    newline = "\n" if newline else ""
    more = " \n... " if msg != "" else ""
    mem, diff = mem_usage
    return (
        f"{newline}{msg}{more}"
        f"Memory usage: current {mem:.2f} GB, difference {diff:+.2f} GB"
    )


def print_memory_usage(msg="", newline=False):
    print(format_memory_usage(get_memory_usage(), msg, newline))
