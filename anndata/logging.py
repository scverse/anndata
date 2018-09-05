import os
import logging

_previous_memory_usage = -1

anndata_logger = logging.getLogger('anndata')
anndata_logger.propagate = False  # Don’t pass log messages on to logging.root and its handler
anndata_logger.setLevel('INFO')
anndata_logger.addHandler(logging.StreamHandler())  # Logs go to stderr
anndata_logger.handlers[-1].setFormatter(logging.Formatter('%(message)s'))
anndata_logger.handlers[-1].setLevel('INFO')


def get_logger(name):
    """Creates a child logger that delegates to anndata_logger instead to logging.root"""
    return anndata_logger.manager.getLogger(name)


def get_memory_usage():
    import psutil
    process = psutil.Process(os.getpid())
    meminfo_attr = ('memory_info' if hasattr(process, 'memory_info')
                    else 'get_memory_info')
    mem = getattr(process, meminfo_attr)()[0] / 2**30  # output in GB
    mem_diff = mem
    global _previous_memory_usage
    if _previous_memory_usage != -1:
        mem_diff = mem - _previous_memory_usage
    _previous_memory_usage = mem
    return mem, mem_diff


def format_memory_usage(mem_usage, msg='', newline=False):
    mem, diff = mem_usage
    string = (('\n' if newline else '')
              + msg + (' \n... ' if msg != '' else '')
              + 'Memory usage: current {:.2f} GB, difference {:+.2f} GB'
              .format(mem, diff))
    return string


def print_memory_usage(msg='', newline=False):
    string = format_memory_usage(get_memory_usage(), msg, newline)
    print(string)
