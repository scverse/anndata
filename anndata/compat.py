import numpy as np
import pandas as pd


# try importing zarr, dask, and zappy
try:
    from zarr.core import Array as ZarrArray
except ImportError:

    class ZarrArray:
        @staticmethod
        def __repr__():
            return 'mock zarr.core.Array'


try:
    from zappy.base import ZappyArray
except ImportError:

    class ZappyArray:
        @staticmethod
        def __repr__():
            return 'mock zappy.base.ZappyArray'


try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        @staticmethod
        def __repr__():
            return 'mock dask.array.core.Array'


def version(package):
    try:
        from importlib.metadata import version
    except ImportError:
        from importlib_metadata import version
    return version(package)


def _from_fixed_length_strings(value):
    """Convert from fixed length strings to unicode.

    For backwards compatability with older h5ad and zarr files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        if isinstance(
            dt_type, tuple
        ):  # vlen strings, could probably match better
            # Fixing issues with h5py v2.10.0, could be cleaner
            if issubclass(np.dtype(dt_type[0]).type, np.string_):
                dt_list[1] = 'U{}'.format(int(dt_type[0][2:]))
                new_dtype.append(tuple(dt_list))
            else:
                dt_list[1] = "O"
                new_dtype.append(tuple(dt_list))
        elif issubclass(np.dtype(dt_type).type, np.string_):
            dt_list[1] = f'U{int(dt_type[2:])}'
            new_dtype.append(tuple(dt_list))
        else:
            new_dtype.append(dt)
    return value.astype(new_dtype)


def _clean_uns(d: dict):
    """Compat function for when categorical keys were stored in uns."""
    k_to_delete = []
    for k, v in d.get("uns", {}).items():
        if k.endswith('_categories'):
            k_stripped = k.replace('_categories', '')
            if isinstance(
                v, (str, int)
            ):  # fix categories with a single category
                v = [v]
            for ann in ['obs', 'var']:
                if k_stripped in d[ann]:
                    d[ann][k_stripped] = pd.Categorical.from_codes(
                        codes=d[ann][k_stripped].values, categories=v
                    )
                    k_to_delete.append(k)
    for k in k_to_delete:
        del d["uns"][k]
