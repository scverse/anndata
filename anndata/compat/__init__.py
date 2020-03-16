from typing import Union, Mapping, MutableMapping

import numpy as np
import pandas as pd


# try importing zarr, dask, and zappy
from packaging import version

try:
    from zarr.core import Array as ZarrArray
except ImportError:

    class ZarrArray:
        @staticmethod
        def __repr__():
            return "mock zarr.core.Array"


try:
    from zappy.base import ZappyArray
except ImportError:

    class ZappyArray:
        @staticmethod
        def __repr__():
            return "mock zappy.base.ZappyArray"


try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        @staticmethod
        def __repr__():
            return "mock dask.array.core.Array"


try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


def pkg_version(package):
    try:
        from importlib.metadata import version as v
    except ImportError:
        from importlib_metadata import version as v
    return version.parse(v(package))


def _from_fixed_length_strings(value):
    """\
    Convert from fixed length strings to unicode.

    For backwards compatability with older h5ad and zarr files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        # could probably match better
        is_annotated = isinstance(dt_type, tuple)
        if is_annotated:
            dt_type = dt_type[0]
        # Fixing issue introduced with h5py v2.10.0, see:
        # https://github.com/h5py/h5py/issues/1307
        if issubclass(np.dtype(dt_type).type, np.string_):
            dt_list[1] = f"U{int(dt_type[2:])}"
        elif is_annotated or np.issubdtype(np.dtype(dt_type), np.str_):
            dt_list[1] = "O"  # Assumption that it’s a vlen str
        new_dtype.append(tuple(dt_list))
    return value.astype(new_dtype)


def _to_fixed_length_strings(value: np.ndarray) -> np.ndarray:
    """\
    Convert variable length strings to fixed length.

    Currently a workaround for
    https://github.com/zarr-developers/zarr-python/pull/422
    """
    new_dtype = []
    for dt_name, (dt_type, dt_offset) in value.dtype.fields.items():
        if dt_type.kind == "O":
            #  Assuming the objects are str
            size = max(len(x.encode()) for x in value.getfield("O", dt_offset))
            new_dtype.append((dt_name, ("U", size)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


def _clean_uns(d: Mapping[str, MutableMapping[str, Union[pd.Series, str, int]]]):
    """
    Compat function for when categorical keys were stored in uns.
    This used to be buggy because when storing categorical columns in obs and var with
    the same column name, only one `<colname>_categories` is retained.
    """
    k_to_delete = set()
    for cats_name, cats in d.get("uns", {}).items():
        if not cats_name.endswith("_categories"):
            continue
        name = cats_name.replace("_categories", "")
        # fix categories with a single category
        if isinstance(cats, (str, int)):
            cats = [cats]
        for ann in ["obs", "var"]:
            if name not in d[ann]:
                continue
            codes: np.ndarray = d[ann][name].values
            # hack to maybe find the axis the categories were for
            if not np.all(codes < len(cats)):
                continue
            d[ann][name] = pd.Categorical.from_codes(codes, cats)
            k_to_delete.add(cats_name)
    for cats_name in k_to_delete:
        del d["uns"][cats_name]
