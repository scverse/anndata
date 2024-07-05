from __future__ import annotations

try:
    from xarray import DataArray
except ImportError:

    class DataArray:
        def __repr__(self) -> str:
            return "mock DataArray"


try:
    import xarray as xr
except ImportError:
    xr = None


try:
    from xarray.backends.zarr import ZarrArrayWrapper
except ImportError:

    class ZarrArrayWrapper:
        def __repr__(self) -> str:
            return "mock ZarrArrayWrapper"


try:
    from xarray.backends import BackendArray
except ImportError:

    class BackendArray:
        def __repr__(self) -> str:
            return "mock BackendArray"


try:
    from xarray import Dataset
except ImportError:

    class Dataset:
        def __repr__(self) -> str:
            return "mock Dataset"
