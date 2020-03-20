from collections.abc import MutableMapping
from typing import Mapping
from warnings import warn

from .._core.access import ElementRef


def _access_warn(key, cur_loc):
    warn(
        f"This location for '{key}' is deprecated. It has been moved to {cur_loc}, "
        "and will not be accesible here in a future version of anndata.",
        FutureWarning,
        stacklevel=3,
    )


class DeprecatedDict(MutableMapping):
    def __init__(self, data, *, deprecated_items: Mapping[str, ElementRef]):
        self.data = data
        self.deprecated_items = deprecated_items

    def __getitem__(self, key):
        if key in self.deprecated_items:
            cur_loc = self.deprecated_items[key]
            _access_warn(key, cur_loc)
            return cur_loc.get()
        else:
            return self.data[key]

    def __setitem__(self, key, value):
        if key in self.deprecated_items:
            cur_loc = self.deprecated_items[key]
            _access_warn(key, cur_loc)
            cur_loc.set(value)
        else:
            self.data[key] = value

    def __delitem__(self, key):
        if key in self.deprecated_items:
            cur_loc = self.deprecated_items[key]
            _access_warn(key, cur_loc)
            cur_loc.delete()
            del self.deprecated_items[key]
        else:
            del self.data[key]

    def __contains__(self, key):
        return key in self.data

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def copy(self) -> dict:
        return self.data.copy()
