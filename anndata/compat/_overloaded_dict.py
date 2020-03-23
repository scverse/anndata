from collections.abc import MutableMapping
from functools import partial
from typing import Any, Callable, Mapping, Optional
from warnings import warn


def _access_warn(key, cur_loc):
    warn(
        f"This location for '{key}' is deprecated. It has been moved to {cur_loc}, "
        "and will not be accesible here in a future version of anndata.",
        FutureWarning,
        stacklevel=4,
    )


class KeyOverload:
    """
    This class contains the information neccesary to overload a key of a dict.

    It's like a descriptor, but for keys of a dict.
    """

    def __init__(
        self,
        key,
        get: Optional[Callable] = None,
        set: Optional[Callable] = None,
        delete: Optional[Callable] = None,
    ):
        self.key = key
        if get is not None:
            self._get = get
        if set is not None:
            self._set = set
        if delete is not None:
            self._delete = delete

    @staticmethod
    def _get(parent, key):
        """Default key getter."""
        return parent.data[key]

    @staticmethod
    def _set(parent, key, value):
        parent.data[key] = value

    @staticmethod
    def _delete(parent, key):
        del parent.data[key]

    @property
    def get(self):
        return partial(self._get, self.parent, self.key)

    @property
    def set(self):
        return partial(self._set, self.parent, self.key)

    @property
    def delete(self):
        return partial(self._delete, self.parent, self.key)


class OverloadedDict(MutableMapping):
    """A dict where some of the keys have been overloaded.

    These keys don't show up in iteration and are not copied.
    """

    def __init__(self, data, *, overloaded_keys: Mapping[Any, KeyOverload]):
        self.data = data
        self.overloaded_keys = overloaded_keys
        for v in overloaded_keys.values():
            v.parent = self

    def __getitem__(self, key):
        if key in self.overloaded_keys:
            return self.overloaded_keys[key].get()
        else:
            return self.data[key]

    def __setitem__(self, key, value):
        if key in self.overloaded_keys:
            self.overloaded_keys[key].set(value)
        else:
            self.data[key] = value

    def __delitem__(self, key):
        if key in self.overloaded_keys:
            self.overloaded_keys[key].delete()
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


def _overloaded_uns(adata):
    def _adjacency_getter(ovld, key, adata: "AnnData"):
        """For overloading:

        >>> mtx = adata.uns["neighbors"]["connectivities"]
        >>> mtx = adata.uns["neighbors"]["distances"]
        """
        _access_warn(key, f".obsp[{key}]")
        return adata.obsp[key]

    def _adjacency_setter(ovld, key, value, adata: "AnnData"):
        """For overloading:

        >>> adata.uns["neighbors"]["connectivities"] = mtx
        >>> adata.uns["neighbors"]["distances"] = mtx
        """
        _access_warn(key, f".obsp[{key}]")
        adata.obsp[key] = value

    def _neighbors_setter(ovld, key, neighbors: Mapping, adata: "AnnData"):
        """For overloading: `adata.uns["neighbors"] = d`."""
        for k in ("distances", "connectivities"):
            if k in neighbors:
                _access_warn(k, f".obsp[{k}]")
                adata.obsp[k] = neighbors.pop(k)
        ovld.data[key] = neighbors

    def _neighbors_getter(ovld, key, adata: "AnnData"):
        """For overloading: `adata.uns["neighbors"]`"""
        return OverloadedDict(
            ovld.data[key],
            overloaded_keys={
                "connectivities": KeyOverload(
                    "connectivities",
                    get=partial(_adjacency_getter, adata=adata),
                    set=partial(_adjacency_setter, adata=adata),
                ),
                "distances": KeyOverload(
                    "distances",
                    get=partial(_adjacency_getter, adata=adata),
                    set=partial(_adjacency_setter, adata=adata),
                ),
            },
        )

    return OverloadedDict(
        adata._uns,
        overloaded_keys={
            "neighbors": KeyOverload(
                "neighbors",
                get=partial(_neighbors_getter, adata=adata),
                set=partial(_neighbors_setter, adata=adata),
            ),
        },
    )
