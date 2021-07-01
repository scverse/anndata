from collections.abc import MutableMapping
from functools import partial
from typing import Any, Callable, List, Mapping, Optional, Union
from warnings import warn
from weakref import proxy


class KeyOverload:
    """
    This class contains the information neccesary to overload a key of a dict.

    It's like a descriptor, but for a key of a dict instead of an attribute.

    Register getter, setter, and deleter methods by passing them at instantiation,
    or assigning them to the `._get`, `._set`, and `._delete` attributes respectivley.
    These functions will be passed the parent `OverloadedDict` and the key as their first
    two arguments. The get and delete methods will be called by the parent with no
    additional arguments, while the setter will be passed the value to set.

    Note that the parent is not set on instantiation. It's currently assumed that's added
    when the parent is constructed.

    Attrs
    -----
    key
        Key in parent dict to overload
    parent
        The parent OverloadedDict this key is attached to.
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
    """A mapping where some of the keys have been overloaded.

    Each overloaded key should be defined as an KeyOverload instance, and can have
    specific getter, settter, and deleter methods. Additionally, overloaded keys don't
    show up in iteration or from `__contains__` calls unless they exist in `.data`.

    Attrs
    -----
    data
        Wrapped mapping.
    overloaded
        Maps from keys to overloaded behaviours.
    """

    data: Mapping
    overloaded: Mapping[Any, KeyOverload]

    def __init__(self, data: Mapping, *, overloaded: Mapping[Any, KeyOverload]):
        self.data = data
        self.overloaded = overloaded
        for v in overloaded.values():
            v.parent = proxy(self)

    def __getitem__(self, key):
        if key in self.overloaded:
            return self.overloaded[key].get()
        else:
            return self.data[key]

    def __setitem__(self, key, value):
        if key in self.overloaded:
            self.overloaded[key].set(value)
        else:
            self.data[key] = value

    def __delitem__(self, key):
        if key in self.overloaded:
            self.overloaded[key].delete()
        else:
            del self.data[key]

    def __contains__(self, key):
        return key in self.data

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return (
            f"OverloadedDict, wrapping:\n\t{self.data!r}\nWith overloaded keys:"
            f"\n\t{list(self.overloaded.keys())}."
        )

    def copy(self) -> dict:
        return self.data.copy()

    def keys(self):
        return self.data.keys()

    def _ipython_key_completions_(self) -> List[str]:
        return list(self.keys())


#######################################
# Handling .uns["neighbors"]
#######################################


def _access_warn(key, cur_loc):
    warn(
        f"This location for '{key}' is deprecated. It has been moved to {cur_loc}, "
        "and will not be accesible here in a future version of anndata.",
        FutureWarning,
        stacklevel=4,
    )


def _adjacency_getter(ovld: OverloadedDict, key, adata: "AnnData"):
    """For overloading:

    >>> mtx = adata.uns["neighbors"]["connectivities"]  # doctest: +SKIP
    >>> mtx = adata.uns["neighbors"]["distances"]  # doctest: +SKIP
    """
    _access_warn(key, f".obsp[{key}]")
    return adata.obsp[key]


def _adjacency_setter(ovld: OverloadedDict, key, value, adata: "AnnData"):
    """For overloading:

    >>> adata.uns["neighbors"]["connectivities"] = mtx  # doctest: +SKIP
    >>> adata.uns["neighbors"]["distances"] = mtx  # doctest: +SKIP
    """
    _access_warn(key, f".obsp[{key}]")
    adata.obsp[key] = value


def _neighbors_setter(ovld: OverloadedDict, key, neighbors: Mapping, adata: "AnnData"):
    """For overloading: `adata.uns["neighbors"] = d`."""
    for k in ("distances", "connectivities"):
        if k in neighbors:
            _access_warn(k, f".obsp[{k}]")
            adata.obsp[k] = neighbors.pop(k)
    ovld.data[key] = neighbors


def _neighbors_getter(ovld: OverloadedDict, key, adata: "AnnData"):
    """For overloading: `adata.uns["neighbors"]`"""
    return OverloadedDict(
        ovld.data[key],
        overloaded={
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


def _overloaded_uns(adata: "AnnData", uns: Union[dict, "DictView"]) -> OverloadedDict:
    return OverloadedDict(
        uns,
        overloaded={
            "neighbors": KeyOverload(
                "neighbors",
                get=partial(_neighbors_getter, adata=adata),
                set=partial(_neighbors_setter, adata=adata),
            ),
        },
    )
