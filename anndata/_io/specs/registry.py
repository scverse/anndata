from __future__ import annotations

from collections.abc import Mapping, Callable, Iterable
from dataclasses import dataclass
from functools import singledispatch, wraps
from types import MappingProxyType
from typing import Any, Union

from anndata.compat import _read_attr
from anndata._types import StorageType, GroupStorageType
from anndata._io.utils import report_write_key_on_error, report_read_key_on_error

# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
# TODO: Should filetype be included in the IOSpec if it changes the encoding? Or does the intent that these things be "the same" overrule that?


@dataclass(frozen=True)
class IOSpec:
    encoding_type: str
    encoding_version: str


# TODO: Should this subclass from LookupError?
class IORegistryError(Exception):
    @classmethod
    def _from_write_parts(cls, dest_type, typ, modifiers) -> IORegistryError:
        msg = f"No method registered for writing {typ} into {dest_type}"
        if modifiers:
            msg += f" with {modifiers}"
        return cls(msg)

    @classmethod
    def _from_read_parts(
        cls,
        method: str,
        registry: Mapping,
        src_typ: StorageType,
        spec: IOSpec,
    ) -> IORegistryError:
        # TODO: Improve error message if type exists, but version does not
        msg = (
            f"No {method} method registered for {spec} from {src_typ}. "
            "You may need to update your installation of anndata."
        )
        return cls(msg)


def write_spec(spec: IOSpec):
    def decorator(func: Callable):
        @wraps(func)
        def wrapper(g, k, *args, **kwargs):
            result = func(g, k, *args, **kwargs)
            g[k].attrs.setdefault("encoding-type", spec.encoding_type)
            g[k].attrs.setdefault("encoding-version", spec.encoding_version)
            return result

        return wrapper

    return decorator


class IORegistry:
    def __init__(self):
        self.read: dict[tuple[type, IOSpec, frozenset[str]], Callable] = {}
        self.read_partial: dict[tuple[type, IOSpec, frozenset[str]], Callable] = {}
        self.write: dict[
            tuple[type, type | tuple[type, str], frozenset[str]], Callable
        ] = {}
        self.write_specs: dict[Union[type, tuple[type, str]], IOSpec] = {}

    def register_write(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        # Record specification for src_type
        if src_type in self.write_specs and (spec != self.write_specs[src_type]):
            # First check for consistency
            current_spec = self.write_specs[src_type]
            raise TypeError(
                "Cannot overwrite IO specifications. Attempted to overwrite encoding "
                f"for {src_type} from {current_spec} to {spec}"
            )
        else:
            self.write_specs[src_type] = spec

        def _register(func):
            self.write[(dest_type, src_type, modifiers)] = write_spec(spec)(func)
            return func

        return _register

    def get_writer(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        modifiers: frozenset[str] = frozenset(),
    ):
        import h5py

        if dest_type is h5py.File:
            dest_type = h5py.Group

        if (dest_type, src_type, modifiers) in self.write:
            return self.write[(dest_type, src_type, modifiers)]
        else:
            raise IORegistryError._from_write_parts(dest_type, src_type, modifiers)

    def has_writer(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        modifiers: frozenset[str],
    ):
        return (dest_type, src_type, modifiers) in self.write

    def register_read(
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read[(src_type, spec, modifiers)] = func
            return func

        return _register

    def get_reader(
        self, src_type: type, spec: IOSpec, modifiers: frozenset[str] = frozenset()
    ):
        if (src_type, spec, modifiers) in self.read:
            return self.read[(src_type, spec, modifiers)]
        else:
            raise IORegistryError._from_read_parts(
                "read", _REGISTRY.read, src_type, spec
            )

    def has_reader(
        self, src_type: type, spec: IOSpec, modifiers: frozenset[str] = frozenset()
    ):
        return (src_type, spec, modifiers) in self.read

    def register_read_partial(
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read_partial[(src_type, spec, modifiers)] = func
            return func

        return _register

    def get_partial_reader(
        self, src_type: type, spec: IOSpec, modifiers: frozenset[str] = frozenset()
    ):
        if (src_type, spec, modifiers) in self.read_partial:
            return self.read_partial[(src_type, spec, modifiers)]
        else:
            raise IORegistryError._from_read_parts(
                "read_partial", _REGISTRY.read_partial, src_type, spec
            )

    def get_spec(self, elem: Any) -> IOSpec:
        if hasattr(elem, "dtype"):
            typ = (type(elem), elem.dtype.kind)
            if typ in self.write_specs:
                return self.write_specs[typ]
        return self.write_specs[type(elem)]


_REGISTRY = IORegistry()


@singledispatch
def proc_spec(spec) -> IOSpec:
    raise NotImplementedError(f"proc_spec not defined for type: {type(spec)}.")


@proc_spec.register(IOSpec)
def proc_spec_spec(spec) -> IOSpec:
    return spec


@proc_spec.register(Mapping)
def proc_spec_mapping(spec) -> IOSpec:
    return IOSpec(**{k.replace("-", "_"): v for k, v in spec.items()})


def get_spec(
    elem: StorageType,
) -> IOSpec:
    return proc_spec(
        {
            k: _read_attr(elem.attrs, k, "")
            for k in ["encoding-type", "encoding-version"]
        }
    )


class Reader:
    def __init__(
        self, registry: IORegistry, callback: Union[Callable, None] = None
    ) -> None:
        self.registry = registry
        self.callback = callback

    @report_read_key_on_error
    def read_elem(
        self,
        elem: StorageType,
        modifiers: frozenset(str) = frozenset(),
    ) -> Any:
        """Read an element from a store. See exported function for more details."""
        from functools import partial

        read_func = self.registry.get_reader(
            type(elem), get_spec(elem), frozenset(modifiers)
        )
        read_func = partial(read_func, _reader=self)
        if self.callback is not None:
            return self.callback(read_func, elem.name, elem, iospec=get_spec(elem))
        else:
            return read_func(elem)


class Writer:
    def __init__(
        self,
        registry: IORegistry,
        callback: Union[
            Callable[
                [
                    GroupStorageType,
                    str,
                    StorageType,
                    dict,
                ],
                None,
            ],
            None,
        ] = None,
    ):
        self.registry = registry
        self.callback = callback

    @report_write_key_on_error
    def write_elem(
        self,
        store: GroupStorageType,
        k: str,
        elem,
        *,
        dataset_kwargs=MappingProxyType({}),
        modifiers=frozenset(),
    ):
        from functools import partial
        from pathlib import PurePath

        dest_type = type(store)
        t = type(elem)

        if elem is None:
            return lambda *_, **__: None

        # Normalize k to abosulte path
        if not PurePath(k).is_absolute():
            k = str(PurePath(store.name) / k)

        if k == "/":
            store.clear()
        elif k in store:
            del store[k]
        if (
            hasattr(elem, "dtype")
            and (dest_type, (t, elem.dtype.kind), modifiers) in self.registry.write
        ):
            write_func = partial(
                self.registry.get_writer(dest_type, (t, elem.dtype.kind), modifiers),
                _writer=self,
            )
        else:
            write_func = partial(
                self.registry.get_writer(dest_type, t, modifiers),
                _writer=self,
            )

        if self.callback is not None:
            return self.callback(
                write_func,
                store,
                k,
                elem,
                dataset_kwargs=dataset_kwargs,
                iospec=self.registry.get_spec(elem),
            )
        else:
            return write_func(store, k, elem, dataset_kwargs=dataset_kwargs)


def read_elem(elem: StorageType) -> Any:
    """
    Read an element from a store.

    Assumes that the element is encoded using the anndata encoding. This function will
    determine the encoded type using the encoding metadata stored in elem's attributes.

    Params
    ------
    elem
        The stored element.
    """
    return Reader(_REGISTRY).read_elem(elem)


def write_elem(
    store: GroupStorageType,
    k: str,
    elem: Any,
    *,
    dataset_kwargs: Mapping = MappingProxyType({}),
) -> None:
    """
    Write an element to a storage group using anndata encoding.

    Params
    ------
    store
        The group to write to.
    k
        The key to write to in the group. Note that absolute paths will be written
        from the root.
    elem
        The element to write. Typically an in-memory object, e.g. an AnnData, pandas
        dataframe, scipy sparse matrix, etc.
    dataset_kwargs
        Keyword arguments to pass to the stores dataset creation function.
        E.g. for zarr this would be `chunks`, `compressor`.
    """
    Writer(_REGISTRY).write_elem(store, k, elem, dataset_kwargs=dataset_kwargs)


# TODO: If all items would be read, just call normal read method
def read_elem_partial(
    elem,
    *,
    items=None,
    indices=(slice(None), slice(None)),
    modifiers: frozenset[str] = frozenset(),
):
    """Read part of an element from an on disk store."""
    return _REGISTRY.get_partial_reader(
        type(elem), get_spec(elem), frozenset(modifiers)
    )(elem, items=items, indices=indices)


@singledispatch
def elem_key(elem) -> str:
    return elem.name


#     raise NotImplementedError()

# @elem_key.register(ZarrGroup)
# @elem_key.register(ZarrArray)
# def _(elem):
#     return elem.name

# @elem_key.register(H5Array)
# @elem_key.register(H5Group)
# def _(elem):
#     re
