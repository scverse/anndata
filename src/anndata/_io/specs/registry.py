from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import singledispatch, wraps
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from anndata._io.utils import report_read_key_on_error, report_write_key_on_error
from anndata.compat import DaskArray, _read_attr

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable

    from anndata._types import GroupStorageType, StorageType


# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
# TODO: Should filetype be included in the IOSpec if it changes the encoding? Or does the intent that these things be "the same" overrule that?
@dataclass(frozen=True)
class IOSpec:
    encoding_type: str
    encoding_version: str


# TODO: Should this subclass from LookupError?
class IORegistryError(Exception):
    @classmethod
    def _from_write_parts(
        cls, dest_type: type, typ: type, modifiers: frozenset[str]
    ) -> IORegistryError:
        msg = f"No method registered for writing {typ} into {dest_type}"
        if modifiers:
            msg += f" with {modifiers}"
        return cls(msg)

    @classmethod
    def _from_read_parts(
        cls,
        method: str,
        registry: Mapping,
        src_typ: type[StorageType],
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
        def wrapper(g: GroupStorageType, k: str, *args, **kwargs):
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
        self.write_specs: dict[type | tuple[type, str] | tuple[type, type], IOSpec] = {}

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
        if isinstance(elem, DaskArray):
            if (typ_meta := (DaskArray, type(elem._meta))) in self.write_specs:
                return self.write_specs[typ_meta]
        elif hasattr(elem, "dtype"):
            if (typ_kind := (type(elem), elem.dtype.kind)) in self.write_specs:
                return self.write_specs[typ_kind]
        return self.write_specs[type(elem)]


_REGISTRY = IORegistry()


@singledispatch
def proc_spec(spec) -> IOSpec:
    raise NotImplementedError(f"proc_spec not defined for type: {type(spec)}.")


@proc_spec.register(IOSpec)
def proc_spec_spec(spec: IOSpec) -> IOSpec:
    return spec


@proc_spec.register(Mapping)
def proc_spec_mapping(spec: Mapping[str, str]) -> IOSpec:
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


def _iter_patterns(
    elem,
) -> Generator[tuple[type, type | str] | tuple[type, type, str], None, None]:
    """Iterates over possible patterns for an element in order of precedence."""
    from anndata.compat import DaskArray

    t = type(elem)

    if isinstance(elem, DaskArray):
        yield (t, type(elem._meta), elem.dtype.kind)
        yield (t, type(elem._meta))
    if hasattr(elem, "dtype"):
        yield (t, elem.dtype.kind)
    yield t


class Reader:
    def __init__(self, registry: IORegistry, callback: Callable | None = None) -> None:
        self.registry = registry
        self.callback = callback

    @report_read_key_on_error
    def read_elem(
        self,
        elem: StorageType,
        modifiers: frozenset[str] = frozenset(),
    ) -> Any:
        """Read an element from a store. See exported function for more details."""
        from functools import partial

        iospec = get_spec(elem)
        read_func = partial(
            self.registry.get_reader(type(elem), iospec, modifiers),
            _reader=self,
        )
        if self.callback is None:
            return read_func(elem)
        return self.callback(read_func, elem.name, elem, iospec=iospec)


class Writer:
    def __init__(self, registry: IORegistry, callback: Callable | None = None):
        self.registry = registry
        self.callback = callback

    def find_writer(self, dest_type: type, elem, modifiers: frozenset[str]):
        for pattern in _iter_patterns(elem):
            if self.registry.has_writer(dest_type, pattern, modifiers):
                return self.registry.get_writer(dest_type, pattern, modifiers)
        # Raises IORegistryError
        return self.registry.get_writer(dest_type, type(elem), modifiers)

    @report_write_key_on_error
    def write_elem(
        self,
        store: GroupStorageType,
        k: str,
        elem: Any,
        *,
        dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
        modifiers: frozenset[str] = frozenset(),
    ):
        from functools import partial
        from pathlib import PurePosixPath

        import h5py

        if isinstance(store, h5py.File):
            store = store["/"]

        dest_type = type(store)

        if elem is None:
            return lambda *_, **__: None

        # Normalize k to absolute path
        if not PurePosixPath(k).is_absolute():
            k = str(PurePosixPath(store.name) / k)

        if k == "/":
            store.clear()
        elif k in store:
            del store[k]

        write_func = partial(
            self.find_writer(dest_type, elem, modifiers),
            _writer=self,
        )

        if self.callback is None:
            return write_func(store, k, elem, dataset_kwargs=dataset_kwargs)
        return self.callback(
            write_func,
            store,
            k,
            elem,
            dataset_kwargs=dataset_kwargs,
            iospec=self.registry.get_spec(elem),
        )


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
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
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
