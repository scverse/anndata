from __future__ import annotations

from abc import abstractmethod, ABC
from collections.abc import Mapping, Callable, Iterable, Set
from functools import singledispatch, wraps
from typing import Any, NamedTuple

from anndata.compat import _read_attr, ZarrArray, ZarrGroup, H5Group, H5Array, Literal
from anndata._io.utils import report_write_key_on_error, report_read_key_on_error

# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
# TODO: Should filetype be included in the IOSpec if it changes the encoding? Or does the intent that these things be "the same" overrule that?


class IOSpec(NamedTuple):
    encoding_type: str
    encoding_version: str


class NoSuchIO(KeyError, ABC):
    def __init__(
        self,
        registry: Literal["read", "read_partial", "write"],
        key: tuple[type, IOSpec | type | tuple[type, str], frozenset[str]],
    ):
        self.registry_name = registry
        self.key = key
        self.registry = self._hierarchical_registry()
        super().__init__(str(self))

    def _hierarchical_registry(
        self,
    ) -> Mapping[type, Mapping[IOSpec | type | tuple[type, str], Set[frozenset[str]]]]:
        reg = {}
        for typ, spec_or_src_type, modifiers in getattr(_REGISTRY, self.registry_name):
            reg.setdefault(typ, {}).setdefault(spec_or_src_type, set()).add(modifiers)
        return reg

    def __str__(self) -> str:
        return f"No such {self.registry_name} function registered: {self._get_msg()}."

    @property
    def typ(self) -> type:
        return self.key[0]

    @property
    def spec_or_dest_type(self) -> IOSpec | type | tuple[type, str]:
        return self.key[1]

    @property
    def modifiers(self) -> frozenset[str]:
        return self.key[2]

    def _get_msg(self) -> str:
        if self.typ not in self.registry:
            return f"Unknown {self.type_desc}"

        if self.spec_or_dest_type not in self.registry[self.typ]:
            return self._get_spec_or_dest_type_msg()

        if self.modifiers not in self.registry[self.typ][self.spec_or_dest_type]:
            desc = self._get_modifiers_desc()
            return f"Unknown modifier set {self.modifiers} for {desc}"

        assert False, f"Don’t create NoSuchIO error from valid key: {self.key}"

    @property
    @abstractmethod
    def type_desc(self) -> str:
        ...

    @abstractmethod
    def _get_spec_or_dest_type_msg(self) -> str:
        ...

    @abstractmethod
    def _get_modifiers_desc(self) -> str:
        ...


class NoSuchWrite(NoSuchIO):
    @property
    def dest_type(self) -> type | tuple[type, str]:
        return self.spec_or_dest_type

    @property
    def type_desc(self) -> str:
        return f"source type {self.typ.__name__}"

    def _get_spec_or_dest_type_msg(self) -> str:
        return f"Destination type {self.dest_type} not found for {self.type_desc}"

    def _get_modifiers_desc(self) -> str:
        return f"{self.type_desc} and destination type {self.dest_type}"


class NoSuchRead(NoSuchIO):
    @property
    def spec(self) -> IOSpec:
        return self.spec_or_dest_type

    @property
    def type_desc(self) -> str:
        return f"destination type {self.typ.__name__}"

    def _get_spec_or_dest_type_msg(self) -> str:
        enc_types = {spec.encoding_type for spec in self.registry[self.typ]}
        if self.spec.encoding_type not in enc_types:
            return (
                f"Unknown encoding type “{self.spec.encoding_type}” "
                f"for {self.type_desc}"
            )
        return (
            f"Unknown encoding version {self.spec.encoding_version} "
            f"for {self.type_desc}’s encoding “{self.spec.encoding_type}”"
        )

    def _get_modifiers_desc(self) -> str:
        return f"{self.type_desc}’s encoding {self.spec}"


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


class IORegistry(object):
    def __init__(self):
        self.read: dict[tuple[type, IOSpec, frozenset[str]], Callable] = {}
        self.read_partial: dict[tuple[type, IOSpec, frozenset[str]], Callable] = {}
        self.write: dict[
            tuple[type, type | tuple[type, str], frozenset[str]], Callable
        ] = {}

    def register_write(
        self,
        dest_type: type,
        typ: type | tuple[type, str],
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.write[(dest_type, typ, modifiers)] = write_spec(spec)(func)
            return func

        return _register

    def get_writer(
        self,
        dest_type: type,
        typ: type | tuple[type, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        import h5py

        if dest_type is h5py.File:
            dest_type = h5py.Group
        modifiers = frozenset(modifiers)

        try:
            return self.write[(dest_type, typ, modifiers)]
        except KeyError:
            raise NoSuchWrite("write", (dest_type, typ, modifiers)) from None

    def has_writer(
        self, dest_type: type, typ: type | tuple[type, str], modifiers: Iterable[str]
    ):
        modifiers = frozenset(modifiers)
        return (dest_type, typ, modifiers) in self.write

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
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        modifiers = frozenset(modifiers)
        try:
            return self.read[(src_type, spec, modifiers)]
        except KeyError:
            raise NoSuchRead("read", (src_type, spec, modifiers)) from None

    def has_reader(
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        modifiers = frozenset(modifiers)
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
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        modifiers = frozenset(modifiers)
        try:
            return self.read_partial[(src_type, spec, modifiers)]
        except KeyError:
            raise NoSuchRead("read_partial", (src_type, spec, modifiers)) from None


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
    elem: H5Array | H5Group | ZarrArray | ZarrGroup,
) -> IOSpec:
    return proc_spec(
        {
            k: _read_attr(elem.attrs, k, "")
            for k in ["encoding-type", "encoding-version"]
        }
    )


@report_write_key_on_error
def write_elem(
    f: H5Group | ZarrGroup,
    k: str,
    elem: Any,
    *args,
    modifiers=frozenset(),
    **kwargs,
):
    """
    Write an element to a disk store using it's anndata encoding.

    Params
    ------
    f
        The store to write to.
    k
        The key to write for this value.
    elem
        The element to write as k to f.
    """
    dest_type = type(f)
    if elem is None:
        return
    t = type(elem)
    if k == "/":
        f.clear()
    elif k in f:
        del f[k]
    if (
        hasattr(elem, "dtype")
        and (dest_type, (t, elem.dtype.kind), modifiers) in _REGISTRY.write
    ):
        _REGISTRY.get_writer(dest_type, (t, elem.dtype.kind), modifiers)(
            f, k, elem, *args, **kwargs
        )
    else:
        _REGISTRY.get_writer(dest_type, t, modifiers)(f, k, elem, *args, **kwargs)


def read_elem(
    elem: H5Array | H5Group | ZarrArray | ZarrGroup,
    modifiers: frozenset[str] = frozenset(),
) -> Any:
    """Read an element from an on disk store."""
    return _REGISTRY.get_reader(type(elem), get_spec(elem), frozenset(modifiers))(elem)


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
