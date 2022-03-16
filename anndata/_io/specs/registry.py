from __future__ import annotations

from collections.abc import Mapping, Callable, Iterable
from functools import singledispatch, wraps
from typing import Any, NamedTuple

from anndata.compat import _read_attr, ZarrArray, ZarrGroup, H5Group, H5Array, Literal
from anndata._io.utils import report_write_key_on_error, report_read_key_on_error

# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
# TODO: Should filetype be included in the IOSpec if it changes the encoding? Or does the intent that these things be "the same" overrule that?


class IOSpec(NamedTuple):
    encoding_type: str
    encoding_version: str


class NoSuchIO(KeyError):
    def __init__(
        self,
        registry: Literal["read", "read_partial", "write"],
        key: tuple[type, IOSpec | type | tuple[type, str], frozenset[str]],
    ):
        self.registry = registry
        self.key = key
        super().__init__(str(self))

    def __str__(self) -> str:
        msg = self._get_msg()
        return f"No such {self.registry} function registered: {msg}."

    def _get_msg(self) -> str:
        registry: Mapping[
            tuple[type, IOSpec | type | tuple[type, str], frozenset[str]], Callable
        ] = getattr(_REGISTRY, self.registry)

        typ, spec_or_src_type, modifiers = self.key
        is_read = self.registry != "write"
        assert isinstance(spec_or_src_type, IOSpec) == is_read

        dir_ = "destination" if is_read else "source"
        desc_type = f"{dir_} type {typ.__name__}"
        types = {typ for typ, _, _ in registry}
        if typ not in types:
            return f"Unknown {desc_type}"

        sosts = {sost for typ, sost, _ in registry if typ in types}
        if spec_or_src_type not in sosts:
            if not is_read:
                return f"Destination type {spec_or_src_type} not found for {desc_type}"
            spec: IOSpec = spec_or_src_type
            enc_types = {spec.encoding_type for spec in sosts}
            if spec.encoding_type not in enc_types:
                return f"Unknown encoding type “{spec.encoding_type}” for {desc_type}"
            return (
                f"Unknown encoding version {spec.encoding_version} "
                f"for {desc_type} encoding “{spec.encoding_type}”"
            )

        modifier_sets = {
            mods for typ, sost, mods in registry if typ in types and sost in sosts
        }
        if modifiers not in modifier_sets:
            # “src type y’s encoding IOSpec(...)” or “dest type x and src type y”
            if is_read:
                desc = f"{desc_type}’s encoding {spec_or_src_type}"
            else:
                desc = f"source type {spec_or_src_type} and {desc_type}"
            return f"Unknown modifier set {modifiers} for {desc}"

        assert False, f"Don’t create NoSuchIO error from valid key: {self.key}"


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

        if (dest_type, typ, modifiers) not in self.write:
            raise TypeError(
                f"No method has been defined for writing {typ} elements to {dest_type}"
            )

        try:
            return self.write[(dest_type, typ, modifiers)]
        except KeyError:
            raise NoSuchIO("write", (dest_type, typ, modifiers)) from None

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
            raise NoSuchIO("read", (src_type, spec, modifiers)) from None

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
            raise NoSuchIO("read_partial", (src_type, spec, modifiers)) from None


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
