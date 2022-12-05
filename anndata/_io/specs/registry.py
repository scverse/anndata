from __future__ import annotations

from abc import abstractmethod, ABC
from collections.abc import Mapping, Callable, Iterable, Set
from functools import singledispatch, wraps
from typing import Any, Literal, NamedTuple, TYPE_CHECKING

if TYPE_CHECKING:
    import h5py
    import zarr

from anndata.compat import _read_attr, ZarrArray, ZarrGroup, H5Group, H5Array
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
    def type0(self) -> type:
        """First type in key: Src type for read, dest type for write"""
        return self.key[0]

    @property
    def spec_or_src_type(self) -> IOSpec | type | tuple[type, str]:
        """IOSpec for read, src type for write"""
        return self.key[1]

    @property
    def modifiers(self) -> frozenset[str]:
        return self.key[2]

    def _get_msg(self) -> str:
        if self.type0 not in self.registry:
            return f"Unknown {self.type0_desc}"

        if self.spec_or_src_type not in self.registry[self.type0]:
            return self._get_spec_or_src_type_msg()

        if self.modifiers not in self.registry[self.type0][self.spec_or_src_type]:
            desc = self._get_modifiers_desc()
            return f"Unknown modifier set {self.modifiers} for {desc}"

        assert False, f"Don’t create NoSuchIO error from valid key: {self.key}"

    @property
    @abstractmethod
    def type0_desc(self) -> str:
        ...

    @abstractmethod
    def _get_spec_or_src_type_msg(self) -> str:
        ...

    @abstractmethod
    def _get_modifiers_desc(self) -> str:
        ...


class NoSuchWrite(NoSuchIO):
    @property
    def src_type(self) -> type | tuple[type, str]:
        return self.spec_or_src_type

    @property
    def type0_desc(self) -> str:
        return f"destination type {self.type0.__name__}"

    def _get_spec_or_src_type_msg(self) -> str:
        return f"Source type {self.src_type} not found for {self.type0_desc}"

    def _get_modifiers_desc(self) -> str:
        return f"{self.type0_desc} and source type {self.src_type}"


class NoSuchRead(NoSuchIO):
    def __str__(self) -> str:
        msg = super().__str__()
        return (
            f"{msg}"
            "You are possibly reading data from a newer anndata version than yours. "
            "You might want to try updating anndata."
        )

    @property
    def spec(self) -> IOSpec:
        return self.spec_or_src_type

    @property
    def type0_desc(self) -> str:
        return f"source type {self.type0.__name__}"

    def _get_spec_or_src_type_msg(self) -> str:
        enc_types = {spec.encoding_type for spec in self.registry[self.type0]}
        if self.spec.encoding_type not in enc_types:
            return (
                f"Unknown encoding type “{self.spec.encoding_type}” "
                f"for {self.type0_desc}"
            )
        return (
            f"Unknown encoding version {self.spec.encoding_version} "
            f"for {self.type0_desc}’s encoding “{self.spec.encoding_type}”"
        )

    def _get_modifiers_desc(self) -> str:
        return f"{self.type0_desc}’s encoding {self.spec}"


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

    def register_write(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

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

        try:
            return self.write[(dest_type, src_type, modifiers)]
        except KeyError:
            raise NoSuchWrite("write", (dest_type, src_type, modifiers)) from None

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
        try:
            return self.read[(src_type, spec, modifiers)]
        except KeyError:
            raise NoSuchRead("read", (src_type, spec, modifiers)) from None

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
    # should be “H5Group | ZarrGroup”, but weirdly that makes Sphinx error
    f: h5py.Group | h5py.File | zarr.hierarchy.Group,
    k: str,
    elem: Any,
    *args,
    modifiers=frozenset(),
    **kwargs,
):
    """
    Write an element to a disk store using its anndata encoding.

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


@report_read_key_on_error
def read_elem(
    elem: H5Array | H5Group | ZarrArray | ZarrGroup,
    *,
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
