from __future__ import annotations

from collections.abc import Mapping
from functools import singledispatch, wraps
from typing import Any, NamedTuple, Tuple, Type, Callable, Union


from anndata.compat import _read_attr, ZarrArray, ZarrGroup, H5Group, H5Array
from anndata._io.utils import report_write_key_on_error, report_read_key_on_error

# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
# TODO: Should filetype be included in the IOSpec if it changes the encoding? Or does the intent that these things be "the same" overrule that?


class IOSpec(NamedTuple):
    encoding_type: str
    encoding_version: str


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
        self.read: Mapping[tuple[str, IOSpec], Callable] = {}
        self.read_partial: Mapping[Tuple[str, IOSpec], Callable] = {}
        self.write: Mapping[Union[Type, Tuple[Type, str]], Callable] = {}

    def register_write(
        self,
        dest_type,
        typ: Union[type, tuple[type, str]],
        spec,
        modifiers: frozenset(str) = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.write[(dest_type, typ, modifiers)] = write_spec(spec)(func)
            return func

        return _register

    def get_writer(self, dest_type, typ, modifiers=frozenset()):
        import h5py

        if dest_type is h5py.File:
            dest_type = h5py.Group
        modifiers = frozenset(modifiers)

        if (dest_type, typ, modifiers) not in self.write:
            raise TypeError(
                f"No method has been defined for writing {typ} elements to {dest_type}"
            )

        return self.write[(dest_type, typ, modifiers)]

    def has_writer(self, dest_type, typ, modifiers):
        modifiers = frozenset(modifiers)
        return (dest_type, typ, modifiers) in self.write

    def register_read(self, src_type, spec, modifiers: frozenset[str] = frozenset()):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read[(src_type, spec, modifiers)] = func
            return func

        return _register

    def get_reader(self, src_type, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return self.read[(src_type, spec, modifiers)]

    def has_reader(self, src_type, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return (src_type, spec, modifiers) in self.read

    def register_read_partial(
        self, src_type, spec, modifiers: frozenset[str] = frozenset()
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read_partial[(src_type, spec, modifiers)] = func
            return func

        return _register

    def get_partial_reader(self, src_type, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return self.read_partial[(src_type, spec, modifiers)]


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
    elem: "Union[h5py.Dataset, h5py.Group, zarr.Group, zarr.Dataset]",
) -> IOSpec:
    return proc_spec(
        {
            k: _read_attr(elem.attrs, k, "")
            for k in ["encoding-type", "encoding-version"]
        }
    )


@report_write_key_on_error
def write_elem(
    f: "Union[H5Group, ZarrGroup]",
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
    elem: Union[H5Array, H5Group, ZarrGroup, ZarrArray],
    modifiers: frozenset(str) = frozenset(),
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
