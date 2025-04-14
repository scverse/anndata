from __future__ import annotations

import inspect
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from functools import partial, singledispatch, wraps
from types import MappingProxyType
from typing import TYPE_CHECKING, Generic, TypeVar

from anndata._io.utils import report_read_key_on_error, report_write_key_on_error
from anndata._types import Read, ReadLazy, _ReadInternal, _ReadLazyInternal
from anndata.compat import DaskArray, ZarrGroup, _read_attr, is_zarr_v2

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable
    from typing import Any

    from anndata._types import (
        GroupStorageType,
        ReadCallback,
        StorageType,
        Write,
        WriteCallback,
        _WriteInternal,
    )
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray
    from anndata.typing import RWAble

    from ..._core.xarray import Dataset2D

    T = TypeVar("T")
    W = TypeVar("W", bound=_WriteInternal)
    LazyDataStructures = DaskArray | Dataset2D | CategoricalArray | MaskedArray


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
        cls, dest_type: type, typ: type | tuple[type, str], modifiers: frozenset[str]
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
    def decorator(func: W) -> W:
        @wraps(func)
        def wrapper(g: GroupStorageType, k: str, *args, **kwargs):
            result = func(g, k, *args, **kwargs)
            g[k].attrs.setdefault("encoding-type", spec.encoding_type)
            g[k].attrs.setdefault("encoding-version", spec.encoding_version)
            return result

        return wrapper

    return decorator


_R = TypeVar("_R", _ReadInternal, _ReadLazyInternal)
R = TypeVar("R", Read, ReadLazy)


class IORegistry(Generic[_R, R]):
    def __init__(self):
        self.read: dict[tuple[type, IOSpec, frozenset[str]], _R] = {}
        self.read_partial: dict[tuple[type, IOSpec, frozenset[str]], Callable] = {}
        self.write: dict[
            tuple[type, type | tuple[type, str], frozenset[str]], _WriteInternal
        ] = {}
        self.write_specs: dict[type | tuple[type, str] | tuple[type, type], IOSpec] = {}

    def register_write(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ) -> Callable[[_WriteInternal[T]], _WriteInternal[T]]:
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        # Record specification for src_type
        if src_type in self.write_specs and (spec != self.write_specs[src_type]):
            # First check for consistency
            current_spec = self.write_specs[src_type]
            msg = (
                "Cannot overwrite IO specifications. Attempted to overwrite encoding "
                f"for {src_type} from {current_spec} to {spec}"
            )
            raise TypeError(msg)
        else:
            self.write_specs[src_type] = spec

        def _register(func):
            self.write[(dest_type, src_type, modifiers)] = write_spec(spec)(func)
            return func

        return _register

    def get_write(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        modifiers: frozenset[str] = frozenset(),
        *,
        writer: Writer,
    ) -> Write:
        import h5py

        if dest_type is h5py.File:
            dest_type = h5py.Group
        if (dest_type, src_type, modifiers) not in self.write:
            raise IORegistryError._from_write_parts(dest_type, src_type, modifiers)
        internal = self.write[(dest_type, src_type, modifiers)]
        return partial(internal, _writer=writer)

    def has_write(
        self,
        dest_type: type,
        src_type: type | tuple[type, str],
        modifiers: frozenset[str],
    ) -> bool:
        return (dest_type, src_type, modifiers) in self.write

    def register_read(
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
        modifiers: Iterable[str] = frozenset(),
    ) -> Callable[[_R], _R]:
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read[(src_type, spec, modifiers)] = func
            return func

        return _register

    def get_read(
        self,
        src_type: type,
        spec: IOSpec,
        modifiers: frozenset[str] = frozenset(),
        *,
        reader: Reader,
    ) -> R:
        if (src_type, spec, modifiers) not in self.read:
            raise IORegistryError._from_read_parts("read", self.read, src_type, spec)  # noqa: EM101
        internal = self.read[(src_type, spec, modifiers)]
        return partial(internal, _reader=reader)

    def has_read(
        self, src_type: type, spec: IOSpec, modifiers: frozenset[str] = frozenset()
    ) -> bool:
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

    def get_partial_read(
        self, src_type: type, spec: IOSpec, modifiers: frozenset[str] = frozenset()
    ):
        if (src_type, spec, modifiers) in self.read_partial:
            return self.read_partial[(src_type, spec, modifiers)]
        name = "read_partial"
        raise IORegistryError._from_read_parts(name, self.read_partial, src_type, spec)

    def get_spec(self, elem: Any) -> IOSpec:
        if isinstance(elem, DaskArray):
            if (typ_meta := (DaskArray, type(elem._meta))) in self.write_specs:
                return self.write_specs[typ_meta]
        elif (
            hasattr(elem, "dtype")
            and (typ_kind := (type(elem), elem.dtype.kind)) in self.write_specs
        ):
            return self.write_specs[typ_kind]
        return self.write_specs[type(elem)]


_REGISTRY: IORegistry[_ReadInternal, Read] = IORegistry()
_LAZY_REGISTRY: IORegistry[_ReadLazyInternal, ReadLazy] = IORegistry()


@singledispatch
def proc_spec(spec) -> IOSpec:
    msg = f"proc_spec not defined for type: {type(spec)}."
    raise NotImplementedError(msg)


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
    def __init__(
        self, registry: IORegistry, callback: ReadCallback | None = None
    ) -> None:
        self.registry = registry
        self.callback = callback

    @report_read_key_on_error
    def read_elem(
        self,
        elem: StorageType,
        modifiers: frozenset[str] = frozenset(),
    ) -> RWAble:
        """Read an element from a store. See exported function for more details."""

        iospec = get_spec(elem)
        read_func: Read = self.registry.get_read(
            type(elem), iospec, modifiers, reader=self
        )
        if self.callback is None:
            return read_func(elem)
        return self.callback(read_func, elem.name, elem, iospec=iospec)


class LazyReader(Reader):
    @report_read_key_on_error
    def read_elem(
        self,
        elem: StorageType,
        modifiers: frozenset[str] = frozenset(),
        chunks: tuple[int, ...] | None = None,
        **kwargs,
    ) -> LazyDataStructures:
        """Read a dask element from a store. See exported function for more details."""

        iospec = get_spec(elem)
        read_func: ReadLazy = self.registry.get_read(
            type(elem), iospec, modifiers, reader=self
        )
        if self.callback is not None:
            msg = "Dask reading does not use a callback. Ignoring callback."
            warnings.warn(msg, stacklevel=2)
        read_params = inspect.signature(read_func).parameters
        for kwarg in kwargs:
            if kwarg not in read_params:
                msg = (
                    f"Keyword argument {kwarg} passed to read_elem_lazy are not supported by the "
                    "registered read function."
                )
                raise ValueError(msg)
        if "chunks" in read_params:
            kwargs["chunks"] = chunks
        return read_func(elem, **kwargs)


class Writer:
    def __init__(self, registry: IORegistry, callback: WriteCallback | None = None):
        self.registry = registry
        self.callback = callback

    def find_write_func(
        self, dest_type: type, elem: Any, modifiers: frozenset[str]
    ) -> Write:
        for pattern in _iter_patterns(elem):
            if self.registry.has_write(dest_type, pattern, modifiers):
                return self.registry.get_write(
                    dest_type, pattern, modifiers, writer=self
                )
        # Raises IORegistryError
        return self.registry.get_write(dest_type, type(elem), modifiers, writer=self)

    @report_write_key_on_error
    def write_elem(
        self,
        store: GroupStorageType,
        k: str,
        elem: RWAble,
        *,
        dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
        modifiers: frozenset[str] = frozenset(),
    ):
        from pathlib import PurePosixPath

        import h5py

        # we allow stores to have a prefix like /uns which are then written to with keys like /uns/foo
        if "/" in k.split(store.name)[-1][1:]:
            msg = "Forward slashes are not allowed in keys."
            raise ValueError(msg)

        if isinstance(store, h5py.File):
            store = store["/"]

        dest_type = type(store)

        # Normalize k to absolute path
        if (isinstance(store, ZarrGroup) and is_zarr_v2()) or (
            isinstance(store, h5py.Group) and not PurePosixPath(k).is_absolute()
        ):
            k = str(PurePosixPath(store.name) / k)

        if k == "/":
            if isinstance(store, ZarrGroup) and not is_zarr_v2():
                from zarr.core.sync import sync

                sync(store.store.clear())
            else:
                store.clear()
        elif k in store:
            del store[k]

        write_func = self.find_write_func(dest_type, elem, modifiers)

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


def read_elem(elem: StorageType) -> RWAble:
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


def read_elem_lazy(
    elem: StorageType, chunks: tuple[int, ...] | None = None, **kwargs
) -> LazyDataStructures:
    """
    Read an element from a store lazily.

    Assumes that the element is encoded using the anndata encoding. This function will
    determine the encoded type using the encoding metadata stored in elem's attributes.


    Parameters
    ----------
    elem
        The stored element.
    chunks, optional
       length `n`, the same `n` as the size of the underlying array.
       Note that the minor axis dimension must match the shape for sparse.
       Defaults to `(1000, adata.shape[1])` for CSR sparse,
       `(adata.shape[0], 1000)` for CSC sparse,
       and the on-disk chunking otherwise for dense.
       Can use `-1` or `None` to indicate use of the size of the corresponding dimension.

    Returns
    -------
        A "lazy" elem

    Examples
    --------

    Setting up our example:

    >>> from scanpy.datasets import pbmc3k
    >>> import tempfile
    >>> import anndata as ad
    >>> import zarr

    >>> tmp_path = tempfile.gettempdir()
    >>> zarr_path = tmp_path + "/adata.zarr"

    >>> adata = pbmc3k()
    >>> adata.layers["dense"] = adata.X.toarray()
    >>> adata.write_zarr(zarr_path)

    Reading a sparse matrix from a zarr store lazily, with custom chunk size and default:

    >>> g = zarr.open(zarr_path)
    >>> adata.X = ad.experimental.read_elem_lazy(g["X"])
    >>> adata.X
    dask.array<make_dask_chunk, shape=(2700, 32738), dtype=float32, chunksize=(1000, 32738), chunktype=scipy.csr_matrix>
    >>> adata.X = ad.experimental.read_elem_lazy(g["X"], chunks=(500, adata.shape[1]))
    >>> adata.X
    dask.array<make_dask_chunk, shape=(2700, 32738), dtype=float32, chunksize=(500, 32738), chunktype=scipy.csr_matrix>

    Reading a dense matrix from a zarr store lazily:

    >>> adata.layers["dense"] = ad.experimental.read_elem_lazy(g["layers/dense"])
    >>> adata.layers["dense"]
    dask.array<from-zarr, shape=(2700, 32738), dtype=float32, chunksize=(169, 2047), chunktype=numpy.ndarray>

    Making a new anndata object from on-disk, with custom chunks:

    >>> adata = ad.AnnData(
    ...     obs=ad.io.read_elem(g["obs"]),
    ...     var=ad.io.read_elem(g["var"]),
    ...     uns=ad.io.read_elem(g["uns"]),
    ...     obsm=ad.io.read_elem(g["obsm"]),
    ...     varm=ad.io.read_elem(g["varm"]),
    ... )
    >>> adata.X = ad.experimental.read_elem_lazy(g["X"], chunks=(500, adata.shape[1]))
    >>> adata.layers["dense"] = ad.experimental.read_elem_lazy(g["layers/dense"])

    We also support using -1 and None as a chunk size to signify the reading the whole axis:

    >>> adata.X = ad.experimental.read_elem_lazy(g["X"], chunks=(500, -1))
    >>> adata.X = ad.experimental.read_elem_lazy(g["X"], chunks=(500, None))
    """
    return LazyReader(_LAZY_REGISTRY).read_elem(elem, chunks=chunks, **kwargs)


def write_elem(
    store: GroupStorageType,
    k: str,
    elem: RWAble,
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
    read_partial = _REGISTRY.get_partial_read(
        type(elem), get_spec(elem), frozenset(modifiers)
    )
    return read_partial(elem, items=items, indices=indices)
