from __future__ import annotations

import inspect
from collections.abc import Mapping
from dataclasses import dataclass
from functools import partial, singledispatch, wraps
from types import MappingProxyType
from typing import TYPE_CHECKING, cast

import numpy as np
import zarr

from anndata._io.utils import (
    IndexedLayout,
    report_read_key_on_error,
    report_write_key_on_error,
)
from anndata._settings import settings
from anndata._types import Read, ReadLazy, _ReadInternal, _ReadLazyInternal
from anndata._warnings import WriteWarning
from anndata._write_compat import WriteCompat, WriteFeatures
from anndata.compat import DaskArray, _read_attr, has_xp

from ...utils import warn

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable
    from typing import Any, Literal

    from anndata._types import (
        ReadCallback,
        StorageType,
        Write,
        WriteCallback,
        _GroupStorageType,
        _WriteInternal,
    )
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray
    from anndata.types import WriteCompatStr
    from anndata.typing import RWAble

    from ..._core.xarray import Dataset2D

    type LazyDataStructures = DaskArray | Dataset2D | CategoricalArray | MaskedArray

    type WriteSrcType = (
        type | tuple[type, str] | tuple[type, type] | tuple[type, type, str]
    )


def to_writeable(x):
    # Convert non-numpy array-API arrays to numpy via DLPack. Array-API arrays
    # that cannot export via DLPack (e.g. pydata/sparse) are left as-is so the
    # registry can dispatch on their concrete type (or raise a clear error).
    if has_xp(x) and not (isinstance(x, np.ndarray) or np.isscalar(x)):
        return np.from_dlpack(x)
    return x


def normalize_nested(obj):
    if isinstance(obj, dict):
        return {k: normalize_nested(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return type(obj)(normalize_nested(v) for v in obj)
    return to_writeable(obj)


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
        cls,
        dest_type: type,
        typ: WriteSrcType,
        *,
        whens: Iterable[WriteFeatures] = (),
    ) -> IORegistryError:
        msg = f"No method registered for writing {typ} into {dest_type}"
        if whens := list(whens):
            reqs = " or ".join(" and ".join(f"`{flag}`" for flag in w) for w in whens)
            msg += f" with the current `compat`: needs flag(s) {reqs}"
            # `WriteCompat` iterates oldest first
            if profile := next(
                (p for p in WriteCompat if any(w in p.features for w in whens)), None
            ):
                msg += f'. Pass `compat="{profile}"` to the write function'
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


def write_spec[W: Callable[..., Any]](spec: IOSpec) -> Callable[[W], W]:
    def decorator(func: W) -> W:
        @wraps(func)
        def wrapper(g: _GroupStorageType, k: str, *args, **kwargs):
            result = func(g, k, *args, **kwargs)
            g[k].attrs.setdefault("encoding-type", spec.encoding_type)
            g[k].attrs.setdefault("encoding-version", spec.encoding_version)
            return result

        # a wrapper can never structurally match the wrapped function’s type
        return cast("W", wrapper)

    return decorator


class _ByFlags[T](dict[WriteFeatures, T]):
    """Values by the :class:`~anndata.WriteFeatures` flags they require."""

    def select(self, compat: WriteFeatures) -> T | None:
        """Pick the value requiring the most of `compat`’s enabled flags."""
        matching = [when for when in self if when in compat]
        if not matching:
            return None
        best = max(matching, key=len)
        if ties := [w for w in matching if len(w) == len(best) and w != best]:
            msg = f"Ambiguous writers for {compat!r}: {[best, *ties]}"
            raise IORegistryError(msg)
        return self[best]


class IORegistry[RI: (_ReadInternal, _ReadLazyInternal), R: (Read, ReadLazy)]:
    read: dict[tuple[type, IOSpec], RI]
    read_partial: dict[tuple[type, IOSpec], Callable]
    write: dict[tuple[type, WriteSrcType], _ByFlags[_WriteInternal]]
    """Writers by destination, source, then by required flags."""
    write_specs: dict[WriteSrcType, _ByFlags[IOSpec]]
    """Specs by source, then by required flags."""

    def __init__(self) -> None:
        self.read = {}
        self.read_partial = {}
        self.write = {}
        self.write_specs = {}

    def register_write[W: _WriteInternal](
        self,
        dest_type: type[StorageType],
        src_type: WriteSrcType,
        spec: IOSpec | Mapping[str, str],
        *,
        when: WriteFeatures = WriteFeatures(0),
    ) -> Callable[[W], W]:
        """Register a writer.

        `when` are :class:`~anndata.WriteFeatures` flags that all need to be enabled
        for this writer to be used. If several writers apply, the one with the most flags wins.
        """
        spec = proc_spec(spec)

        # Record specification for src_type
        specs = self.write_specs.setdefault(src_type, _ByFlags())
        if when in specs and (spec != specs[when]):
            # First check for consistency
            msg = (
                "Cannot overwrite IO specifications. Attempted to overwrite encoding "
                f"for {src_type} from {specs[when]} to {spec}"
            )
            raise TypeError(msg)
        specs[when] = spec

        def _register(func: W) -> W:
            writers = self.write.setdefault((dest_type, src_type), _ByFlags())
            writers[when] = write_spec(spec)(func)
            return func

        return _register

    def get_write(
        self,
        dest_type: type,
        src_type: WriteSrcType,
        *,
        writer: Writer,
    ) -> Write:
        import h5py

        if dest_type is h5py.File:
            dest_type = h5py.Group
        candidates = self.write.get((dest_type, src_type), _ByFlags())
        if (internal := candidates.select(writer.compat)) is None:
            raise IORegistryError._from_write_parts(
                dest_type, src_type, whens=candidates.keys()
            )
        return partial(internal, _writer=writer)

    def has_write(self, dest_type: type, src_type: WriteSrcType) -> bool:
        return (dest_type, src_type) in self.write

    def writeable_types(
        self, compat: WriteFeatures, store_type: Literal["h5", "zarr"] | None = None
    ) -> set[WriteSrcType]:
        """Source types that have a writer for `compat` (and `store_type`, if given)."""
        return {
            src_type
            for (dest_type, src_type), writers in self.write.items()
            if store_type is None or store_type in dest_type.__module__
            if writers.select(compat) is not None
        }

    def register_read(
        self,
        src_type: type,
        spec: IOSpec | Mapping[str, str],
    ) -> Callable[[RI], RI]:
        spec = proc_spec(spec)

        def _register(func: RI) -> RI:
            self.read[src_type, spec] = func
            return func

        return _register

    def get_read(
        self,
        src_type: type,
        spec: IOSpec,
        *,
        reader: Reader,
    ) -> R:
        if (src_type, spec) not in self.read:
            raise IORegistryError._from_read_parts("read", self.read, src_type, spec)  # noqa: EM101
        # the registry cannot express that `RI`’s `_reader` type matches `reader`’s
        internal: Callable[..., Any] = self.read[src_type, spec]
        return partial(internal, _reader=reader)

    def has_read(self, src_type: type, spec: IOSpec) -> bool:
        return (src_type, spec) in self.read

    def register_read_partial(self, src_type: type, spec: IOSpec | Mapping[str, str]):
        spec = proc_spec(spec)

        def _register(func):
            self.read_partial[src_type, spec] = func
            return func

        return _register

    def get_partial_read(self, src_type: type, spec: IOSpec):
        if (src_type, spec) in self.read_partial:
            return self.read_partial[src_type, spec]
        name = "read_partial"
        raise IORegistryError._from_read_parts(name, self.read_partial, src_type, spec)

    def get_spec(
        self, elem: StorageType, compat: WriteFeatures = WriteCompat.DEFAULT.features
    ) -> IOSpec:
        patterns: list[WriteSrcType] = []
        if isinstance(elem, DaskArray):
            patterns.append((DaskArray, type(elem._meta)))
        elif hasattr(elem, "dtype"):
            patterns.append((type(elem), elem.dtype.kind))
        for pattern in [*patterns, *_iter_patterns(elem)]:
            specs = self.write_specs.get(pattern, _ByFlags())
            if (spec := specs.select(compat)) is not None:
                return spec
        raise KeyError(type(elem))


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


def get_spec(elem: StorageType) -> IOSpec:
    return proc_spec({
        k: _read_attr(elem.attrs, k, "") for k in ["encoding-type", "encoding-version"]
    })


def _iter_patterns(elem: RWAble) -> Generator[WriteSrcType]:
    """Iterates over possible patterns for an element in order of precedence."""
    from anndata.compat import DaskArray

    t = type(elem)

    if isinstance(elem, DaskArray):
        yield (t, type(elem._meta), elem.dtype.kind)
        yield (t, type(elem._meta))
    # Array API dtypes don’t have guaranteed attributes
    if isinstance(elem, np.ndarray):
        yield (t, elem.dtype.kind)
    yield t
    # accessor subclasses (e.g. from other packages) are written by their closest base’s writer
    from anndata.acc import AdRef, RefAcc

    yield from (b for b in t.__mro__[1:] if b in {AdRef, RefAcc})


class Reader:
    def __init__(
        self, registry: IORegistry, callback: ReadCallback | None = None
    ) -> None:
        self.registry = registry
        self.callback = callback

    @report_read_key_on_error
    def read_elem(self, elem: StorageType) -> RWAble:
        """Read an element from a store. See exported function for more details."""

        iospec = get_spec(elem)
        read_func: Read = self.registry.get_read(type(elem), iospec, reader=self)
        if self.callback is None:
            return read_func(elem)
        return self.callback(read_func, elem.name, elem, iospec=iospec)


class LazyReader(Reader):
    @report_read_key_on_error
    def read_elem(
        self,
        elem: StorageType,
        chunks: tuple[int | None, ...] | None = None,
        **kwargs,
    ) -> LazyDataStructures:
        """Read a dask element from a store. See exported function for more details."""

        iospec = get_spec(elem)
        read_func: ReadLazy = self.registry.get_read(type(elem), iospec, reader=self)
        if self.callback is not None:
            msg = "Dask reading does not use a callback. Ignoring callback."
            warn(msg, UserWarning)
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
    def __init__(
        self,
        registry: IORegistry,
        callback: WriteCallback | None = None,
        *,
        compat: WriteCompat | WriteCompatStr | WriteFeatures = WriteCompat.DEFAULT,
    ):
        self.registry = registry
        self.callback = callback
        self.compat = WriteFeatures(compat)

    def find_write_func(self, dest_type: type, elem: RWAble) -> Write:
        for pattern in _iter_patterns(elem):
            if self.registry.has_write(dest_type, pattern):
                return self.registry.get_write(dest_type, pattern, writer=self)
        # Raises IORegistryError
        return self.registry.get_write(dest_type, type(elem), writer=self)

    @report_write_key_on_error
    def write_elem(
        self,
        store: _GroupStorageType,
        k: str,
        elem: RWAble,
        *,
        dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    ):
        from pathlib import PurePosixPath

        import h5py

        from anndata._io.zarr import is_group_consolidated

        if isinstance(store, h5py.File):
            store = store["/"]
        elif is_group_consolidated(store, strict=False):
            msg = "Cannot overwrite/edit a store with consolidated metadata"
            raise ValueError(msg)

        if k == "/":
            if store.name != "/":
                msg = f"'/' is not in the subpath of {store.name!r}"
                raise ValueError(msg)

            if isinstance(store, zarr.Group):
                from zarr.core.sync import sync

                sync(store.store.clear())
            else:
                store.clear()
        else:
            # we allow stores to have a prefix like /uns which are then written to with keys like /uns/foo
            if k.startswith("/"):
                k = str(PurePosixPath(k).relative_to(store.name, walk_up=False))

            _check_child_key(store, k)
            if k in store:
                del store[k]

        # Normalize array-API (e.g., JAX/CuPy) even if not AnnData
        elem = normalize_nested(elem)

        write_func = self.find_write_func(type(store), elem)

        if self.callback is not None:
            write_func = partial(
                self.callback,
                write_func,
                iospec=self.registry.get_spec(elem, self.compat),
            )
        write_func(store, k, elem, dataset_kwargs=dataset_kwargs)
        _add_to_column_order(store, k)
        _warn_unindexed_child(store, k)


def _check_child_key(store: _GroupStorageType, k: str) -> None:
    hint = 'Write the whole mapping or dataframe with `compat="0.14"` to allow any string as key'
    if not k:  # would refer to `store` itself
        msg = f"Empty keys are not allowed in {type(store)}. {hint}"
        raise ValueError(msg)

    # Apart from this code, we also ban keys containing slashes in `write_adata`/`write_h5ad`
    # for AnnData elements other than `obs`, `var`, and `uns`.
    if "/" in k:
        if isinstance(store, zarr.Group) or settings.disallow_forward_slash_in_h5ad:
            msg = f"Forward slashes are not allowed in keys in {type(store)}. {hint}"
            raise ValueError(msg)
        msg = "Forward slashes will be written differently in a future anndata version"
        warn(msg, FutureWarning)


def _add_to_column_order(store: _GroupStorageType, k: str) -> None:
    """Register a column written into an existing dataframe group in its column order."""
    if (
        k != "/"  # writing the dataframe itself to the root
        and store.attrs.get("encoding-type") == "dataframe"
        and "column-order" in store.attrs  # not in the 0.3.0 layout
        and k != store.attrs["_index"]
        and k not in (order := list(cast("Iterable[str]", store.attrs["column-order"])))
    ):
        store.attrs["column-order"] = [*order, k]


def _warn_unindexed_child(store: _GroupStorageType, k: str) -> None:
    """Warn when writing a child that readers of an indexed group will ignore."""
    if (
        k == "/"
        or (layout := IndexedLayout.of(store)) is None
        or k == layout.names_key
        or (
            (prefix := f"{layout.prefix}_") == k[: len(prefix)]
            and k[len(prefix) :].isdecimal()
        )
    ):
        return
    msg = (
        f"Wrote {k!r} into {store.name!r}, which stores its names in "
        f"{layout.names_attr!r} or {layout.names_key!r}, so readers will ignore it. "
        f"Write the whole {layout.encoding_type} instead, "
        "or use `AnnData.write_*(…, use=…)` for elements of an AnnData store."
    )
    warn(msg, WriteWarning)


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
    elem: StorageType, chunks: tuple[int | None, ...] | None = None, **kwargs
) -> LazyDataStructures:
    """
    Read an element from a store lazily.

    Assumes that the element is encoded using the anndata encoding. This function will
    determine the encoded type using the encoding metadata stored in elem's attributes.


    Parameters
    ----------
    elem
        The stored element.
    chunks
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

    ..
        TODO: remove “SKIP” once https://github.com/zarr-developers/zarr-python/issues/3602 becomes minimum zarr (3.1.6)

    >>> adata.layers["dense"] = ad.experimental.read_elem_lazy(g["layers/dense"])
    >>> adata.layers["dense"]  # doctest: +SKIP
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
    store: _GroupStorageType,
    k: str,
    elem: RWAble,
    *,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    compat: WriteCompat | WriteCompatStr | WriteFeatures = WriteCompat.DEFAULT,
) -> None:
    """
    Write an element to a storage group using anndata encoding.

    Params
    ------
    store
        The group to write to.
    k
        The key to write into the group.
        If the group is the root, set `k` to `"/"` to write directly into it.
        Passing an absolute path referring to a direct child of the group is also allowed.
    elem
        The element to write. Typically an in-memory object, e.g. an AnnData, pandas
        dataframe, scipy sparse matrix, etc.
    dataset_kwargs
        Keyword arguments to pass to the stores dataset creation function.
        E.g. for zarr this would be `chunks`, `compressor`.
    compat
        Which anndata version’s write behavior to target, see :class:`~anndata.WriteCompat`.
    """
    Writer(_REGISTRY, compat=compat).write_elem(
        store, k, elem, dataset_kwargs=dataset_kwargs
    )


# TODO: If all items would be read, just call normal read method
def read_elem_partial(
    elem,
    *,
    items=None,
    indices=(slice(None), slice(None)),
):
    """Read part of an element from an on disk store."""
    read_partial = _REGISTRY.get_partial_read(type(elem), get_spec(elem))
    return read_partial(elem, items=items, indices=indices)
