from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from functools import WRAPPER_ASSIGNMENTS, cache, wraps
from itertools import pairwise
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal, NamedTuple, cast

import numpy as np
import pandas as pd

from .._core.sparse_dataset import BaseCompressedSparseDataset
from ..utils import warn

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Iterable, Iterator
    from typing import Any, Literal

    import h5py
    import zarr
    from pandas.core.dtypes.dtypes import BaseMaskedDtype

    from .. import AnnData
    from .._core.raw import Raw
    from .._types import (
        StorageType,
        _ArrayStorageType,
        _GroupStorageType,
        _WriteInternal,
    )
    from .._write_compat import WriteConfig
    from ..acc import AdAcc, AdRef, DFAcc, MapAcc, RefAcc
    from ..typing import RWAble
    from .specs.registry import IOSpec, Writer

    Storage = StorageType | BaseCompressedSparseDataset


# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


# Could be numba’d if it returned tuples instead of slices
def idx_chunks_along_axis(shape: tuple[int, ...], axis: int, chunk_size: int):
    """\
    Gives indexer tuples chunked along an axis.

    Params
    ------
    shape
        Shape of array to be chunked
    axis
        Axis to chunk along
    chunk_size
        Size of chunk along axis

    Returns
    -------
    An iterator of tuples for indexing into an array of passed shape.
    """
    total = shape[axis]
    cur = 0
    mutable_idx = [slice(None) for i in range(len(shape))]
    while cur + chunk_size < total:
        mutable_idx[axis] = slice(cur, cur + chunk_size)
        yield tuple(mutable_idx)
        cur += chunk_size
    mutable_idx[axis] = slice(cur, None)
    yield tuple(mutable_idx)


def is_float(string):
    """\
    Check whether string is float.

    See also
    --------
    http://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string):
    """Check whether string is integer."""
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string):
    """Check whether string is boolean."""
    if string == "True":
        return True, True
    elif string == "False":
        return True, False
    else:
        return False, False


def convert_string(string):
    """Convert string to int, float or bool."""
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == "None":
        return None
    else:
        return string


def check_key(key):
    """Checks that passed value is a valid h5py key.

    Should convert it if there is an obvious conversion path, error otherwise.
    """
    typ = type(key)
    if issubclass(typ, str):
        return str(key)
    # TODO: Should I try to decode bytes? It's what h5py would do,
    # but it will be read out as a str.
    # elif issubclass(typ, bytes):
    # return key
    else:
        msg = f"{key} of type {typ} is an invalid key. Should be str."
        raise TypeError(msg)


@cache
def pandas_nullable_dtype(dtype: np.dtype) -> BaseMaskedDtype:
    """Infer nullable dtype from numpy dtype.

    There is no public pandas API for this, so this is the cleanest way.
    See <https://github.com/pandas-dev/pandas/issues/63608>
    """
    try:
        from pandas.core.dtypes.dtypes import BaseMaskedDtype
    except ImportError:
        pass
    else:
        if hasattr(BaseMaskedDtype, "from_numpy_dtype"):
            return BaseMaskedDtype.from_numpy_dtype(dtype)

    array_type: type[pd.arrays.BooleanArray | pd.arrays.IntegerArray]
    match dtype.kind:
        case "b":
            array_type = pd.arrays.BooleanArray
        case "i" | "u":
            array_type = pd.arrays.IntegerArray
        case _:
            raise NotImplementedError
    return array_type(np.ones(1, dtype), np.ones(1, bool)).dtype


def write_use(
    f: _GroupStorageType,
    adata: AnnData,
    use: AdRef | RefAcc | MapAcc | AdAcc | Sequence[AdRef | RefAcc | MapAcc | AdAcc],
    *,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    compat: WriteConfig,
) -> None:
    """Write only the parts of `adata` referenced by `use` into the existing anndata group `f`."""
    from .specs import write_elem

    # resolve everything first so an invalid accessor fails before anything is written
    uses = use if isinstance(use, Sequence) else [use]
    _check_store(f, adata, use)
    targets = [(u, *_descend(f, adata, u)) for u in uses]
    for _, g, inner in targets:
        if g is not f:
            _check_store(g, inner, use)
    # `g` is the nested group, but `u` resolves its own scope, so pass it the root `adata`
    for g, parent, key, elem in [
        (g, *_resolve_use(g, adata, u)) for u, g, _ in targets
    ]:
        if (
            parent and parent[0] not in g
        ):  # fixed attributes like `layers` may be missing
            write_elem(g, parent[0], {}, compat=compat)
        group = g
        for name in parent:  # `_resolve_use` checked that the rest exists
            group = cast("_GroupStorageType", dict(mapping_items(group, {name}))[name])
        if (layout := IndexedLayout.of(group)) is None:
            write_elem(group, key, elem, dataset_kwargs=dataset_kwargs, compat=compat)
            continue
        # the group stores names in an array, so add `key` there if it’s new
        names = layout.read_names(group)
        # skip the index, whose name might equal a column’s (e.g. `""`)
        i_idx = group.attrs.get("df:index")
        i = next(
            (j for j, n in enumerate(names) if n == key and j != i_idx), len(names)
        )
        write_elem(
            group,
            layout.child_key(i),
            elem,
            dataset_kwargs=dataset_kwargs,
            compat=compat,
        )
        if i == len(names):
            layout.write_names(group, [*names, key])
            if "df:shape" in group.attrs:
                n_rows, _ = cast("Sequence[int]", group.attrs["df:shape"])
                group.attrs["df:shape"] = [n_rows, i + 1]


def _check_store(f: _GroupStorageType, adata: AnnData, use: object) -> None:
    """Check that `f` is an anndata group with the same names as `adata`."""
    if (
        f.attrs.get("encoding-type"),
        f.attrs.get("encoding-version"),
    ) != ("anndata", "0.1.0"):
        msg = "Cannot write a subelement of an AnnData object to a store that is not an AnnData store."
        raise ValueError(msg)
    for dim in ("obs", "var"):
        on_disk = read_df_index(cast("_GroupStorageType", f[dim]))
        if isinstance(on_disk, pd.RangeIndex):  # AnnData turns these into strings
            on_disk = on_disk.astype(str)
        if not on_disk.equals(getattr(adata, f"{dim}_names")):
            msg = f"Cannot write {use}: {dim}_names in the store do not match the AnnData object."
            raise ValueError(msg)


def _descend(
    f: _GroupStorageType, adata: AnnData, use: AdRef | RefAcc | MapAcc | AdAcc
) -> tuple[_GroupStorageType, AnnData]:
    """Get the anndata group and object `use` refers into."""
    from .. import AnnData
    from ..acc import AdAcc, AdRef, ItemAcc, MapAcc, RefAcc

    match use:
        case AdRef(acc=RefAcc(scope=scope)) | RefAcc(scope=scope) | MapAcc(scope=scope):
            pass
        case AdAcc(scope=RefAcc(scope=scope)):  # same scope as `A.uns["n"]`
            pass
        case _:  # bare `A`, which `_resolve_use` rejects
            scope = None
    if scope is None:
        return f, adata
    f, _ = _descend(f, adata, scope)
    if not isinstance(scope, ItemAcc):  # e.g. a mudata modality
        msg = f"Cannot write {use}: {scope!r} has no known location in the store."
        raise ValueError(msg)
    group = _require_group(f, ("uns", *scope.path), "anndata", use)
    if not isinstance(inner := scope.get(adata), AnnData):
        msg = f"Cannot write {use}: {scope!r} is a {type(inner).__name__}, not an AnnData object."
        raise ValueError(msg)
    return group, inner


def _require_group(
    f: _GroupStorageType, path: Sequence[str], encoding: str, use: object
) -> _GroupStorageType:
    """Get the group at `path` in `f` if it has `encoding`, as only fixed attributes get created."""
    group = f
    for name in path:
        if (child := dict(mapping_items(group, {name})).get(name)) is None:
            break
        group = cast("_GroupStorageType", child)
    else:
        if group.attrs.get("encoding-type") == encoding:
            return group
    msg = f"Cannot write {use}: there is no {encoding!r} element at {'/'.join(path)!r} in the store. Write it first."
    raise ValueError(msg)


def _frame_col(
    group: _GroupStorageType,
    path: tuple[str, ...],
    adata: AnnData,
    acc: DFAcc,
    col: str,
) -> tuple[tuple[str, ...], str, RWAble]:
    """Get the parent path, key, and element to write column `col` of the frame `group` at `path`."""
    idx = df_children(group).index
    # in the 0.2.0 layout, `_index` is the child key of an unnamed index
    if (idx is not None and col == idx.name) or (
        col == "_index" and IndexedLayout.of(group) is None
    ):
        msg = f"Cannot write column {col!r}: it is the name of the {'/'.join(path)} index."
        raise ValueError(msg)
    # `._values` is what `write_dataframe` writes for each column
    return path, col, acc.frame(adata)[col]._values


def _resolve_use(  # noqa: PLR0911, PLR0912
    f: _GroupStorageType, adata: AnnData, use: AdRef | RefAcc | MapAcc | AdAcc
) -> tuple[tuple[str, ...], str, RWAble]:
    """Get the parent path (names, not child keys), key, and element to write for `use`."""
    from ..acc import (
        AdAcc,
        AdRef,
        DFAcc,
        GraphAcc,
        GraphMapAcc,
        ItemAcc,
        LayerAcc,
        LayerMapAcc,
        MetaAcc,
        MultiAcc,
        MultiMapAcc,
        RefAcc,
        UnsMapAcc,
    )

    match use:
        case DFAcc(frame_acc=RefAcc() as frame):  # `.DF`
            parent, key, _ = _resolve_use(f, adata, frame)
            return parent, key, use.get(adata)  # `.get` checks the cast
        case AdAcc(scope=RefAcc() as item):  # `.AD`
            parent, key, _ = _resolve_use(f, adata, item)
            return parent, key, use.container(adata)  # `.container` checks the cast
        case ItemAcc(path=path):
            for i in range(1, len(path)):  # `uns` itself is a fixed attribute
                _require_group(f, ("uns", *path[:i]), "dict", use)
            return ("uns", *path[:-1]), path[-1], use.get(adata)
        case UnsMapAcc():
            return (), "uns", dict(use.get(adata))
        case AdRef(acc=MetaAcc(dim=dim) as acc, idx=str() as col):
            return _frame_col(
                cast("_GroupStorageType", f[dim]), (dim,), adata, acc, col
            )
        case AdRef(acc=DFAcc(frame_acc=RefAcc() as frame) as acc, idx=str() as col):
            # `.DF`: write into the frame where it is in the store
            frame_parent, frame_key, _ = _resolve_use(f, adata, frame)
            path = (*frame_parent, frame_key)
            group = _require_group(f, path, "dataframe", use)
            return _frame_col(group, path, adata, acc, col)
        case AdRef():
            msg = f"Cannot write {use}: only data frame columns can be written as single vectors."
            raise ValueError(msg)
        case LayerAcc(k=k):
            return (
                () if k is None else ("layers",),
                "X" if k is None else k,
                use.get(adata),
            )
        case MetaAcc(dim=dim):
            return (), dim, use.get(adata)
        case MultiAcc(dim=dim, k=k):
            return (f"{dim}m",), k, use.get(adata)
        case GraphAcc(dim=dim, k=k):
            return (f"{dim}p",), k, use.get(adata)
        case LayerMapAcc():  # `adata.layers` contains `X` as `None`
            return (
                (),
                "layers",
                {k: v for k, v in use.get(adata).items() if k is not None},
            )
        case MultiMapAcc(dim=dim):
            return (), f"{dim}m", dict(use.get(adata))
        case GraphMapAcc(dim=dim):
            return (), f"{dim}p", dict(use.get(adata))
        case _:
            msg = f"Cannot write {use!r}."
            raise TypeError(msg)


def read_attribute(*args, **kwargs):
    from .specs import read_elem

    msg = "This internal function has been deprecated, please use read_elem instead"
    warn(msg, FutureWarning)
    return read_elem(*args, **kwargs)


def write_attribute(*args, **kwargs):
    from .specs import write_elem

    msg = "This internal function has been deprecated, please use write_elem instead"
    warn(msg, FutureWarning)
    return write_elem(*args, **kwargs)


# -------------------------------------------------------------------------------
# Errors handling
# -------------------------------------------------------------------------------
# TODO: Is there a consistent way to do this which just modifies the previously
# thrown error? Could do a warning?


class AnnDataReadError(OSError):
    """Error caused while trying to read in AnnData."""


def _get_display_path(store: Storage) -> str:
    """Return an absolute path of an element (always starts with “/”)."""
    group = store.group if isinstance(store, BaseCompressedSparseDataset) else store
    path = group.name or "??"  # can be None
    return f"/{path.removeprefix('/')}"


def add_key_note(
    e: BaseException, store: Storage, path: str, key: str, op: Literal["read", "writ"]
) -> None:
    if any(
        f"Error raised while {op}ing key" in note
        for note in getattr(e, "__notes__", [])
    ):
        return

    dir = "to" if op == "writ" else "from"
    msg = f"Error raised while {op}ing key {key!r} of {type(store)} {dir} {path}"
    e.add_note(msg)


def report_read_key_on_error(func):
    """\
    A decorator for hdf5/zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> import numpy as np
    >>> @report_read_key_on_error
    ... def read_arr(group):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr", mode="w")
    >>> z["X"] = np.array([1, 2, 3])
    >>> read_arr(z["X"])  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        from anndata._io.specs import Reader

        # Figure out signature (method vs function) by going through args
        for arg in args:
            if not isinstance(arg, Reader):
                store = cast("Storage", arg)
                break
        else:
            msg = "No element found in args."
            raise ValueError(msg)
        try:
            return func(*args, **kwargs)
        except Exception as e:
            path, key = _get_display_path(store).rsplit("/", 1)
            add_key_note(e, store, path or "/", key, "read")
            raise

    return func_wrapper


def report_write_key_on_error(func):
    """\
    A decorator for hdf5/zarr element writing which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_write_key_on_error
    ... def write_arr(group, key, val):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr", mode="w")
    >>> X = [1, 2, 3]
    >>> write_arr(z, "X", X)  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        __tracebackhide__ = True
        from anndata._io.specs import Writer

        # Figure out signature (method vs function) by going through args
        for arg, _key in pairwise(args):
            key = _key
            if not isinstance(arg, Writer):
                store = cast("Storage", arg)
                break
        else:
            msg = "No element found in args."
            raise ValueError(msg)
        try:
            return func(*args, **kwargs)
        except Exception as e:
            path = _get_display_path(store)
            add_key_note(e, store, path, key, "writ")
            raise

    return func_wrapper


def _check_has_no_slash_key(attr: str, elem: object, *, compat: WriteConfig) -> None:
    """Only attempt to write slash keys where people rely on it for backwards compatibility."""

    if attr in {"obs", "var", "uns", "raw"}:
        return  # `write_elem` checks these against `settings.disallow_forward_slash_in_h5ad`
    if compat.indexed_keys:
        return  # keys are stored in an array, not as child names
    assert isinstance(elem, Mapping)
    if any("/" in k for k in elem if k not in {"/", None}):
        msg = f"Forward slashes are not allowed in keys in {attr}"
        raise ValueError(msg)


@dataclass(frozen=True)
class _IndexedLayout:
    encoding_type: str
    encoding_version: str
    names_key: str
    """Child key of the names array."""
    names_attr: str
    """Attribute holding the names instead of the names array."""
    prefix: str
    """Child keys are `{prefix}_{i}`."""


NAMES_ATTR_MAX_BYTES = 4096
"""Larger names arrays don’t go into attributes, which are loaded whenever the group is opened."""


class IndexedLayout(_IndexedLayout, Enum):
    """A layout storing child names in an attribute or array and naming children by position."""

    DATAFRAME = ("dataframe", "0.3.0", "columns", "df:columns", "column")
    # no convention defines `dict:keys`, it mirrors the dataframe convention’s `df:columns`
    DICT = ("dict", "0.2.0", "keys", "dict:keys", "item")

    @classmethod
    def of(cls, group: _GroupStorageType) -> IndexedLayout | None:
        """Get the layout `group` uses, if any."""
        spec = (group.attrs.get("encoding-type"), group.attrs.get("encoding-version"))
        return next(
            (l for l in cls if (l.encoding_type, l.encoding_version) == spec), None
        )

    @property
    def spec(self) -> IOSpec:
        from .specs.registry import IOSpec

        return IOSpec(self.encoding_type, self.encoding_version)

    def child_key(self, i: int) -> str:
        return f"{self.prefix}_{i}"

    def read_names(self, group: _GroupStorageType) -> list[str]:
        from .specs import read_elem

        names = group.attrs.get(self.names_attr)
        if names is None:
            names = read_elem(group[self.names_key])
        elif self.names_key in group:
            msg = (
                f"Group {group.name!r} stores names both in attribute "
                f"{self.names_attr!r} and array {self.names_key!r}"
            )
            raise ValueError(msg)
        return [str(n) for n in cast("Iterable[str]", names)]

    def write_names(self, group: _GroupStorageType, names: Sequence[str]) -> None:
        """Store `names` in an attribute (saving a request) if small enough, else in an array.

        HDF5 always gets the array, as it has no requests to save but a limit on attribute size.
        """
        import zarr

        from .specs import write_elem

        names = list(names)
        if isinstance(group, zarr.Group) and (
            len(json.dumps(names)) <= NAMES_ATTR_MAX_BYTES
        ):
            if self.names_key in group:
                del group[self.names_key]
            group.attrs[self.names_attr] = names
            return
        if self.names_attr in group.attrs:
            del group.attrs[self.names_attr]
        # object dtype makes this a `string-array`, not a `nullable-string-array` group
        write_elem(group, self.names_key, np.array(names, dtype=object))


def write_indexed(
    group: _GroupStorageType,
    layout: IndexedLayout,
    items: Iterable[tuple[str, RWAble]],
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any],
) -> int:
    """Write `items` as numbered children and their names via :meth:`IndexedLayout.write_names`.

    Returns the number of items written.
    """
    names: list[str] = []
    for i, (name, value) in enumerate(items):
        _writer.write_elem(
            group, layout.child_key(i), value, dataset_kwargs=dataset_kwargs
        )
        names.append(name)
    layout.write_names(group, names)
    return len(names)


def mapping_items(
    group: _GroupStorageType, keys: Collection[str] | None = None
) -> Iterator[tuple[str, StorageType]]:
    """Iterate over the keys and children of a `dict` encoded group.

    If `keys` is given, only open the children with these keys that exist.
    """
    if (layout := IndexedLayout.of(group)) is None:
        if keys is None:
            return iter(dict(group).items())
        return ((k, group[k]) for k in group if k in keys)  # store order, like below
    names = layout.read_names(group)
    return (
        (k, group[layout.child_key(i)])
        for i, k in enumerate(names)
        if keys is None or k in keys
    )


class DfChild[N: str | None](NamedTuple):
    """A dataframe’s index or column."""

    name: N
    """Its name, which is `None` for an unnamed index."""
    key: str
    """The key of the child holding its data."""


@dataclass(frozen=True)
class DfChildren:
    index: DfChild[str | None] | None
    """The index, which is optional from v0.3.0 on."""
    columns: list[DfChild[str]]


def df_children(group: _GroupStorageType) -> DfChildren:
    """Get the index and columns of a `dataframe` encoded group."""
    from ..compat import _read_attr

    if (layout := IndexedLayout.of(group)) is None:
        idx = cast("str", _read_attr(group.attrs, "_index"))
        col_order = cast("Iterable[str]", _read_attr(group.attrs, "column-order"))
        return DfChildren(
            DfChild(None if idx == "_index" else idx, idx),
            [DfChild(c, c) for c in col_order],
        )
    names = layout.read_names(group)
    cols = [DfChild(n, layout.child_key(i)) for i, n in enumerate(names)]
    if (i_idx := group.attrs.get("df:index")) is None:
        return DfChildren(None, cols)
    idx_name, idx_key = cols.pop(int(cast("int", i_idx)))
    return DfChildren(DfChild(idx_name or None, idx_key), cols)


def read_df_index(group: _GroupStorageType) -> pd.Index:
    """Read a `dataframe` encoded group’s index, or a `RangeIndex` if it has none."""
    from .specs import read_elem

    if (idx := df_children(group).index) is None:
        n_rows, _ = cast("Sequence[int]", group.attrs["df:shape"])
        return pd.RangeIndex(int(n_rows))
    return pd.Index(cast("np.ndarray", read_elem(cast("StorageType", group[idx.key]))))


# -------------------------------------------------------------------------------
# Common h5ad/zarr stuff
# -------------------------------------------------------------------------------


def _read_legacy_raw[M: (dict[str, Any], Raw)](
    f: zarr.Group | h5py.Group,
    modern_raw: M | None,
    read_df: Callable,
    read_attr: Callable,
    *,
    attrs: Collection[str] = ("X", "var", "varm"),
) -> M | dict[str, Any]:
    """\
    Backwards compat for reading legacy raw.
    Makes sure that no modern raw group coexists with legacy raw.* groups.
    """
    if modern_raw:
        if any(k.startswith("raw.") for k in f):
            what = f"File {f.filename}" if hasattr(f, "filename") else "Store"
            msg = f"{what} has both legacy and current raw formats."
            raise ValueError(msg)
        return modern_raw

    raw = {}
    if "X" in attrs and "raw.X" in f:
        raw["X"] = read_attr(f["raw.X"])
    if "var" in attrs and "raw.var" in f:
        raw["var"] = read_df(f["raw.var"])  # Backwards compat
    if "varm" in attrs and "raw.varm" in f:
        raw["varm"] = read_attr(f["raw.varm"])
    return raw


def zero_dim_array_as_scalar[S: StorageType, T: np.ndarray | _ArrayStorageType](
    func: _WriteInternal[S, T],
) -> _WriteInternal[S, T]:
    """\
    A decorator for write_elem implementations of arrays where zero-dimensional arrays need special handling.
    """

    @wraps(func, assigned=(*WRAPPER_ASSIGNMENTS, "__defaults__", "__kwdefaults__"))
    def func_wrapper(
        f: S,
        k: str,
        elem: T,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None:
        if elem.shape == ():
            _writer.write_elem(f, k, elem[()], dataset_kwargs=dataset_kwargs)
        else:
            func(f, k, elem, _writer=_writer, dataset_kwargs=dataset_kwargs)

    return func_wrapper


def no_write_dataset_2d(write):
    def raise_error_if_dataset_2d_present(store, adata, *args, **kwargs):
        from anndata.experimental.backed._compat import has_dataset_2d

        if has_dataset_2d(adata):
            msg = (
                "Writing AnnData objects with a Dataset2D not supported yet. "
                "Please use `ds.to_memory` to bring the dataset into memory. "
                "Note that if you have generated this object by concatenating several `AnnData` objects"
                "the original types may be lost."
            )
            raise NotImplementedError(msg)
        return write(store, adata, *args, **kwargs)

    return raise_error_if_dataset_2d_present
