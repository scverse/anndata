from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from functools import WRAPPER_ASSIGNMENTS, cache, wraps
from itertools import pairwise
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
import pandas as pd

from .._core.sparse_dataset import BaseCompressedSparseDataset
from ..utils import warn

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    import h5py
    import zarr
    from pandas.core.dtypes.dtypes import BaseMaskedDtype

    from .. import AnnData
    from .._types import (
        ReadCallback,
        StorageType,
        _ArrayStorageType,
        _GroupStorageType,
        _WriteInternal,
    )
    from ..acc import AdRef, MapAcc, RefAcc
    from ..typing import RWAble
    from .specs.registry import Writer

    Storage = StorageType | BaseCompressedSparseDataset


# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


# Could be numba’d if it returned tuples instead of slices
def idx_chunks_along_axis(shape: tuple, axis: int, chunk_size: int):
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


def write_subelems(
    f: _GroupStorageType,
    adata: AnnData,
    elem: AdRef | RefAcc | MapAcc | Sequence[AdRef | RefAcc | MapAcc],
    *,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> None:
    """Write only the parts of `adata` referenced by `elem` into the existing anndata group `f`."""
    from .specs import read_elem, write_elem

    if (
        f.attrs.get("encoding-type"),
        f.attrs.get("encoding-version"),
    ) != ("anndata", "0.1.0"):
        msg = "Cannot write a subelement of an AnnData object to a store that is not an AnnData store."
        raise ValueError(msg)
    for dim in ("obs", "var"):
        group = cast("_GroupStorageType", f[dim])
        on_disk = read_elem(
            cast("StorageType", group[cast("str", group.attrs["_index"])])
        )
        if not pd.Index(cast("np.ndarray", on_disk)).equals(
            getattr(adata, f"{dim}_names")
        ):
            msg = f"Cannot write {elem}: {dim}_names in the store do not match the AnnData object."
            raise ValueError(msg)

    # resolve everything first so an invalid accessor fails before anything is written
    elems = elem if isinstance(elem, Sequence) else [elem]
    for parent, key, value in [_resolve_subelem(f, adata, e) for e in elems]:
        if parent is not None and parent not in f:
            write_elem(f, parent, {})
        write_elem(
            f if parent is None else cast("_GroupStorageType", f[parent]),
            key,
            value,
            dataset_kwargs=dataset_kwargs,
        )


def _resolve_subelem(  # noqa: PLR0911
    f: _GroupStorageType, adata: AnnData, elem: AdRef | RefAcc | MapAcc
) -> tuple[str | None, str, RWAble]:
    """Get the parent group, key, and element to write for `elem`."""
    from ..acc import (
        AdRef,
        GraphAcc,
        GraphMapAcc,
        LayerAcc,
        LayerMapAcc,
        MetaAcc,
        MultiAcc,
        MultiMapAcc,
    )

    match elem:
        case AdRef(acc=MetaAcc(dim=dim) as acc, idx=str() as col):
            if col in {"_index", cast("_GroupStorageType", f[dim]).attrs["_index"]}:
                msg = f"Cannot write column {col!r}: it is the name of the {dim} index."
                raise ValueError(msg)
            # `._values` is what `write_dataframe` writes for each column
            return dim, col, acc.get(adata)[col]._values
        case AdRef():
            msg = f"Cannot write {elem}: only `obs`/`var` columns can be written as single vectors."
            raise ValueError(msg)
        case LayerAcc(k=k):
            return (
                None if k is None else "layers",
                "X" if k is None else k,
                elem.get(adata),
            )
        case MetaAcc(dim=dim):
            return None, dim, elem.get(adata)
        case MultiAcc(dim=dim, k=k):
            return f"{dim}m", k, elem.get(adata)
        case GraphAcc(dim=dim, k=k):
            return f"{dim}p", k, elem.get(adata)
        case LayerMapAcc():  # `adata.layers` contains `X` as `None`
            return (
                None,
                "layers",
                {k: v for k, v in elem.get(adata).items() if k is not None},
            )
        case MultiMapAcc(dim=dim):
            return None, f"{dim}m", dict(elem.get(adata))
        case GraphMapAcc(dim=dim):
            return None, f"{dim}p", dict(elem.get(adata))
        case _:
            msg = f"Cannot write {elem!r}."
            raise TypeError(msg)


def filter_callback(
    f: _GroupStorageType,
    filter: AdRef | Sequence[AdRef],
    callback: ReadCallback,
) -> ReadCallback:
    """Wrap `callback` to only read the obs/var subset selected by the boolean `obs`/`var` columns in `filter`."""
    from .._settings import settings
    from ..acc import AdRef, MetaAcc
    from .specs import read_elem
    from .specs.registry import read_elem_partial

    masks: dict[str, np.ndarray] = {}
    for ref in [filter] if isinstance(filter, AdRef) else filter:
        match ref:
            case AdRef(acc=MetaAcc(dim=dim), idx=str() as col):
                mask = read_elem(cast("_GroupStorageType", f[dim])[col])
            case _:
                msg = (
                    f"Cannot filter by {ref!r}: only `obs`/`var` columns are supported."
                )
                raise ValueError(msg)
        if not pd.api.types.is_bool_dtype(mask):
            msg = f"Cannot filter by {ref!r}: it is not a boolean column."
            raise ValueError(msg)
        masks[dim] = masks.get(dim, True) & np.asarray(mask, dtype=bool)
    obs, var = (
        np.flatnonzero(masks[dim]) if dim in masks else slice(None)
        for dim in ("obs", "var")
    )

    def cb(read_func, elem_name: str, elem: StorageType, iospec):
        match elem_name.strip("/").split("/"):
            case ["X"] | ["layers", _]:
                idx = (obs, var)
            case ["obs"] | ["obsm", _] | ["raw", "X"]:
                idx = (obs, slice(None))
            case ["var"] | ["varm", _]:
                idx = (var, slice(None))
            case ["obsp", _]:
                idx = (obs, obs)
            case ["varp", _]:
                idx = (var, var)
            case [""]:
                adata = callback(read_func, elem_name, elem, iospec=iospec)
                if settings.remove_unused_categories:  # like `adata[...]` does
                    for df in (adata.obs, adata.var):
                        type(adata)._remove_unused_categories(df, df, adata.uns)
                    for df in (*adata.obsm.values(), *adata.varm.values()):
                        if isinstance(df, pd.DataFrame):
                            type(adata)._remove_unused_categories(df, df, {})
                return adata
            case _:
                return callback(read_func, elem_name, elem, iospec=iospec)
        return read_elem_partial(elem, indices=idx)

    return cb


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


def _check_has_no_slash_key(attr: str, elem: object) -> None:
    """Only attempt to write slash keys where people rely on it for backwards compatibility."""
    if attr in {"obs", "var", "uns", "raw"}:
        return  # separate check for `settings.disallow_forward_slash_in_h5ad` is done in `write_elem`
    assert isinstance(elem, Mapping)
    if any("/" in k for k in elem if k not in {"/", None}):
        msg = f"Forward slashes are not allowed in keys in {attr}"
        raise ValueError(msg)


# -------------------------------------------------------------------------------
# Common h5ad/zarr stuff
# -------------------------------------------------------------------------------


def _read_legacy_raw(
    f: zarr.Group | h5py.Group,
    modern_raw,  # TODO: type
    read_df: Callable,
    read_attr: Callable,
    *,
    attrs=("X", "var", "varm"),
) -> dict:
    """\
    Backwards compat for reading legacy raw.
    Makes sure that no modern raw group coexists with legacy raw.* groups.
    """
    # a 0-obs `AnnData` is falsy, so only treat empty mappings as missing
    if modern_raw is not None and (not isinstance(modern_raw, Mapping) or modern_raw):
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
