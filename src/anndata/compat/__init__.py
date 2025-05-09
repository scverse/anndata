from __future__ import annotations

from codecs import decode
from collections.abc import Mapping
from functools import cache, partial, singledispatch
from importlib.util import find_spec
from types import EllipsisType
from typing import TYPE_CHECKING, TypeVar
from warnings import warn

import h5py
import numpy as np
import pandas as pd
import scipy
from packaging.version import Version
from zarr import Array as ZarrArray  # noqa: F401
from zarr import Group as ZarrGroup

if TYPE_CHECKING:
    from typing import Any

#############################
# scipy sparse array comapt #
#############################


CSMatrix = scipy.sparse.csr_matrix | scipy.sparse.csc_matrix
CSArray = scipy.sparse.csr_array | scipy.sparse.csc_array


class Empty:
    pass


Index1D = slice | int | str | np.int64 | np.ndarray | pd.Series
IndexRest = Index1D | EllipsisType
Index = (
    IndexRest
    | tuple[Index1D, IndexRest]
    | tuple[IndexRest, Index1D]
    | tuple[Index1D, Index1D, EllipsisType]
    | tuple[EllipsisType, Index1D, Index1D]
    | tuple[Index1D, EllipsisType, Index1D]
    | CSMatrix
    | CSArray
)
H5Group = h5py.Group
H5Array = h5py.Dataset
H5File = h5py.File


#############################
# Optional deps
#############################
@cache
def is_zarr_v2() -> bool:
    import zarr
    from packaging.version import Version

    return Version(zarr.__version__) < Version("3.0.0")


if is_zarr_v2():
    msg = "anndata will no longer support zarr v2 in the near future. Please prepare to upgrade to zarr>=3."
    warn(msg, DeprecationWarning, stacklevel=2)


if find_spec("awkward") or TYPE_CHECKING:
    import awkward  # noqa: F401
    from awkward import Array as AwkArray
else:

    class AwkArray:
        @staticmethod
        def __repr__():
            return "mock awkward.highlevel.Array"


if find_spec("zappy") or TYPE_CHECKING:
    from zappy.base import ZappyArray
else:

    class ZappyArray:
        @staticmethod
        def __repr__():
            return "mock zappy.base.ZappyArray"


if TYPE_CHECKING:
    # type checkers are confused and can only see …core.Array
    from dask.array.core import Array as DaskArray
elif find_spec("dask"):
    from dask.array import Array as DaskArray
else:

    class DaskArray:
        @staticmethod
        def __repr__():
            return "mock dask.array.core.Array"


if find_spec("xarray") or TYPE_CHECKING:
    import xarray
    from xarray import DataArray as XDataArray
    from xarray import Dataset as XDataset
    from xarray.backends import BackendArray as XBackendArray
    from xarray.backends.zarr import ZarrArrayWrapper as XZarrArrayWrapper
else:
    xarray = None

    class XDataArray:
        def __repr__(self) -> str:
            return "mock DataArray"

    class XDataset:
        def __repr__(self) -> str:
            return "mock Dataset"

    class XZarrArrayWrapper:
        def __repr__(self) -> str:
            return "mock ZarrArrayWrapper"

    class XBackendArray:
        def __repr__(self) -> str:
            return "mock BackendArray"


# https://github.com/scverse/anndata/issues/1749
def is_cupy_importable() -> bool:
    try:
        import cupy  # noqa: F401
    except ImportError:
        return False
    return True


if is_cupy_importable() or TYPE_CHECKING:
    from cupy import ndarray as CupyArray
    from cupyx.scipy.sparse import csc_matrix as CupyCSCMatrix
    from cupyx.scipy.sparse import csr_matrix as CupyCSRMatrix
    from cupyx.scipy.sparse import spmatrix as CupySparseMatrix

    try:
        import dask.array as da
    except ImportError:
        pass
    else:
        da.register_chunk_type(CupyCSRMatrix)
        da.register_chunk_type(CupyCSCMatrix)
else:

    class CupySparseMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.spmatrix"

    class CupyCSRMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.csr_matrix"

    class CupyCSCMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.csc_matrix"

    class CupyArray:
        @staticmethod
        def __repr__():
            return "mock cupy.ndarray"


if find_spec("legacy_api_wrap") or TYPE_CHECKING:
    from legacy_api_wrap import legacy_api  # noqa: TID251

    old_positionals = partial(legacy_api, category=FutureWarning)
else:

    def old_positionals(*old_positionals):
        return lambda func: func


#############################
# IO helpers
#############################


@singledispatch
def _read_attr(attrs: Mapping, name: str, default: Any | None = Empty):
    if default is Empty:
        return attrs[name]
    else:
        return attrs.get(name, default=default)


@_read_attr.register(h5py.AttributeManager)
def _read_attr_hdf5(
    attrs: h5py.AttributeManager, name: str, default: Any | None = Empty
):
    """
    Read an HDF5 attribute and perform all necessary conversions.

    At the moment, this only implements conversions for string attributes, other types
    are passed through. String conversion is needed compatibility with other languages.
    For example Julia's HDF5.jl writes string attributes as fixed-size strings, which
    are read as bytes by h5py.
    """
    if name not in attrs and default is not Empty:
        return default
    attr = attrs[name]
    attr_id = attrs.get_id(name)
    dtype = h5py.check_string_dtype(attr_id.dtype)
    if dtype is None or dtype.length is None:  # variable-length string, no problem
        return attr
    elif len(attr_id.shape) == 0:  # Python bytestring
        return attr.decode("utf-8")
    else:  # NumPy array
        return [decode(s, "utf-8") for s in attr]


def _from_fixed_length_strings(value):
    """\
    Convert from fixed length strings to unicode.

    For backwards compatibility with older h5ad and zarr files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        # could probably match better
        is_annotated = isinstance(dt_type, tuple)
        if is_annotated:
            dt_type = dt_type[0]
        # Fixing issue introduced with h5py v2.10.0, see:
        # https://github.com/h5py/h5py/issues/1307
        if issubclass(np.dtype(dt_type).type, np.bytes_):
            dt_list[1] = f"U{int(dt_type[2:])}"
        elif is_annotated or np.issubdtype(np.dtype(dt_type), np.str_):
            dt_list[1] = "O"  # Assumption that it’s a vlen str
        new_dtype.append(tuple(dt_list))
    return value.astype(new_dtype)


def _decode_structured_array(
    arr: np.ndarray, *, dtype: np.dtype | None = None, copy: bool = False
) -> np.ndarray:
    """
    h5py 3.0 now reads all strings as bytes. There is a helper method which can convert these to strings,
    but there isn't anything for fields of structured dtypes.

    Params
    ------
    arr
        An array with structured dtype
    dtype
        dtype of the array. This is checked for h5py string data types.
        Passing this is allowed for cases where array may have been processed by another function before hand.
    """
    if copy:
        arr = arr.copy()
    if dtype is None:
        dtype = arr.dtype
    # codecs.decode is 2x slower than this lambda, go figure
    decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
    for k, (dt, _) in dtype.fields.items():
        check = h5py.check_string_dtype(dt)
        if check is not None and check.encoding == "utf-8":
            decode(arr[k], out=arr[k])
    return arr


def _to_fixed_length_strings(value: np.ndarray) -> np.ndarray:
    """\
    Convert variable length strings to fixed length.

    Currently a workaround for
    https://github.com/zarr-developers/zarr-python/pull/422
    """
    new_dtype = []
    for dt_name, (dt_type, dt_offset) in value.dtype.fields.items():
        if dt_type.kind == "O":
            #  Assuming the objects are str
            size = max(len(x.encode()) for x in value.getfield("O", dt_offset))
            new_dtype.append((dt_name, ("U", size)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


Group_T = TypeVar("Group_T", bound=ZarrGroup | h5py.Group)


# TODO: This is a workaround for https://github.com/scverse/anndata/issues/874
# See https://github.com/h5py/h5py/pull/2311#issuecomment-1734102238 for why this is done this way.
def _require_group_write_dataframe(
    f: Group_T, name: str, df: pd.DataFrame, *args, **kwargs
) -> Group_T:
    if len(df.columns) > 5_000 and isinstance(f, H5Group):
        # actually 64kb is the limit, but this should be a conservative estimate
        return f.create_group(name, *args, track_order=True, **kwargs)
    return f.require_group(name, *args, **kwargs)


#############################
# Dealing with uns
#############################


def _clean_uns(adata: AnnData):  # noqa: F821
    """
    Compat function for when categorical keys were stored in uns.
    This used to be buggy because when storing categorical columns in obs and var with
    the same column name, only one `<colname>_categories` is retained.
    """
    k_to_delete = set()
    for cats_name, cats in adata.uns.items():
        if not cats_name.endswith("_categories"):
            continue
        name = cats_name.replace("_categories", "")
        # fix categories with a single category
        if isinstance(cats, str | int):
            cats = [cats]
        for ann in [adata.obs, adata.var]:
            if name not in ann:
                continue
            codes: np.ndarray = ann[name].values
            # hack to maybe find the axis the categories were for
            if not np.all(codes < len(cats)):
                continue
            ann[name] = pd.Categorical.from_codes(codes, cats)
            k_to_delete.add(cats_name)
    for cats_name in k_to_delete:
        del adata.uns[cats_name]


def _move_adj_mtx(d):
    """
    Read-time fix for moving adjacency matrices from uns to obsp
    """
    n = d.get("uns", {}).get("neighbors", {})
    obsp = d.setdefault("obsp", {})

    for k in ("distances", "connectivities"):
        if (
            (k in n)
            and isinstance(n[k], scipy.sparse.spmatrix | np.ndarray)
            and len(n[k].shape) == 2
        ):
            msg = (
                f"Moving element from .uns['neighbors'][{k!r}] to .obsp[{k!r}].\n\n"
                "This is where adjacency matrices should go now."
            )
            # 5: caller -> 4: legacy_api_wrap -> 3: `AnnData.__init__` -> 2: `_init_as_actual` → 1: here
            warn(msg, FutureWarning, stacklevel=5)
            obsp[k] = n.pop(k)


def _find_sparse_matrices(d: Mapping, n: int, keys: tuple, paths: list):
    """Find paths to sparse matrices with shape (n, n)."""
    for k, v in d.items():
        if isinstance(v, Mapping):
            _find_sparse_matrices(v, n, (*keys, k), paths)
        elif scipy.sparse.issparse(v) and v.shape == (n, n):
            paths.append((*keys, k))
    return paths


def _transpose_by_block(dask_array: DaskArray) -> DaskArray:
    import dask.array as da

    b = dask_array.blocks
    b_raveled = b.ravel()
    block_layout = np.zeros(b.shape, dtype=object)

    for i in range(block_layout.size):
        block_layout.flat[i] = b_raveled[i].map_blocks(
            lambda x: x.T, chunks=b_raveled[i].chunks[::-1]
        )

    return da.block(block_layout.T.tolist())


def _safe_transpose(x):
    """Safely transpose x

    This is a workaround for: https://github.com/scipy/scipy/issues/19161
    """

    if isinstance(x, DaskArray) and scipy.sparse.issparse(x._meta):
        return _transpose_by_block(x)
    else:
        return x.T


def _map_cat_to_str(cat: pd.Categorical) -> pd.Categorical:
    if Version(pd.__version__) >= Version("2.1"):
        # Argument added in pandas 2.1
        return cat.map(str, na_action="ignore")
    else:
        return cat.map(str)
