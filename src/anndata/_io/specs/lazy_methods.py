from __future__ import annotations

from contextlib import contextmanager
from functools import partial, singledispatch
from pathlib import Path
from typing import TYPE_CHECKING, overload

import h5py
import numpy as np
import pandas as pd
import zarr
from scipy import sparse

import anndata as ad
from anndata._core.file_backing import filename, get_elem_name
from anndata._core.xarray import Dataset2D, requires_xarray
from anndata.compat import DaskArray, XDataset, XVariable, pandas_as_str

from .registry import _LAZY_REGISTRY, IOSpec, read_elem

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping, Sequence
    from typing import Any, Final, Literal

    from anndata.abc import CSCDataset, CSRDataset
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ..._types import StorageType
    from ...compat import CSArray, CSMatrix
    from .registry import LazyDataStructures, LazyReader

    BlockInfo = Mapping[
        None,
        dict[str, Sequence[tuple[int, int]]],
    ]


@overload
@contextmanager
def maybe_open_h5(
    path: Path | h5py.File, /, elem_name: str
) -> Generator[h5py.File]: ...
@overload
@contextmanager  # D actually accepts anything, but that’d confuse the type checker
def maybe_open_h5[D: (zarr.Group, CSRDataset, CSCDataset)](
    obj: D, /, elem_name: str
) -> Generator[D]: ...
@contextmanager
def maybe_open_h5(
    path_or_other: Path | h5py.File | zarr.Group | CSRDataset | CSCDataset,
    /,
    elem_name: str,
) -> Generator[h5py.File | zarr.Group | CSRDataset | CSCDataset]:
    if not isinstance(path_or_other, Path):
        yield path_or_other
        return
    file = h5py.File(path_or_other, "r")
    try:
        yield file[elem_name]
    finally:
        file.close()


_DEFAULT_STRIDE = 1000


def compute_chunk_layout_for_axis_size(
    chunk_axis_size: int, full_axis_size: int
) -> tuple[int, ...]:
    # A zero-length axis must be expressed as a single zero-length chunk:
    # dask rejects empty chunk tuples ("Empty tuples are not allowed in
    # chunks"), and ``np.divmod(0, 0)`` would otherwise divide by zero.
    if full_axis_size == 0:
        return (0,)
    n_strides, rest = np.divmod(full_axis_size, chunk_axis_size)
    chunk = (chunk_axis_size,) * n_strides
    if rest > 0:
        chunk += (rest,)
    return chunk


def make_dask_chunk(
    path_or_sparse_dataset: Path | CSRDataset | CSCDataset,
    elem_name: str,
    block_info: BlockInfo | None = None,
) -> CSMatrix | CSArray | np.typing.NDArray:
    if block_info is None:
        msg = "Block info is required"
        raise ValueError(msg)
    # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
    # https://github.com/scverse/anndata/issues/1105
    with maybe_open_h5(path_or_sparse_dataset, elem_name) as f:
        return _compute_chunk(f, block_info)


def _compute_chunk(
    f: h5py.File | CSRDataset | CSCDataset, block_info: BlockInfo
) -> CSMatrix | CSArray | np.typing.NDArray:
    # See https://github.com/scverse/anndata/pull/2005 for why
    # should_cache_indptr is False.
    # The purpose of caching the indptr was when the dataset is reused
    # which is in general the case but is not here.  Hence
    # caching it on every access to the dataset here is quite costly.
    # `maybe_open_h5` yields a dense `Dataset`, a sparse `Group` or a sparse dataset
    mtx = (
        ad.io.sparse_dataset(f, should_cache_indptr=False)
        if isinstance(f, h5py.Group)
        else f
    )
    idx = tuple(
        slice(start, stop) for start, stop in block_info[None]["array-location"]
    )
    rv = mtx[idx]
    if TYPE_CHECKING:  # annotation bug: indexing with slices never returns a scalar
        assert not isinstance(rv, int | float)
    return rv


@singledispatch
def get_chunksize(obj) -> tuple[int, ...]:
    if hasattr(obj, "chunks"):
        return obj.chunks
    msg = "object of type {type(obj)} has no recognized chunks"
    raise ValueError(msg)


@_LAZY_REGISTRY.register_read(h5py.Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(h5py.Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(zarr.Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(zarr.Group, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: h5py.Group | zarr.Group,
    *,
    _reader: LazyReader,
    # the reader registry fixes this signature; only `tuple[int, int]` is accepted
    chunks: tuple[int | None, ...] | None = None,
) -> DaskArray:
    import dask.array as da

    path_or_sparse_dataset: Path | CSRDataset | CSCDataset
    if isinstance(elem, h5py.Group):
        path_or_sparse_dataset = Path(filename(elem))
        dtype = elem["data"].dtype
    else:
        path_or_sparse_dataset = ad.io.sparse_dataset(elem, should_cache_indptr=False)
        dtype = path_or_sparse_dataset.dtype
    elem_name = get_elem_name(elem)
    attrs: Mapping[str, Any] = elem.attrs
    shape: tuple[int, int] = tuple(attrs["shape"])
    is_csc: bool = attrs["encoding-type"] == "csc_matrix"

    stride: int = _DEFAULT_STRIDE
    major_dim, minor_dim = (1, 0) if is_csc else (0, 1)
    if chunks is not None:
        if len(chunks) != 2:
            msg = "`chunks` must be a tuple of two integers"
            raise ValueError(msg)
        if chunks[minor_dim] not in {shape[minor_dim], -1, None}:
            msg = (
                "Only the major axis can be chunked. "
                f"Try setting chunks to {((-1, _DEFAULT_STRIDE) if is_csc else (_DEFAULT_STRIDE, -1))}"
            )
            raise ValueError(msg)
        major_chunk = chunks[major_dim]
        stride = major_chunk if major_chunk not in {None, -1} else shape[major_dim]

    shape_minor, shape_major = shape if is_csc else shape[::-1]
    chunks_major = compute_chunk_layout_for_axis_size(stride, shape_major)
    chunks_minor = (shape_minor,)
    chunk_layout = (
        (chunks_minor, chunks_major) if is_csc else (chunks_major, chunks_minor)
    )
    memory_format = sparse.csc_matrix if is_csc else sparse.csr_matrix
    make_chunk = partial(make_dask_chunk, path_or_sparse_dataset, elem_name)
    da_mtx = da.map_blocks(
        make_chunk,
        dtype=dtype,
        chunks=chunk_layout,
        meta=memory_format((0, 0), dtype=dtype),
    )
    return da_mtx


def resolve_chunks(
    elem: h5py.Dataset | zarr.Array, chunks_arg: tuple[int | None, ...] | None
) -> tuple[int, ...]:
    shape = tuple(elem.shape)
    if chunks_arg is not None:
        # None and -1 on a given axis indicate that one should use the shape
        # in `dask`'s semantics.
        return tuple(
            c if c not in {None, -1} else s
            for c, s in zip(chunks_arg, shape, strict=True)
        )
    elif elem.chunks is None:  # h5 unchunked
        return tuple(min(_DEFAULT_STRIDE, s) for s in shape)
    return elem.chunks


# TODO: `map_blocks` of a string array in h5py is so insanely slow on benchmarking that in the case someone has
# a pure string annotation (not categoricals! or nullables strings!), it's probably better to pay the memory penalty.
# In the long run, it might be good to figure out what exactly is going on here but for now, this will do.
@_LAZY_REGISTRY.register_read(h5py.Dataset, IOSpec("string-array", "0.2.0"))
def read_h5_string_array(
    elem: h5py.Dataset,
    *,
    _reader: LazyReader,
    chunks: tuple[int | None, ...] | None = None,
) -> DaskArray:
    import dask.array as da

    return da.from_array(read_elem(elem), chunks=resolve_chunks(elem, chunks))


@_LAZY_REGISTRY.register_read(h5py.Dataset, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem: h5py.Dataset,
    *,
    _reader: LazyReader,
    chunks: tuple[int | None, ...] | None = None,
) -> DaskArray:
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name: str = elem.name
    shape = tuple(elem.shape)
    dtype = elem.dtype
    resolved_chunks = resolve_chunks(elem, chunks)

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_size(resolved_chunks[i], shape[i])
        for i in range(len(shape))
    )

    make_chunk = partial(make_dask_chunk, path, elem_name)
    return da.map_blocks(
        make_chunk, dtype=dtype, chunks=chunk_layout, meta=np.array([])
    )


@_LAZY_REGISTRY.register_read(zarr.Array, IOSpec("string-array", "0.2.0"))
@_LAZY_REGISTRY.register_read(zarr.Array, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem: zarr.Array,
    *,
    _reader: LazyReader,
    chunks: tuple[int | None, ...] | None = None,
) -> DaskArray:
    import dask.array as da

    return da.from_zarr(elem, chunks=chunks)


def _gen_xarray_dict_iterator_from_elems(
    elem_dict: Mapping[str, LazyDataStructures | pd.Index],
    dim_name: str,
    index: pd.Index,
) -> Generator[tuple[str, XVariable], None, None]:
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ...compat import xarray as xr

    for k, v in elem_dict.items():
        if isinstance(v, DaskArray) and k != dim_name:
            variable = xr.Variable([dim_name], data=v)
        elif isinstance(v, CategoricalArray | MaskedArray) and k != dim_name:
            variable = xr.Variable(
                [dim_name],
                data=xr.core.indexing.LazilyIndexedArray(v),
                attrs={
                    "base_path_or_zarr_group": v.base_path_or_zarr_group,
                    "elem_name": v.elem_name,
                    "is_nullable_string": (
                        isinstance(v, MaskedArray)
                        and isinstance(v.dtype, pd.StringDtype | np.dtypes.StringDType)
                    ),
                },
            )
        elif k == dim_name:
            variable = xr.Variable([dim_name], data=index)
        else:
            msg = f"Could not read {k}: {v} from into xarray Dataset2D"
            raise ValueError(msg)
        yield k, variable


DUMMY_RANGE_INDEX_KEY = "_anndata_dummy_range_index"


def _read_index(elem: StorageType) -> pd.Index:
    values = read_elem(elem)
    assert isinstance(values, np.ndarray | pd.api.extensions.ExtensionArray)
    return pd.Index(values)


@_LAZY_REGISTRY.register_read(zarr.Group, IOSpec("dataframe", "0.2.0"))
@_LAZY_REGISTRY.register_read(h5py.Group, IOSpec("dataframe", "0.2.0"))
@requires_xarray
def read_dataframe(
    elem: h5py.Group | zarr.Group,
    *,
    _reader: LazyReader,
    use_range_index: bool = False,
    chunks: tuple[int | None, ...] | None = None,
) -> Dataset2D:
    attrs: Mapping[str, Any] = elem.attrs
    # going through dask for reading into memory the index doesn't make sense, hence the ternary.
    elem_dict = {
        k: _reader.read_elem(elem[k], chunks=chunks)
        if (use_range_index and k == attrs["_index"]) or k != attrs["_index"]
        else _read_index(elem[k])
        for k in [*attrs["column-order"], attrs["_index"]]
    }
    if pd.api.types.is_string_dtype(elem_dict[attrs["_index"]]) and not use_range_index:
        elem_dict[attrs["_index"]] = pandas_as_str(elem_dict[attrs["_index"]])
    # If we use a range index, the coord axis needs to have the special dim name
    # which is used below as well.
    if not use_range_index:
        dim_name = attrs["_index"]
        # no sense in reading this in multiple times since xarray requires an in-memory index
        index = elem_dict[dim_name]
    else:
        dim_name = DUMMY_RANGE_INDEX_KEY
        index = pd.RangeIndex(len(elem_dict[attrs["_index"]])).astype("str")
    elem_xarray_dict = dict(
        _gen_xarray_dict_iterator_from_elems(elem_dict, dim_name, index)
    )
    if use_range_index:
        elem_xarray_dict[DUMMY_RANGE_INDEX_KEY] = XVariable(
            [DUMMY_RANGE_INDEX_KEY],
            data=index,
        )
    ds = Dataset2D(XDataset(elem_xarray_dict))
    ds.is_backed = True
    # We ensure the indexing_key attr always points to the true index
    # so that the roundtrip works even for the `use_range_index` `True` case
    ds.true_index_dim = attrs["_index"]
    return ds


@_LAZY_REGISTRY.register_read(zarr.Group, IOSpec("categorical", "0.2.0"))
@_LAZY_REGISTRY.register_read(h5py.Group, IOSpec("categorical", "0.2.0"))
@requires_xarray
def read_categorical(
    elem: h5py.Group | zarr.Group,
    *,
    _reader: LazyReader,
    chunks: tuple[int | None, ...] | None = None,
) -> CategoricalArray:
    from anndata.experimental.backed._lazy_arrays import CategoricalArray

    del chunks  # ignored when reading groups

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, h5py.Group) else elem
    )
    elem_name = get_elem_name(elem)
    attrs: Mapping[str, Any] = elem.attrs
    return CategoricalArray(
        codes=elem["codes"],
        categories=elem["categories"],
        ordered=attrs["ordered"],
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


@requires_xarray
def read_nullable(
    elem: h5py.Group | zarr.Group,
    *,
    encoding_type: Literal[
        "nullable-integer", "nullable-boolean", "nullable-string-array"
    ],
    _reader: LazyReader,
    chunks: tuple[int | None, ...] | None = None,
) -> MaskedArray:
    from anndata.experimental.backed._lazy_arrays import MaskedArray

    del chunks  # ignored when reading groups

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, h5py.Group) else elem
    )
    elem_name = get_elem_name(elem)
    if encoding_type == "nullable-string-array" and isinstance(elem, h5py.Group):
        # HDF5 stores strings as bytes; use .astype("T") to decode on access
        # h5py recommends .astype("T") over .asstr() when using numpy ≥2
        values = elem["values"].astype("T")
    else:
        values = elem["values"]
    return MaskedArray(
        values=values,
        mask=elem["mask"],
        dtype_str=encoding_type,
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


_NULLABLE_ENCODING_TYPES: Final = (
    "nullable-integer",
    "nullable-boolean",
    "nullable-string-array",
)

for encoding_type in _NULLABLE_ENCODING_TYPES:
    for group_type in (zarr.Group, h5py.Group):
        _LAZY_REGISTRY.register_read(group_type, IOSpec(encoding_type, "0.1.0"))(
            partial(read_nullable, encoding_type=encoding_type)
        )
