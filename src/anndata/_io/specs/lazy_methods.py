from __future__ import annotations

import re
from contextlib import contextmanager
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, overload

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import anndata as ad
from anndata._core.file_backing import filename, get_elem_name
from anndata.abc import CSCDataset, CSRDataset
from anndata.compat import DaskArray, H5Array, H5Group, ZarrArray, ZarrGroup

from .registry import _LAZY_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping, Sequence
    from typing import Literal, ParamSpec, TypeVar

    from anndata.experimental.backed._compat import DataArray, Dataset2D
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ...compat import H5File, SpArray
    from .registry import DaskReader, LazyDataStructures, LazyReader

    BlockInfo = Mapping[
        Literal[None],
        dict[str, Sequence[tuple[int, int]]],
    ]

    P = ParamSpec("P")
    R = TypeVar("R")
    D = TypeVar("D")


@overload
@contextmanager
def maybe_open_h5(
    path_or_other: Path, elem_name: str
) -> Generator[H5File, None, None]: ...
@overload
@contextmanager
def maybe_open_h5(path_or_other: D, elem_name: str) -> Generator[D, None, None]: ...
@contextmanager
def maybe_open_h5(
    path_or_other: H5File | D, elem_name: str
) -> Generator[H5File | D, None, None]:
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
    n_strides, rest = np.divmod(full_axis_size, chunk_axis_size)
    chunk = (chunk_axis_size,) * n_strides
    if rest > 0:
        chunk += (rest,)
    return chunk


def make_dask_chunk(
    path_or_sparse_dataset: Path | D,
    elem_name: str,
    block_info: BlockInfo | None = None,
) -> sparse.csr_matrix | sparse.csc_matrix | SpArray:
    if block_info is None:
        msg = "Block info is required"
        raise ValueError(msg)
    # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
    # https://github.com/scverse/anndata/issues/1105
    with maybe_open_h5(path_or_sparse_dataset, elem_name) as f:
        mtx = ad.io.sparse_dataset(f) if isinstance(f, H5Group) else f
        idx = tuple(
            slice(start, stop) for start, stop in block_info[None]["array-location"]
        )
        chunk = mtx[idx]
    return chunk


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: H5Group | ZarrGroup,
    *,
    _reader: DaskReader,
    chunks: tuple[int, ...] | None = None,  # only tuple[int, int] is supported here
) -> DaskArray:
    import dask.array as da

    path_or_sparse_dataset = (
        Path(filename(elem))
        if isinstance(elem, H5Group)
        else ad.io.sparse_dataset(elem)
    )
    elem_name = get_elem_name(elem)
    shape: tuple[int, int] = tuple(elem.attrs["shape"])
    if isinstance(path_or_sparse_dataset, CSRDataset | CSCDataset):
        dtype = path_or_sparse_dataset.dtype
    else:
        dtype = elem["data"].dtype
    is_csc: bool = elem.attrs["encoding-type"] == "csc_matrix"

    stride: int = _DEFAULT_STRIDE
    major_dim, minor_dim = (1, 0) if is_csc else (0, 1)
    if chunks is not None:
        if len(chunks) != 2:
            raise ValueError("`chunks` must be a tuple of two integers")
        if chunks[minor_dim] not in {shape[minor_dim], -1, None}:
            raise ValueError(
                "Only the major axis can be chunked. "
                f"Try setting chunks to {((-1, _DEFAULT_STRIDE) if is_csc else (_DEFAULT_STRIDE, -1))}"
            )
        stride = (
            chunks[major_dim]
            if chunks[major_dim] not in {None, -1}
            else shape[major_dim]
        )

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


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("string-array", "0.2.0"))
def read_h5_string_array(
    elem: H5Array,
    *,
    _reader: LazyReader,
    chunks: tuple[int, int] | None = None,
) -> DaskArray:
    import dask.array as da

    from anndata._io.h5ad import read_dataset

    return da.from_array(
        read_dataset(elem),
        chunks=chunks if chunks is not None else (_DEFAULT_STRIDE,) * len(elem.shape),
    )


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem: H5Array, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name: str = elem.name
    shape = tuple(elem.shape)
    dtype = elem.dtype
    chunks: tuple[int, ...] = (
        tuple(
            c if c not in {None, -1} else s for c, s in zip(chunks, shape, strict=True)
        )
        if chunks is not None
        else tuple(min(_DEFAULT_STRIDE, s) for s in shape)
    )

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_size(chunks[i], shape[i])
        for i in range(len(shape))
    )

    make_chunk = partial(make_dask_chunk, path, elem_name)
    return da.map_blocks(
        make_chunk, dtype=dtype, chunks=chunk_layout, meta=np.array([])
    )


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("string-array", "0.2.0"))
@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem: ZarrArray, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    chunks: tuple[int, ...] = chunks if chunks is not None else elem.chunks
    import dask.array as da

    return da.from_zarr(elem, chunks=chunks)


DUMMY_RANGE_INDEX_KEY = "_anndata_dummy_range_index"


def _gen_xarray_dict_iterator_from_elems(
    elem_dict: dict[str, LazyDataStructures],
    index_label: str,
    index_key: str,
    index: np.NDArray,
) -> Generator[tuple[str, DataArray], None, None]:
    from anndata.experimental.backed._compat import DataArray
    from anndata.experimental.backed._compat import xarray as xr
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    for k, v in elem_dict.items():
        data_array_name = k
        if isinstance(v, DaskArray) and k != index_key:
            data_array = DataArray(v, coords=[index], dims=[index_label], name=k)
        elif isinstance(v, CategoricalArray | MaskedArray) and k != index_key:
            variable = xr.Variable(
                data=xr.core.indexing.LazilyIndexedArray(v), dims=[index_label]
            )
            data_array = DataArray(
                variable,
                coords=[index],
                dims=[index_label],
                name=k,
                attrs={
                    "base_path_or_zarr_group": v.base_path_or_zarr_group,
                    "elem_name": v.elem_name,
                },
            )
        elif k == index_key:
            data_array = DataArray(
                index, coords=[index], dims=[index_label], name=index_label
            )
            data_array_name = index_label
        else:
            raise ValueError(f"Could not read {k}: {v} from into xarray Dataset2D")
        yield data_array_name, data_array
    if index_key == DUMMY_RANGE_INDEX_KEY:
        yield (
            index_label,
            DataArray(index, coords=[index], dims=[index_label], name=index_label),
        )


@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.2.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.2.0"))
def read_dataframe(
    elem: H5Group | ZarrGroup,
    *,
    _reader: LazyReader,
    use_range_index: bool = False,
) -> Dataset2D:
    from anndata.experimental.backed._compat import Dataset2D

    elem_dict = {
        k: _reader.read_elem(elem[k])
        for k in [*elem.attrs["column-order"], elem.attrs["_index"]]
    }
    elem_name = get_elem_name(elem)
    # Determine whether we can use label based indexing i.e., is the elem `obs` or `var`
    obs_var_matches = re.findall(r"(obs|var)", elem_name)
    if not len(obs_var_matches) == 1:
        label_based_indexing_key = "index"
    else:
        label_based_indexing_key = f"{obs_var_matches[0]}_names"
    # If we are not using a range index, the underlying on disk label for the index
    # could be different than {obs,var}_names - otherwise we use a dummy value.
    if not use_range_index:
        index_label = label_based_indexing_key
        index_key = elem.attrs["_index"]
        # no sense in reading this in multiple times
        index = elem_dict[index_key].compute()
    else:
        index_label = DUMMY_RANGE_INDEX_KEY
        index_key = DUMMY_RANGE_INDEX_KEY
        index = pd.RangeIndex(len(elem_dict[elem.attrs["_index"]]))
    elem_xarray_dict = dict(
        _gen_xarray_dict_iterator_from_elems(elem_dict, index_label, index_key, index)
    )
    ds = Dataset2D(elem_xarray_dict, attrs={"indexing_key": label_based_indexing_key})
    if use_range_index:
        return ds.rename_vars({elem.attrs["_index"]: label_based_indexing_key})
    return ds


@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("categorical", "0.2.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("categorical", "0.2.0"))
def read_categorical(
    elem: H5Group | ZarrGroup,
    *,
    _reader: LazyReader,
) -> CategoricalArray:
    from anndata.experimental.backed._lazy_arrays import CategoricalArray

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, H5Group) else elem
    )
    elem_name = get_elem_name(elem)
    return CategoricalArray(
        codes=elem["codes"],
        categories=elem["categories"],
        ordered=elem.attrs["ordered"],
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


def read_nullable(
    elem: H5Group | ZarrGroup,
    *,
    encoding_type: Literal["nullable-integer", "nullable-boolean"],
    _reader: LazyReader,
) -> MaskedArray:
    from anndata.experimental.backed._lazy_arrays import MaskedArray

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, H5Group) else elem
    )
    elem_name = get_elem_name(elem)
    return MaskedArray(
        values=elem["values"],
        mask=elem["mask"] if "mask" in elem else None,
        dtype_str=encoding_type,
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("nullable-integer", "0.1.0"))(
    partial(read_nullable, encoding_type="nullable-integer")
)
_LAZY_REGISTRY.register_read(H5Group, IOSpec("nullable-integer", "0.1.0"))(
    partial(read_nullable, encoding_type="nullable-integer")
)
_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("nullable-boolean", "0.1.0"))(
    partial(read_nullable, encoding_type="nullable-boolean")
)
_LAZY_REGISTRY.register_read(H5Group, IOSpec("nullable-boolean", "0.1.0"))(
    partial(read_nullable, encoding_type="nullable-boolean")
)
