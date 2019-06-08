import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Union, MutableMapping

import pandas as pd
import math
import numpy as np
from scipy.sparse import issparse

from .. import AnnData
from .. import h5py
from ..compat import PathLike, fspath
from ..logging import get_logger


logger = get_logger(__name__)


def write_csvs(dirname: PathLike, adata: AnnData, skip_data: bool = True, sep: str = ','):
    """See :meth:`~anndata.AnnData.write_csvs`.
    """
    dirname = Path(dirname)
    if dirname.suffix == '.csv':
        dirname = dirname.with_suffix('')
    logger.info("writing '.csv' files to %s", dirname)
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=True)
    dir_uns = dirname / 'uns'
    if not dir_uns.is_dir():
        dir_uns.mkdir(parents=True, exist_ok=True)
    d = dict(
        obs=adata._obs,
        var=adata._var,
        obsm=adata._obsm.to_df(),
        varm=adata._varm.to_df(),
    )
    if not skip_data:
        d['X'] = pd.DataFrame(
            adata._X.toarray() if issparse(adata._X) else adata._X)
    d_write = {**d, **adata._uns}
    not_yet_raised_sparse_warning = True
    for key, value in d_write.items():
        if issparse(value):
            if not_yet_raised_sparse_warning:
                warnings.warn('Omitting to write sparse annotation.')
                not_yet_raised_sparse_warning = False
            continue
        filename = dirname
        if key not in {'X', 'var', 'obs', 'obsm', 'varm'}:
            filename = dir_uns
        filename /= '{}.csv'.format(key)
        df = value
        if not isinstance(value, pd.DataFrame):
            value = np.array(value)
            if np.ndim(value) == 0:
                value = value[None]
            try:
                df = pd.DataFrame(value)
            except Exception as e:
                warnings.warn('Omitting to write {!r}.'.format(key), type(e))
                continue
        df.to_csv(
            filename, sep=sep,
            header=key in {'obs', 'var', 'obsm', 'varm'},
            index=key in {'obs', 'var'},
        )


def write_loom(filename: PathLike, adata: AnnData, write_obsm_varm: bool = False):
    filename = Path(filename)
    row_attrs = {k: np.array(v) for k, v in adata.var.to_dict('list').items()}
    row_attrs['var_names'] = adata.var_names.values
    col_attrs = {k: np.array(v) for k, v in adata.obs.to_dict('list').items()}
    col_attrs['obs_names'] = adata.obs_names.values

    if adata.X is None:
        raise ValueError('loompy does not accept empty matrices as data')

    if write_obsm_varm:
        for key in adata.obsm.keys():
            col_attrs[key] = adata.obsm[key]
        for key in adata.varm.keys():
            row_attrs[key] = adata.varm[key]
    else:
        if len(adata.obsm.keys()) > 0 or len(adata.varm.keys()) > 0:
            logger.warning(
                'The loom file will lack these fields:\n{}\n'
                'Use write_obsm_varm=True to export multi-dimensional annotations'
                .format(adata.obsm.keys() + adata.varm.keys()))

    layers = {'': adata.X.T}
    for key in adata.layers.keys():
        layers[key] = adata.layers[key].T

    from loompy import create
    if filename.exists():
        filename.unlink()
    create(fspath(filename), layers, row_attrs=row_attrs, col_attrs=col_attrs)


def write_zarr(store: Union[MutableMapping, PathLike], adata: AnnData, **kwargs):
    if isinstance(store, Path):
        store = str(store)
    d = adata._to_dict_fixed_width_arrays(var_len_str=False)
    import zarr
    f = zarr.open(store, mode='w')
    for key, value in d.items():
        _write_key_value_to_zarr(f, key, value, **kwargs)


def _write_key_value_to_zarr(f, key, value, **kwargs):
    if isinstance(value, Mapping):
        for k, v in value.items():
            if not isinstance(k, str):
                warnings.warn('dict key {} transformed to str upon writing to zarr,'
                              'using string keys is recommended'
                              .format(k))
            _write_key_value_to_zarr(f, key + '/' + str(k), v, **kwargs)
        return

    def preprocess_writing(value):
        if value is None:
            return value
        elif issparse(value):
            return value
        elif isinstance(value, dict):
            # old hack for storing dicts, is never reached
            # in the current implementation, can be removed in the future
            value = np.array([str(value)])
        else:
            # make sure value is an array
            value = np.array(value)
            # hm, why that?
            if value.ndim == 0: value = np.array([value])
        # make sure string format is chosen correctly
        if value.dtype.kind == 'U': value = value.astype(np.string_)
        return value

    value = preprocess_writing(value)

    # for some reason, we need the following for writing string arrays
    if key in f.keys() and value is not None: del f[key]

    # ignore arrays with empty dtypes
    if value is None or not value.dtype.descr:
        return
    try:
        if key in set(f.keys()):
            import zarr
            is_valid_group = isinstance(f[key], zarr.hierarchy.Group) \
                             and f[key].shape == value.shape \
                             and f[key].dtype == value.dtype
            if not is_valid_group and not issparse(value):
                f[key][()] = value
                return
            else:
                del f[key]
        #f.create_dataset(key, data=value, **kwargs)
        if key != 'X' and 'chunks' in kwargs:  # TODO: make this more explicit
            del kwargs['chunks']
        import numcodecs  # TODO: only set object_codec for objects
        ds = f.create_dataset(key, shape=value.shape,
                                 dtype=value.dtype, object_codec=numcodecs.JSON(), **kwargs)
        _write_in_zarr_chunks(ds, key, value)
    except TypeError:
        # try writing as byte strings
        try:
            if value.dtype.names is None:
                if key in set(f.keys()):
                    if (f[key].shape == value.shape
                            and f[key].dtype == value.dtype):
                        f[key][()] = value.astype('S')
                        return
                    else:
                        del f[key]
                #f.create_dataset(key, data=value.astype('S'), **kwargs)
                ds = f.create_dataset(key, shape=value.astype('S').shape,
                                         dtype=value.astype('S').dtype, **kwargs)
                _write_in_zarr_chunks(ds, key, value.astype('S'))
            else:
                new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                             for dt in value.dtype.descr]
                if key in set(f.keys()):
                    if (f[key].shape == value.shape
                            and f[key].dtype == value.dtype):
                        f[key][()] = value.astype(new_dtype)
                        return
                    else:
                        del f[key]
                #f.create_dataset(
                #    key, data=value.astype(new_dtype), **kwargs)
                ds = f.create_dataset(key, shape=value.astype(new_dtype).shape,
                                      dtype=value.astype(new_dtype).dtype, **kwargs)
                _write_in_zarr_chunks(ds, key, value.astype(new_dtype))
        except Exception as e:
            warnings.warn('Could not save field with key = "{}" '
                          'to hdf5 file: {}'.format(key, e))


def _get_chunk_indices(za):
    # TODO: does zarr provide code for this?
    """
    Return all the indices (coordinates) for the chunks in a zarr array, even empty ones.
    """
    return [(i, j) for i in range(int(math.ceil(float(za.shape[0])/za.chunks[0])))
            for j in range(int(math.ceil(float(za.shape[1])/za.chunks[1])))]


def _write_in_zarr_chunks(za, key, value):
    if key != 'X':
        za[:] = value # don't chunk metadata
    else:
        for ci in _get_chunk_indices(za):
            s0, e0 = za.chunks[0] * ci[0], za.chunks[0] * (ci[0] + 1)
            s1, e1 = za.chunks[1] * ci[1], za.chunks[1] * (ci[1] + 1)
            print(ci, s0, e1, s1, e1)
            if issparse(value):
                za[s0:e0, s1:e1] = value[s0:e0, s1:e1].todense()
            else:
                za[s0:e0, s1:e1] = value[s0:e0, s1:e1]


def _write_h5ad(filename: PathLike, adata: AnnData, force_dense: bool = False, **kwargs):
    filename = Path(filename)
    if filename.suffix not in ('.h5', '.h5ad'):
        raise ValueError("Filename needs to end with '.h5ad'.")
    if adata.isbacked:
        # close so that we can reopen below
        adata.file.close()
    # create directory if it doesn't exist
    dirname = filename.parent
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=True)
    d = adata._to_dict_fixed_width_arrays()
    # we're writing to a different location than the backing file
    # - load the matrix into the memory...
    if adata.isbacked and filename != adata.filename:
        d['X'] = adata.X[:]
    # need to use 'a' if backed, otherwise we loose the backed objects
    with h5py.File(filename, 'a' if adata.isbacked else 'w', force_dense=force_dense) as f:
        for key, value in d.items():
            _write_key_value_to_h5(f, key, value, **kwargs)
    if adata.isbacked:
        adata.file.open(filename, 'r+')


def _write_key_value_to_h5(f, key, value, **kwargs):
    if isinstance(value, Mapping):
        for k, v in value.items():
            if not isinstance(k, str):
                warnings.warn(
                    'dict key {} transformed to str upon writing to h5,'
                    'using string keys is recommended'
                    .format(k)
                )
            _write_key_value_to_h5(f, key + '/' + str(k), v, **kwargs)
        return

    def preprocess_writing(value):
        if value is None or issparse(value):
            return value
        else:
            value = np.array(value)  # make sure value is an array
            if value.ndim == 0: value = np.array([value])  # hm, why that?
        # make sure string format is chosen correctly
        if value.dtype.kind == 'U': value = value.astype(h5py.special_dtype(vlen=str))
        return value

    value = preprocess_writing(value)

    # FIXME: for some reason, we need the following for writing string arrays
    if key in f.keys() and value is not None: del f[key]

    # ignore arrays with empty dtypes
    if value is None or not value.dtype.descr:
        return
    try:
        if key in set(f.keys()):
            is_valid_group = (
                isinstance(f[key], h5py.Group)
                and f[key].shape == value.shape
                and f[key].dtype == value.dtype
                and not isinstance(f[key], h5py.SparseDataset)
            )
            if not is_valid_group and not issparse(value):
                f[key][()] = value
                return
            else:
                del f[key]
        f.create_dataset(key, data=value, **kwargs)
    except TypeError:
        try:
            if value.dtype.names is None:
                dt = h5py.special_dtype(vlen=str)
                f.create_dataset(key, data=value, dtype=dt, **kwargs)
            else:  # try writing composite datatypes with byte strings
                new_dtype = [
                    (dt[0], 'S{}'.format(int(dt[1][2:])*4))
                    for dt in value.dtype.descr
                ]
                if key in set(f.keys()):
                    if (
                        f[key].shape == value.shape
                        and f[key].dtype == value.dtype
                    ):
                        f[key][()] = value.astype(new_dtype)
                        return
                    else:
                        del f[key]
                f.create_dataset(
                    key, data=value.astype(new_dtype), **kwargs)
        except Exception as e:
            warnings.warn(
                'Could not save field with key = {!r} to hdf5 file: {}'
                .format(key, e)
            )
