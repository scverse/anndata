"""Main class and helper functions.
"""
from enum import Enum
from collections import OrderedDict
from collections.abc import MutableMapping
from functools import reduce
from pathlib import Path
from typing import Any, Union, Optional
from typing import Iterable, Sized, Sequence, Mapping
from typing import Tuple, List
from copy import deepcopy
import warnings

import numpy as np
from numpy import ma
import pandas as pd
from pandas.core.index import RangeIndex
from pandas.api.types import is_string_dtype, is_categorical
from scipy import sparse
from scipy.sparse import issparse
from natsort import natsorted


# try importing zarr
try:
    from zarr.core import Array as ZarrArray
except ImportError:
    class ZarrArray:
        @staticmethod
        def __rep__():
            return 'mock zarr.core.Array'

# try importing zappy
try:
    from zappy.base import ZappyArray
except ImportError:
    class ZappyArray:
        @staticmethod
        def __rep__():
            return 'mock zappy.base.ZappyArray'


from .alignedmapping import (
    AxisArraysBase, AxisArrays,
    # PairwiseArraysBase, PairwiseArrays,
    LayersBase, Layers
)
from .. import h5py
from .views import ArrayView, SparseCSRView, SparseCSCView, DictView, DataFrameView

from .. import utils
from ..utils import Index, get_n_items_idx, convert_to_dict, unpack_index
from ..logging import anndata_logger as logger
from ..compat import PathLike, warn_flatten


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sparse.spmatrix
    ZarrArry = ZarrArray
    ZappyArry = ZappyArray

    @classmethod
    def classes(cls):
        print(ZarrArray)
        return tuple(c.value for c in cls.__members__.values())


# for backwards compat
def _find_corresponding_multicol_key(key, keys_multicol):
    """Find the corresponding multicolumn key."""
    for mk in keys_multicol:
        if key.startswith(mk) and 'of' in key:
            return mk
    return None


# for backwards compat
def _gen_keys_from_multicol_key(key_multicol, n_keys):
    """Generates single-column keys from multicolumn key."""
    keys = [('{}{:03}of{:03}')
            .format(key_multicol, i+1, n_keys) for i in range(n_keys)]
    return keys


def df_to_records_fixed_width(df, var_len_str=True):
    uns = {}  # unstructured dictionary for storing categories
    names = ['index']
    if is_string_dtype(df.index):
        if var_len_str:
            index = df.index.values.astype(h5py.special_dtype(vlen=str))
        else:
            max_len_index = 0 if 0 in df.shape else df.index.map(len).max()
            index = df.index.values.astype('S{}'.format(max_len_index))
    else:
        index = df.index.values
    arrays = [index]
    for k in df.columns:
        names.append(k)
        if is_string_dtype(df[k]) and not is_categorical(df[k]):
            if var_len_str:
                arrays.append(df[k].values.astype(h5py.special_dtype(vlen=str)))
            else:
                lengths = df[k].map(len)
                if is_categorical(lengths): lengths = lengths.cat.as_ordered()
                arrays.append(df[k].values.astype('S{}'.format(lengths.max())))
        elif is_categorical(df[k]):
            uns[k + '_categories'] = df[k].cat.categories
            arrays.append(df[k].cat.codes)
        else:
            arrays.append(df[k].values)
    formats = [v.dtype for v in arrays]
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}), uns


def _check_2d_shape(X):
    """Check shape of array or sparse matrix.

    Assure that X is always 2D: Unlike numpy we always deal with 2D arrays.
    """
    if X.dtype.names is None and len(X.shape) != 2:
        raise ValueError('X needs to be 2-dimensional, not '
                         '{}-dimensional.'.format(len(X.shape)))


def _normalize_index(index, names):
    if not isinstance(names, RangeIndex):
        assert names.dtype != float and names.dtype != int, \
            'Don’t call _normalize_index with non-categorical/string names'

    # the following is insanely slow for sequences, we replaced it using pandas below
    def name_idx(i):
        if isinstance(i, str):
            # `where` returns an 1-tuple (1D array) of found indices
            i_found = np.where(names == i)[0]
            if len(i_found) == 0:  # returns array of length 0 if nothing is found
                raise IndexError(
                    'Key "{}" is not valid observation/variable name/index.'
                    .format(i))
            i = i_found[0]
        return i

    if isinstance(index, slice):
        start = name_idx(index.start)
        stop = name_idx(index.stop)
        # string slices can only be inclusive, so +1 in that case
        if isinstance(index.stop, str):
            stop = None if stop is None else stop + 1
        step = index.step
        return slice(start, stop, step)
    elif isinstance(index, (np.integer, int, str)):
        return name_idx(index)
    elif isinstance(index, (Sequence, np.ndarray, pd.Index)):
        # here, we replaced the implementation based on name_idx with this
        # incredibly faster one
        positions = pd.Series(index=names, data=range(len(names)))
        positions = positions[index]
        if positions.isnull().values.any():
            not_found = positions.index[positions.isnull().values]
            raise KeyError(
                "Values {}, from {}, are not valid obs/ var names or indices."
                .format(list(not_found), list(positions.index)))
        return positions.values
    else:
        raise IndexError('Unknown index {!r} of type {}'
                         .format(index, type(index)))


def _gen_dataframe(anno, length, index_names):
    if isinstance(anno, pd.DataFrame):
        anno = anno.copy()
        if not is_string_dtype(anno.index):
            logger.warning('Transforming to str index.')
            anno.index = anno.index.astype(str)
        return anno
    if anno is None or len(anno) == 0:
        _anno = pd.DataFrame(index=RangeIndex(0, length, name=None).astype(str))
    else:
        for index_name in index_names:
            if index_name in anno:
                _anno = pd.DataFrame(
                    anno, index=anno[index_name],
                    columns=[k for k in anno.keys() if k != index_name])
                break
        else:
            _anno = pd.DataFrame(anno, index=RangeIndex(0, length, name=None).astype(str))
    return _anno


class AnnDataFileManager:
    """Backing file manager for AnnData.
    """

    def __init__(
        self,
        adata: 'AnnData',
        filename: Optional[PathLike] = None,
        filemode: Optional[str] = None,
    ):
        self._adata = adata
        self.filename = filename
        self._filemode = filemode
        self._file = None
        if filename:
            self.open()

    def __repr__(self) -> str:
        if self.filename is None:
            return 'Backing file manager: no file is set.'
        else:
            return 'Backing file manager of file {}.'.format(self.filename)

    def __getitem__(self, key: str) -> Union[h5py.Group, h5py.Dataset, h5py.SparseDataset]:
        return self._file[key]

    def __setitem__(self, key: str, value: Union[h5py.Group, h5py.Dataset, h5py.SparseDataset]):
        self._file[key] = value

    def __delitem__(self, key: str):
        del self._file[key]

    @property
    def filename(self) -> Path:
        return self._filename

    @filename.setter
    def filename(self, filename: Optional[PathLike]):
        self._filename = None if filename is None else Path(filename)

    def open(
        self,
        filename: Optional[PathLike] = None,
        filemode: Optional[str] = None,
    ):
        if filename is not None:
            self.filename = filename
        if filemode is not None:
            self._filemode = filemode
        if self.filename is None:
            raise ValueError(
                'Cannot open backing file if backing not initialized.')
        self._file = h5py.File(self.filename, self._filemode, force_dense=True)

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        if self._file is not None:
            self._file.close()

    def _to_memory_mode(self):
        """Close the backing file, forget filename, *do* change to memory mode."""
        self._adata.__X = self._adata.X[()]
        self._file.close()
        self._file = None
        self._filename = None

    @property
    def isopen(self) -> bool:
        """State of backing file."""
        if self._file is None:
            return False
        # try accessing the id attribute to see if the file is open
        return bool(self._file.id)


class Raw:
    def __init__(
        self,
        adata: Optional['AnnData'] = None,
        X: Union[np.ndarray, sparse.spmatrix, None] = None,
        var: Optional[AxisArraysBase] = None,
        varm: Optional[AxisArraysBase] = None,
    ):
        self._adata = adata
        self._n_obs = adata.n_obs
        if X is not None:
            self._X = X
            self._var = var
            self._varm = varm
        else:
            self._X = None if adata.isbacked else adata.X.copy()
            self._var = adata.var.copy()
            self._varm = adata.varm.copy()

    @property
    def X(self):
        if self._adata.isbacked:
            if not self._adata.file.isopen: self._adata.file.open()
            X = self._adata.file['raw.X']
            if self._adata.isview: return X[self._adata._oidx, self._adata._vidx]
            else: return X
        else:
            if self.n_obs == 1 and self.n_vars == 1:
                warn_flatten()
                return self._X[0, 0]
            elif self.n_obs == 1 or self.n_vars == 1:
                warn_flatten()
                X = self._X
                if issparse(self._X): X = self._X.toarray()
                return X.flatten()
            else:
                return self._X

    @property
    def shape(self):
        return self.X.shape

    @property
    def var(self):
        return self._var

    @property
    def n_vars(self):
        return self._var.shape[0]

    @property
    def n_obs(self):
        return self._n_obs

    @property
    def varm(self):
        return self._varm

    @property
    def var_names(self):
        return self.var.index

    @property
    def obs_names(self):
        return self._adata.obs_names

    def __getitem__(self, index):
        oidx, vidx = self._normalize_indices(index)
        if self._adata is not None or not self._adata.isbacked: X = self._X[oidx, vidx]
        else: X = self._adata.file['raw.X'][oidx, vidx]
        if isinstance(vidx, (int, np.integer)):
            # To preserve two dimensional shape
            vidx = slice(vidx, vidx + 1, 1)
        var = self._var.iloc[vidx]
        new = Raw(self._adata, X=X, var=var)
        if self._varm is not None:
            new._varm = self._varm._view(self, vidx)
        return new

    def copy(self):
        return Raw(self._adata, X=self._X.copy(), var=self._var.copy(),
                   varm=None if self._varm is None else self._varm.copy())

    def _normalize_indices(self, packed_index):
        # deal with slicing with pd.Series
        if isinstance(packed_index, pd.Series):
            packed_index = packed_index.values
        if isinstance(packed_index, tuple):
            if len(packed_index) != 2:
                raise IndexDimError(len(packed_index))
            if isinstance(packed_index[1], pd.Series):
                packed_index = packed_index[0], packed_index[1].values
            if isinstance(packed_index[0], pd.Series):
                packed_index = packed_index[0].values, packed_index[1]
        obs, var = unpack_index(packed_index)
        obs = _normalize_index(obs, self._adata.obs_names)
        var = _normalize_index(var, self.var_names)
        return obs, var

    def var_vector(self, k: str) -> np.ndarray:
        """
        Convenience function for returning a 1 dimensional ndarray of values
        from `.X` or `.var`.

        Made for convenience, not performance. Intentionally permissive about
        arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in `.obs_names` or `.var.columns`.

        Returns
        -------
        A one dimensional nd array, with values for each var in the same order
        as `.var_names`.
        """
        if k in self.var:
            return self.var[k].values
        else:
            a = self[k, :].X
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)

    def obs_vector(self, k: str) -> np.ndarray:
        """
        Convenience function for returning a 1 dimensional ndarray of values
        from `.X`.

        Made for convenience, not performance. Intentionally permissive about
        arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in `.var_names` or `.obs.columns`. If `use_raw`,
            value should be in `.raw.var_names` instead of `.var_names`.

        Returns
        -------
        A one dimensional nd array, with values for each obs in the same order
        as `.obs_names`.
        """
        a = self[:, k].X
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)


INDEX_DIM_ERROR_MSG = 'You tried to slice an AnnData(View) object with an' \
    '{}-dimensional index, but only 2 dimensions exist in such an object.'
INDEX_DIM_ERROR_MSG_1D = '\nIf you tried to slice cells using adata[cells, ], ' \
    'be aware that Python (unlike R) uses adata[cells, :] as slicing syntax.'


class IndexDimError(IndexError):
    def __init__(self, n_dims):
        msg = INDEX_DIM_ERROR_MSG.format(n_dims)
        if n_dims == 1:
            msg += INDEX_DIM_ERROR_MSG_1D
        super().__init__(msg)


class AnnData(metaclass=utils.DeprecationMixinMeta):
    """An annotated data matrix.

    :class:`~anndata.AnnData` stores a data matrix :attr:`X` together with annotations
    of observations :attr:`obs`, variables :attr:`var` and unstructured annotations :attr:`uns`.

    .. figure:: https://falexwolf.de/img/scanpy/anndata.svg
       :width: 350px

    An :class:`~anndata.AnnData` object ``adata`` can be sliced like a pandas
    dataframe, for instance, ``adata_subset = adata[:, list_of_variable_names]``.
    :class:`~anndata.AnnData`'s basic structure is similar to R's ExpressionSet
    [Huber15]_. If setting an ``.h5ad``-formatted HDF5 backing file ``.filename``,
    data remains on the disk but is automatically loaded into memory if needed.
    See this `blog post`_ for more details.

    .. _blog post: http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/

    Parameters
    ----------
    X
        A #observations × #variables data matrix. A view of the data is used if the
        data type matches, otherwise, a copy is made.
    obs
        Key-indexed one-dimensional observations annotation of length #observations.
    var
        Key-indexed one-dimensional variables annotation of length #variables.
    uns
        Key-index unstructured annotation.
    obsm
        Key-indexed multi-dimensional observations annotation of length #observations.
        If passing a :class:`~numpy.ndarray`, it needs to have a structured datatype.
    varm
        Key-indexed multi-dimensional variables annotation of length #variables.
        If passing a :class:`~numpy.ndarray`, it needs to have a structured datatype.
    dtype
        Data type used for storage.
    shape
        Shape tuple (#observations, #variables). Can only be provided if ``X`` is ``None``.
    filename
        Name of backing file. See :class:`anndata.h5py.File`.
    filemode
        Open mode of backing file. See :class:`anndata.h5py.File`.
    layers
        Dictionary with keys as layers' names and values as matrices of the same dimensions as X.

    See Also
    --------
    read_h5ad
    read_csv
    read_excel
    read_hdf
    read_loom
    read_zarr
    read_mtx
    read_text
    read_umi_tools

    Notes
    -----
    Multi-dimensional annotations are stored in :attr:`obsm` and :attr:`varm`.

    Indexing into an AnnData object with a numeric is supposed to be positional,
    like pandas’ :attr:`~pandas.DataFrame.iloc` accessor, while indexing with a string/categorical is
    supposed to behave like :attr:`~pandas.DataFrame.loc`.

    If the unstructured annotations :attr:`uns` contain a sparse matrix of shape
    :attr:`n_obs` × :attr:`n_obs`, these are sliced when calling ``[]``.

    A data matrix is flattened if either :attr:`n_obs` or :attr:`n_vars` is 1, so that
    numpy's slicing behavior is reproduced::

        adata = AnnData(np.ones((2, 2)))
        adata[:, 0].X == adata.X[:, 0]

    :class:`~anndata.AnnData` stores observations (samples) of variables
    (features) in the rows of a matrix. This is the convention of the modern
    classics of statistics [Hastie09]_ and machine learning [Murphy12]_, the
    convention of dataframes both in R and Python and the established statistics
    and machine learning packages in Python (statsmodels_, scikit-learn_).

    .. _statsmodels: http://www.statsmodels.org/stable/index.html
    .. _scikit-learn: http://scikit-learn.org/
    """

    _BACKED_ATTRS = ['X', 'raw.X']

    # backwards compat
    _H5_ALIASES = {
        'X': {'X', '_X', 'data', '_data'},
        'obs': {'obs', '_obs', 'smp', '_smp'},
        'var': {'var', '_var'},
        'uns': {'uns'},
        'obsm': {'obsm', '_obsm', 'smpm', '_smpm'},
        'varm': {'varm', '_varm'},
        'layers': {'layers', '_layers'},
    }

    _H5_ALIASES_NAMES = {
        'obs': {'obs_names', 'smp_names', 'row_names', 'index'},
        'var': {'var_names', 'col_names', 'index'},
    }

    def __init__(
        self,
        X: Optional[Union[np.ndarray, sparse.spmatrix, pd.DataFrame]] = None,
        obs: Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]] = None,
        var: Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]] = None,
        uns: Optional[Mapping[str, Any]] = None,
        obsm: Optional[Union[np.ndarray, Mapping[str, Sequence[Any]]]] = None,
        varm: Optional[Union[np.ndarray, Mapping[str, Sequence[Any]]]] = None,
        layers: Optional[Mapping[str, Union[np.ndarray, sparse.spmatrix]]] = None,
        raw: Optional[Raw] = None,
        dtype: Union[np.dtype, str] = 'float32',
        shape: Optional[Tuple[int, int]] = None,
        filename: Optional[PathLike] = None,
        filemode: Optional[str] = None,
        asview: bool = False,
        *,
        # obsp: Optional[Union[np.ndarray, Mapping[str, Sequence[Any]]]] = None,
        # varp: Optional[Union[np.ndarray, Mapping[str, Sequence[Any]]]] = None,
        oidx: Index = None,
        vidx: Index = None
    ):
        if asview:
            if not isinstance(X, AnnData):
                raise ValueError('`X` has to be an AnnData object.')
            self._init_as_view(X, oidx, vidx)
        else:
            self._init_as_actual(
                X=X, obs=obs, var=var, uns=uns,
                obsm=obsm, varm=varm,
                # obsp=obsp, varp=varp,
                raw=raw,
                layers=layers,
                dtype=dtype, shape=shape,
                filename=filename, filemode=filemode)

    def _init_as_view(self, adata_ref: 'AnnData', oidx: Index, vidx: Index):
        if adata_ref.isbacked and adata_ref.isview:
            raise ValueError(
                'Currently, you cannot index repeatedly into a backed AnnData, '
                'that is, you cannot make a view of a view.')
        self._isview = True
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx
        # the file is the same as of the reference object
        self.file = adata_ref.file
        # views on attributes of adata_ref
        oidx_normalized, vidx_normalized = oidx, vidx
        if isinstance(oidx, (int, np.integer)): oidx_normalized = slice(oidx, oidx+1, 1)
        if isinstance(vidx, (int, np.integer)): vidx_normalized = slice(vidx, vidx+1, 1)
        obs_sub = adata_ref.obs.iloc[oidx_normalized]
        var_sub = adata_ref.var.iloc[vidx_normalized]
        self._obsm = adata_ref.obsm._view(self, oidx_normalized)
        self._varm = adata_ref.varm._view(self, vidx_normalized)
        # TODO: Should this be normalized
        self._layers = adata_ref.layers._view(self, (oidx, vidx))
        # self._obsp = adata_ref.obsp._view(self, oidx_normalized)
        # self._varp = adata_ref.varp._view(self, vidx_normalized)
        # hackish solution here, no copy should be necessary
        uns_new = deepcopy(self._adata_ref._uns)
        # need to do the slicing before setting the updated self._n_obs, self._n_vars
        self._n_obs = self._adata_ref.n_obs  # use the original n_obs here
        self._slice_uns_sparse_matrices_inplace(uns_new, self._oidx)
        # fix _n_obs, _n_vars
        if isinstance(oidx, slice):
            self._n_obs = get_n_items_idx(obs_sub.index, adata_ref.n_obs)
        elif isinstance(oidx, (int, np.integer)):
            self._n_obs = 1
        elif isinstance(oidx, Sized):
            self._n_obs = get_n_items_idx(oidx, adata_ref.n_obs)
        else:
            raise KeyError('Unknown Index type')
        if isinstance(vidx, slice):
            self._n_vars = get_n_items_idx(var_sub.index, adata_ref.n_vars)
        elif isinstance(vidx, (int, np.integer)):
            self._n_vars = 1
        elif isinstance(vidx, Sized):
            self._n_vars = get_n_items_idx(vidx, adata_ref.n_vars)
        else:
            raise KeyError('Unknown Index type')
        # fix categories
        self._remove_unused_categories(adata_ref.obs, obs_sub, uns_new)
        self._remove_unused_categories(adata_ref.var, var_sub, uns_new)
        # set attributes
        self._obs = DataFrameView(obs_sub, view_args=(self, 'obs'))
        self._var = DataFrameView(var_sub, view_args=(self, 'var'))
        self._uns = DictView(uns_new, view_args=(self, 'uns'))
        # set data
        if self.isbacked:
            self._X = None
        else:
            self._init_X_as_view()

        # set raw, easy, as it's immutable anyways...
        if adata_ref._raw is not None:
            # slicing along variables axis is ignored
            self._raw = adata_ref.raw[oidx]
        else:
            self._raw = None

    def _init_X_as_view(self):
        if self._adata_ref.X is None:
            self._X = None
            return
        X = self._adata_ref._X[self._oidx, self._vidx]
        if isinstance(X, sparse.csr_matrix):
            self._X = SparseCSRView(X, view_args=(self, 'X'))
        elif isinstance(X, sparse.csc_matrix):
            self._X = SparseCSCView(X, view_args=(self, 'X'))
        elif issparse(X):
            raise ValueError('View on non-csr/csc sparse matrices not implemented.')
        elif isinstance(X, ZappyArray): # ZappyArray acts as a view itself
            self._X = X
        else:
            shape = (
                get_n_items_idx(self._oidx, self._adata_ref.n_obs),
                get_n_items_idx(self._vidx, self._adata_ref.n_vars)
            )
            if np.isscalar(X):
                X = X.view()
            self._X = ArrayView(X.reshape(shape), view_args=(self, 'X'))

    def _init_as_actual(
            self, X=None, obs=None, var=None, uns=None,
            obsm=None, varm=None,
            # varp=None, obsp=None,
            raw=None, layers=None,
            dtype='float32', shape=None,
            filename=None, filemode=None):
        from ..readwrite.read import _read_args_from_h5ad

        # view attributes
        self._isview = False
        self._adata_ref = None
        self._oidx = None
        self._vidx = None

        # ----------------------------------------------------------------------
        # various ways of initializing the data
        # ----------------------------------------------------------------------

        # init from file
        if filename is not None:
            if any((X, obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If initializing from `filename`, '
                    'no further arguments may be passed.')
            self.file = AnnDataFileManager(self, filename, filemode)
            X, obs, var, uns, obsm, varm, layers, raw = _read_args_from_h5ad(self, mode=filemode)
            if X is not None:
                # this is not a function that a user would use, hence it's fine to set the dtype
                dtype = X.dtype.name
        else:
            self.file = AnnDataFileManager(self, None)

            # init from AnnData
            if isinstance(X, AnnData):
                if any((obs, var, uns, obsm, varm)):
                    raise ValueError(
                        'If `X` is a dict no further arguments must be provided.')
                X, obs, var, uns, obsm, varm, layers, raw = X._X, X.obs, X.var, X.uns, X.obsm, X.varm, X.layers, X.raw

            # init from DataFrame
            elif isinstance(X, pd.DataFrame):
                if obs is None:
                    obs = pd.DataFrame(index=X.index)
                else:
                    if not X.index.equals(obs.index):
                        raise ValueError(
                            'Index of obs must match index of X.')
                if var is None:
                    var = pd.DataFrame(index=X.columns)
                else:
                    if not X.columns.equals(var.index):
                        raise ValueError(
                            'Index of var must match columns of X.')
                X = X.values

        # ----------------------------------------------------------------------
        # actually process the data
        # ----------------------------------------------------------------------

        # check data type of X
        if X is not None:
            for s_type in StorageType:
                if isinstance(X, s_type.value):
                    break
            else:
                class_names = ', '.join(c.__name__ for c in StorageType.classes())
                raise ValueError('`X` needs to be of one of {}, not {}.'
                                 .format(class_names, type(X)))
            if shape is not None:
                raise ValueError('`shape` needs to be `None` is `X` is not `None`.')
            _check_2d_shape(X)
            # if type doesn't match, a copy is made, otherwise, use a view
            if issparse(X) or isinstance(X, ma.MaskedArray):
                # TODO: maybe use view on data attribute of sparse matrix
                #       as in readwrite.read_10x_h5
                if X.dtype != np.dtype(dtype): X = X.astype(dtype)
            elif isinstance(X, ZarrArray):
                X = X.astype(dtype)
            else:  # is np.ndarray
                X = X.astype(dtype, copy=False)
            # data matrix and shape
            self._X = X
            self._n_obs, self._n_vars = self._X.shape
        else:
            self._X = None
            self._n_obs = len([] if obs is None else obs)
            self._n_vars = len([] if var is None else var)
            # check consistency with shape
            if shape is not None:
                if self._n_obs == 0:
                    self._n_obs = shape[0]
                else:
                    if self._n_obs != shape[0]:
                        raise ValueError('`shape` is inconsistent with `obs`')
                if self._n_vars == 0:
                    self._n_vars = shape[1]
                else:
                    if self._n_vars != shape[1]:
                        raise ValueError('`shape` is inconsistent with `var`')

        # annotations
        self._obs = _gen_dataframe(obs, self._n_obs,
                                   ['obs_names', 'row_names', 'smp_names'])
        self._var = _gen_dataframe(var, self._n_vars, ['var_names', 'col_names'])

        # unstructured annotations
        self._uns = uns or OrderedDict()

        # TODO: Think about consequences of making obsm a group in hdf
        self._obsm = AxisArrays(self, 0, vals=convert_to_dict(obsm))
        self._varm = AxisArrays(self, 1, vals=convert_to_dict(varm))

        # self._obsp = PairwiseArrays(self, 0, vals=convert_to_dict(obsp))
        # self._varp = PairwiseArrays(self, 1, vals=convert_to_dict(varp))

        self._check_dimensions()
        self._check_uniqueness()

        # raw
        if raw is None:
            self._raw = None
        else:
            if isinstance(raw, Raw):
                self._raw = raw
            else:
                # is dictionary from reading the file, nothing that is meant for a user
                shape = self.file['raw.X'].shape if self.isbacked else raw['X'].shape

                self._raw = Raw(
                    self,
                    X=raw['X'],
                    var=_gen_dataframe(raw['var'], shape[1], ['var_names', 'col_names']),
                    varm=raw['varm'] if 'varm' in raw else None)

        # clean up old formats
        self._clean_up_old_format(uns)

        # layers
        self._layers = Layers(self, layers)

    def __sizeof__(self) -> int:
        size = 0
        for attr in ['_X', '_obs', '_var', '_uns', '_obsm', '_varm']:
            s = getattr(self, attr).__sizeof__()
            size += s
        return size

    def _gen_repr(self, n_obs, n_vars) -> str:
        if self.isbacked:
            backed_at = 'backed at \'{}\''.format(self.filename)
        else:
            backed_at = ''
        descr = (
            'AnnData object with n_obs × n_vars = {} × {} {}'
            .format(n_obs, n_vars, backed_at))
        for attr in ['obs', 'var', 'uns', 'obsm', 'varm', 'layers']:  # 'obsp', 'varp',
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += '\n    {}: {}'.format(attr, str(list(keys))[1:-1])
        return descr

    def __repr__(self) -> str:
        if self.isview:
            return 'View of ' + self._gen_repr(self.n_obs, self.n_vars)
        else:
            return self._gen_repr(self.n_obs, self.n_vars)

    @property
    def shape(self) -> Tuple[int, int]:
        r"""Shape of data matrix (:attr:`n_obs`, :attr:`n_vars`)."""
        return self.n_obs, self.n_vars

    @property
    def X(self) -> Optional[Union[np.ndarray, sparse.spmatrix, ArrayView]]:
        """Data matrix of shape :attr:`n_obs` × :attr:`n_vars`."""
        if self.isbacked:
            if not self.file.isopen: self.file.open()
            X = self.file['X']
            if self.isview:
                X = X[self._oidx, self._vidx]
            return X
        else:
            if self.n_obs == 1 and self.n_vars == 1:
                warn_flatten()
                return self._X[0, 0]
            elif self.n_obs == 1 or self.n_vars == 1:
                warn_flatten()
                X = self._X
                if issparse(self._X): X = self._X.toarray()
                return X.flatten()
            else:
                return self._X

    @X.setter
    def X(self, value: Optional[Union[np.ndarray, sparse.spmatrix]]):
        if value is None:
            if self.isview:
                raise ValueError('Copy the view before setting the data matrix to `None`.')
            if self.isbacked:
                raise ValueError('Not implemented.')
            self._X = None
            return
        var_get = self.n_vars == 1 and self.n_obs == len(value)
        obs_get = self.n_obs == 1 and self.n_vars == len(value)
        if var_get or obs_get or self.shape == value.shape:
            if self.isbacked:
                if self.isview:
                    self.file['X'][self._oidx, self._vidx] = value
                else:
                    self._set_backed('X', value)
            else:
                if self.isview:
                    # exit the view if we go from sparse to dense
                    if (
                        issparse(value) and not issparse(self._adata_ref._X)
                        or not issparse(value) and issparse(self._adata_ref._X)
                    ):
                        self._init_as_actual(self.copy())
                        self._X = value
                    else:
                        self._adata_ref._X[self._oidx, self._vidx] = value
                        self._init_X_as_view()
                else:
                    self._X = value
        else:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self.shape))

    @property
    def layers(self) -> LayersBase:
        """Dictionary-like object with values of the same dimensions as :attr:`X`.

        Layers in AnnData have API similar to loompy :ref:`loomlayers`.
        Return the layer named ``"unspliced"``::

            adata.layers["unspliced"]

        Create or replace the ``"spliced"`` layer::

            adata.layers["spliced"] = ...

        Assign the 10th column of layer ``"spliced"`` to the variable a::

            a = adata.layers["spliced"][:, 10]

        Delete the ``"spliced"`` layer::

            del adata.layers["spliced"]

        Return layers’ names::

            adata.layers.keys()

        .. warning::

           If AnnData is a view, setting subsets of layers modifies the original data.
        """
        return self._layers

    @property
    def raw(self) -> Raw:
        """Store raw version of :attr:`X` and :attr:`var` as ``.raw.X`` and ``.raw.var``.

        The :attr:`raw` attribute is initialized with the current content of an object by setting::

            adata.raw = adata

        Its content can be deleted by setting it back to ``None``::

            adata.raw = None

        Upon slicing an AnnData object along the observations (row) axis,
        :attr:`raw` is also sliced. Slicing an AnnData object along the variables
        (columns) axis, leaves :attr:`raw` unaffected. Note that you can call::

             adata.raw[:, 'orig_variable_name'].X

        to retrieve the data associated with a variable that might have been
        filtered out or "compressed away" in :attr:`X`.
        """
        return self._raw

    @raw.setter
    def raw(self, value: Optional['AnnData']):
        if not (isinstance(value, AnnData) or value is None):
            raise ValueError(
                'Can only init raw attribute with an AnnData object or `None`.')
        if value is None:
            self._raw = None
        else:
            if self.isview:
                self._init_as_actual(self.copy())
            self._raw = Raw(value)

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        return self._n_obs

    @property
    def n_vars(self) -> int:
        """Number of variables/features."""
        return self._n_vars

    @property
    def obs(self) -> pd.DataFrame:
        """One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    @obs.setter
    def obs(self, value: pd.DataFrame):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        utils.warn_no_string_index(value.index)
        if self.isview:
            self._init_as_actual(self.copy())
        self._obs = value

    @property
    def var(self) -> pd.DataFrame:
        """One-dimensional annotation of variables/ features (`pd.DataFrame`)."""
        return self._var

    @var.setter
    def var(self, value: pd.DataFrame):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_vars:
            raise ValueError('Length does not match.')
        utils.warn_no_string_index(value.index)
        if self.isview:
            self._init_as_actual(self.copy())
        self._var = value

    @property
    def uns(self) -> MutableMapping:
        """Unstructured annotation (ordered dictionary)."""
        return self._uns

    @uns.setter
    def uns(self, value: MutableMapping):
        if not isinstance(value, MutableMapping):
            raise ValueError('Only mutable mapping types (e.g. dict) are allowed for `.uns`.')
        if self.isview:
            self._init_as_actual(self.copy())
        self._uns = value

    @property
    def obsm(self) -> AxisArraysBase:
        """Multi-dimensional annotation of observations (mutable structured :class:`~numpy.ndarray`).

        Stores for each key, a two or higher-dimensional :class:`np.ndarray` of length
        ``n_obs``. Is sliced with ``data`` and ``obs`` but behaves otherwise like a
        :class:`dict`.
        """
        return self._obsm

    @obsm.setter
    def obsm(self, value):
        obsm = AxisArrays(self, 0, vals=convert_to_dict(value))
        if self.isview:
            self._init_as_actual(self.copy())
        self._obsm = obsm

    @property
    def varm(self) -> AxisArraysBase:
        """Multi-dimensional annotation of variables/ features (mutable structured :class:`~numpy.ndarray`).

        Stores for each key, a two or higher-dimensional :class:`~numpy.ndarray` of length
        ``n_vars``. Is sliced with ``data`` and ``var`` but behaves otherwise like a
        :class:`dict`.
        """
        return self._varm

    @varm.setter
    def varm(self, value):
        varm = AxisArrays(self, 1, vals=convert_to_dict(value))
        if self.isview:
            self._init_as_actual(self.copy())
        self._varm = varm

    # @property
    # def obsp(self) -> PairwiseArraysBase:
    #     """Pairwise annotation of observations, a mutable mapping with array-like values.

    #     Stores for each key, a two or higher-dimensional :class:`np.ndarray` whose
    #     first two dimensions are of length ``n_obs``. Is sliced with ``data`` and
    #     ``obs`` but behaves otherwise like a :class:`dict`.
    #     """
    #     return self._obsp

    # @obsp.setter
    # def obsp(self, value):
    #     obsp = PairwiseArrays(self, 0, vals=convert_to_dict(value))
    #     if self.isview:
    #         self._init_as_actual(self.copy())
    #     self._obsp = obsp

    # @property
    # def varp(self) -> PairwiseArraysBase:
    #     """Pairwise annotation of observations, a mutable mapping with array-like values.

    #     Stores for each key, a two or higher-dimensional :class:`np.ndarray` whose
    #     first two dimensions are of length ``n_var``. Is sliced with ``data`` and
    #     ``var`` but behaves otherwise like a :class:`dict`.
    #     """
    #     return self._varp

    # @varp.setter
    # def varp(self, value):
    #     varp = PairwiseArrays(self, 1, vals=convert_to_dict(value))
    #     if self.isview:
    #         self._init_as_actual(self.copy())
    #     self._varp = varp

    @property
    def obs_names(self) -> pd.Index:
        """Names of observations (alias for ``.obs.index``)."""
        return self.obs.index

    @obs_names.setter
    def obs_names(self, names: Sequence[str]):
        utils.warn_no_string_index(names)
        self._obs.index = names
        if not self._obs.index.is_unique:
            utils.warn_names_duplicates('obs')

    @property
    def var_names(self) -> pd.Index:
        """Names of variables (alias for ``.var.index``)."""
        return self._var.index

    @var_names.setter
    def var_names(self, names: Sequence[str]):
        utils.warn_no_string_index(names)
        self._var.index = names
        if not self._var.index.is_unique:
            utils.warn_names_duplicates('var')

    def obs_keys(self) -> List[str]:
        """List keys of observation annotation :attr:`obs`."""
        return self._obs.keys().tolist()

    def var_keys(self) -> List[str]:
        """List keys of variable annotation :attr:`var`."""
        return self._var.keys().tolist()

    def obsm_keys(self) -> List[str]:
        """List keys of observation annotation :attr:`obsm`."""
        return list(self._obsm.keys())

    def varm_keys(self) -> List[str]:
        """List keys of variable annotation :attr:`varm`."""
        return list(self._varm.keys())

    def uns_keys(self) -> List[str]:
        """List keys of unstructured annotation."""
        return sorted(list(self._uns.keys()))

    @property
    def isbacked(self) -> bool:
        """``True`` if object is backed on disk, ``False`` otherwise."""
        return self.filename is not None

    @property
    def isview(self) -> bool:
        """``True`` if object is view of another AnnData object, ``False`` otherwise."""
        return self._isview

    @property
    def filename(self) -> Optional[PathLike]:
        """Change to backing mode by setting the filename of a ``.h5ad`` file.

        - Setting the filename writes the stored data to disk.
        - Setting the filename when the filename was previously another name
          moves the backing file from the previous file to the new file. If you
          want to copy the previous file, use ``copy(filename='new_filename')``.
        """
        return self.file.filename

    @filename.setter
    def filename(self, filename: Optional[PathLike]):
        # convert early for later comparison
        filename = None if filename is None else Path(filename)
        # change from backing-mode back to full loading into memory
        if filename is None:
            if self.filename is not None:
                self.file._to_memory_mode()
            else:
                # both filename and self.filename are None
                # do nothing
                return
        else:
            if self.filename is not None:
                if self.filename != filename:
                    # write the content of self to the old file
                    # and close the file
                    self.write()
                    self.filename.rename(filename)
                else:
                    # do nothing
                    return
            else:
                # change from memory to backing-mode
                # write the content of self to disk
                self.write(filename, force_dense=True)
            # open new file for accessing
            self.file.open(filename, 'r+')
            # as the data is stored on disk, we can safely set self._X to None
            self._X = None

    def _set_backed(self, attr, value):
        if (not isinstance(self.file[attr], h5py.SparseDataset)
                and not issparse(value)):
            self.file[attr] = value
        else:
            del self.file[attr]
            self.file._file.create_dataset(attr, data=value)

    def _normalize_indices(self, index: Optional[Index]):
        # deal with tuples of length 1
        if isinstance(index, tuple) and len(index) == 1:
            index = index[0]
        # deal with pd.Series
        if isinstance(index, pd.Series):
            index = index.values
        if isinstance(index, tuple):
            if len(index) > 2:
                raise ValueError(
                    'AnnData can only be sliced in rows and columns.')
            # deal with pd.Series
            if isinstance(index[1], pd.Series):
                index = index[0], index[1].values
            if isinstance(index[0], pd.Series):
                index = index[0].values, index[1]

            no_slice = not any(isinstance(i, slice) for i in index)
            both_scalars = all(isinstance(i, (int, str, type(None))) for i in index)
            if no_slice and not both_scalars:
                raise NotImplementedError(
                    'Slicing with two indices at the same time is not yet implemented. '
                    'As a workaround, do row and column slicing succesively.')
            # Speed up and error prevention for boolean indices (Don’t convert to integer indices)
            # Needs to be refactored once we support a tuple of two arbitrary index types
            if any(isinstance(i, np.ndarray) and i.dtype == bool for i in index):
                return index
        obs, var = unpack_index(index)
        obs = _normalize_index(obs, self.obs_names)
        var = _normalize_index(var, self.var_names)
        return obs, var

    # TODO: this is not quite complete...
    def __delitem__(self, index: Index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if not self.isbacked:
            del self._X[obs, var]
        else:
            X = self.file['X']
            del X[obs, var]
            self._set_backed('X', X)
        if var == slice(None):
            del self._obs.iloc[obs, :]
        if obs == slice(None):
            del self._var.iloc[var, :]

    def __getitem__(self, index: Index) -> 'AnnData':
        """Returns a sliced view of the object."""
        return self._getitem_view(index)

    def _getitem_view(self, index: Index) -> 'AnnData':
        oidx, vidx = self._normalize_indices(index)
        return AnnData(self, oidx=oidx, vidx=vidx, asview=True)

    def _remove_unused_categories(self, df_full, df_sub, uns):
        from pandas.api.types import is_categorical
        for k in df_full:
            if is_categorical(df_full[k]):
                all_categories = df_full[k].cat.categories
                df_sub[k].cat.remove_unused_categories(inplace=True)
                # also correct the colors...
                if k + '_colors' in uns:
                    # this is a strange hack...
                    if np.array(uns[k + '_colors']).ndim == 0:
                        uns[k + '_colors'] = np.array(uns[k + '_colors'])[None]
                    uns[k + '_colors'] = np.array(uns[k + '_colors'])[
                        np.where(np.in1d(
                            all_categories, df_sub[k].cat.categories))[0]]

    def rename_categories(self, key: str, categories: Sequence[Any]):
        """Rename categories of annotation ``key`` in
        :attr:`obs`, :attr:`var` and :attr:`uns`.

        Only supports passing a list/array-like ``categories`` argument.

        Besides calling ``self.obs[key].cat.categories = categories`` -
        similar for :attr:`var` - this also renames categories in unstructured
        annotation that uses the categorical annotation ``key``.

        Parameters
        ----------
        key
             Key for observations or variables annotation.
        categories
             New categories, the same number as the old categories.
        """
        if isinstance(categories, Mapping):
            raise ValueError('Only list-like `categories` is supported.')
        if key in self.obs:
            old_categories = self.obs[key].cat.categories.tolist()
            self.obs[key].cat.rename_categories(categories, inplace=True)
        elif key in self.var:
            old_categories = self.var[key].cat.categories.tolist()
            self.var[key].cat.rename_categories(categories, inplace=True)
        else:
            raise ValueError('{} is neither in `.obs` nor in `.var`.'
                             .format(key))
        # this is not a good solution
        # but depends on the scanpy conventions for storing the categorical key
        # as `groupby` in the `params` slot
        for k1, v1 in self.uns.items():
            if isinstance(v1, Mapping):
                if 'params' in v1 and 'groupby' in v1['params']:
                    if v1['params']['groupby'] == key:
                        for k2, v2 in v1.items():
                            # picks out the recarrays that are named according to the old
                            # categories
                            if isinstance(v2, np.ndarray) and v2.dtype.names is not None:
                                if list(v2.dtype.names) == old_categories:
                                    self.uns[k1][k2].dtype.names = categories
                                else:
                                    logger.warning(
                                        'Omitting {}/{} as old categories do not match.'
                                        .format(k1, k2))

    def strings_to_categoricals(self, df: Optional[pd.DataFrame] = None):
        """Transform string annotations to categoricals.

        Only affects string annotations that lead to less categories than the
        total number of observations.

        Params
        ------
        df
            If ``df`` is ``None``, modifies both :attr:`obs` and :attr:`.var`,
            otherwise modifies ``df`` inplace.

        Notes
        -----
        Turns the view of an :class:`~anndata.AnnData` into an actual
        :class:`~anndata.AnnData`.
        """
        dont_modify = False  # only necessary for backed views
        if df is None:
            dfs = [self.obs, self.var]
            if self.isview:
                if not self.isbacked:
                    self._init_as_actual(self.copy())
                else:
                    dont_modify = True
        else:
            dfs = [df]
        for df in dfs:
            string_cols = [
                key for key in df.columns
                if is_string_dtype(df[key])
                and not is_categorical(df[key])
            ]
            for key in string_cols:
                # make sure we only have strings (could be that there are
                # np.nans (float), -666, '-666', for instance)
                c = df[key].astype('U')
                # make a categorical
                c = pd.Categorical(c, categories=natsorted(np.unique(c)))
                if len(c.categories) < len(c):
                    if dont_modify:
                        raise RuntimeError(
                            'Please call `.strings_to_categoricals()` on full AnnData, not on this view. '
                            'You might encounter this error message while copying or writing to disk.'
                        )
                    df[key] = c
                    logger.info('... storing {!r} as categorical'.format(key))

    _sanitize = strings_to_categoricals  # backwards compat

    def _slice_uns_sparse_matrices_inplace(self, uns, oidx):
        # slice sparse spatrices of n_obs × n_obs in self.uns
        if not (isinstance(oidx, slice) and
                oidx.start is None and oidx.step is None and oidx.stop is None):
            for k, v in uns.items():
                # treat nested dicts
                if isinstance(v, Mapping):
                    self._slice_uns_sparse_matrices_inplace(v, oidx)
                if isinstance(v, sparse.spmatrix) and v.shape == (
                        self.n_obs, self.n_obs):
                    uns[k] = v.tocsc()[:, oidx].tocsr()[oidx, :]

    def _inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as ``adata = adata[:, index]``, but inplace.
        """
        adata_subset = self[:, index].copy()
        self._init_as_actual(adata_subset, dtype=self._X.dtype)

    def _inplace_subset_obs(self, index):
        """Inplace subsetting along variables dimension.

        Same as ``adata = adata[index, :]``, but inplace.
        """
        adata_subset = self[index].copy()
        self._init_as_actual(adata_subset, dtype=self._X.dtype)

    def __setitem__(self, index: Index, val: Union[int, float, np.ndarray, sparse.spmatrix]):
        if self.isview:
            raise ValueError('Object is view and cannot be accessed with `[]`.')
        obs, var = self._normalize_indices(index)
        if not self.isbacked:
            self._X[obs, var] = val
        else:
            X = self.file['X']
            X[obs, var] = val
            self._set_backed('X', X)

    def __len__(self) -> int:
        return self.shape[0]

    def transpose(self) -> 'AnnData':
        """Transpose whole object.

        Data matrix is transposed, observations and variables are interchanged.
        """
        if not self.isbacked: X = self._X
        else: X = self.file['X']
        if self.isview:
            raise ValueError(
                'You\'re trying to transpose a view of an `AnnData`, which is currently not implemented. '
                'Call `.copy()` before transposing.')
        layers = {k:(v.T.tocsr() if sparse.isspmatrix_csr(v) else v.T) for (k, v) in self.layers.items()}
        if sparse.isspmatrix_csr(X):
            return AnnData(X.T.tocsr(), self._var, self._obs, self._uns,
                           self._varm.flipped(), self._obsm.flipped(),
                           filename=self.filename, layers=layers, dtype=self.X.dtype.name)
        return AnnData(X.T, self._var, self._obs, self._uns,
                       self._varm.flipped(), self._obsm.flipped(),
                       filename=self.filename, layers=layers, dtype=self.X.dtype.name)

    T = property(transpose)

    def to_df(self) -> pd.DataFrame:
        """Generate shallow :class:`~pandas.DataFrame`.

        The data matrix :attr:`X` is returned as
        :class:`~pandas.DataFrame`, where :attr:`obs_names` initializes the
        index, and :attr:`var_names` the columns.

        * No annotations are maintained in the returned object.
        * The data matrix is densified in case it is sparse.
        """
        if issparse(self._X):
            X = self._X.toarray()
        else:
            X = self._X
        return pd.DataFrame(X, index=self.obs_names, columns=self.var_names)

    def _get_X(self, use_raw=False, layer=None):
        """
        Convenience method for getting expression values with common arguments and error handling.
        """
        is_layer = layer is not None
        if use_raw and is_layer:
            raise ValueError(
                "Cannot use expression from both layer and raw. You provided:"
                "'use_raw={}' and 'layer={}'".format(use_raw, layer)
            )
        if is_layer:
            return self.layers[layer]
        elif use_raw:
            if self.raw is None:
                raise ValueError("This AnnData doesn't have a value in `.raw`.")
            return self.raw.X
        else:
            return self.X

    def obs_vector(
        self, k: str, *, layer: Optional[str] = None
    ) -> np.ndarray:
        """
        Convenience function for returning a 1 dimensional ndarray of values
        from `.X`, `.layers[k]`, or `.obs`.

        Made for convenience, not performance. Intentionally permissive about
        arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in `.var_names` or `.obs.columns`.
        layer
            What layer values should be returned from. If `None`, `.X` is used.

        Returns
        -------
        A one dimensional nd array, with values for each obs in the same order
        as `.obs_names`.
        """
        if layer == "X":
            if "X" in self.layers:
                pass
            else:
                warnings.warn(
                    "In a future version of AnnData, access to `.X` by passing"
                    " `layer='X'` will be removed. Instead pass `layer=None`.",
                    FutureWarning
                )
                layer = None

        if k in self.obs:
            return self.obs[k].values
        else:
            idx = self._normalize_indices((slice(None), k))
            a = self._get_X(layer=layer)[idx]
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)

    def var_vector(
        self, k, *, layer: Optional[str] = None
    ) -> np.ndarray:
        """
        Convenience function for returning a 1 dimensional ndarray of values
        from `.X`, `.layers[k]`, or `.obs`.

        Made for convenience, not performance. Intentionally permissive about
        arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in `.obs_names` or `.var.columns`.
        layer
            What layer values should be returned from. If `None`, `.X` is used.

        Returns
        -------
        A one dimensional nd array, with values for each var in the same order
        as `.var_names`.
        """
        if layer == "X":
            if "X" in self.layers:
                pass
            else:
                warnings.warn(
                    "In a future version of AnnData, access to `.X` by passing"
                    " `layer='X'` will be removed. Instead pass `layer=None`.",
                    FutureWarning
                )
                layer = None

        if k in self.var:
            return self.var[k].values
        else:
            idx = self._normalize_indices((k, slice(None)))
            a = self._get_X(layer=layer)[idx]
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)

    @utils.deprecated("obs_vector")
    def _get_obs_array(self, k, use_raw=False, layer=None):
        """Get an array from the layer (default layer='X') along the observation dimension by first looking up
        obs.keys and then var.index."""
        if not use_raw or k in self.obs.columns:
            return self.obs_vector(k=k, layer=layer)
        else:
            return self.raw.obs_vector(k)

    @utils.deprecated("var_vector")
    def _get_var_array(self, k, use_raw=False, layer=None):
        """Get an array from the layer (default layer='X') along the variables dimension by first looking up
        ``var.keys`` and then ``obs.index``."""
        if not use_raw or k in self.var.columns:
            return self.var_vector(k=k, layer=layer)
        else:
            return self.raw.var_vector(k)

    def copy(self, filename: Optional[PathLike] = None) -> 'AnnData':
        """Full copy, optionally on disk."""
        if not self.isbacked:
            return AnnData(self._X.copy() if self._X is not None else None,
                           self._obs.copy(),
                           self._var.copy(),
                           # deepcopy on DictView does not work and is unnecessary
                           # as uns was copied already before
                           self._uns.copy() if isinstance(self._uns, DictView) else deepcopy(self._uns),
                           self._obsm.copy(), self._varm.copy(),
                           raw=None if self._raw is None else self._raw.copy(),
                           layers=dict(self.layers),
                           dtype=self._X.dtype.name if self._X is not None else 'float32')
        else:
            if filename is None:
                raise ValueError(
                    'To copy an AnnData object in backed mode, '
                    'pass a filename: `.copy(filename=\'myfilename.h5ad\')`.')
            if self.isview:
                self.write(filename)
            else:
                from shutil import copyfile
                copyfile(self.filename, filename)
            return AnnData(filename=filename)

    def concatenate(
        self, *adatas: 'AnnData',
        join: str = 'inner',
        batch_key: str = 'batch',
        batch_categories: Sequence[Any] = None,
        index_unique: Optional[str] = '-'
    ) -> 'AnnData':
        """Concatenate along the observations axis.

        The :attr:`uns`, :attr:`varm` and :attr:`obsm` attributes are ignored.

        Currently, this works only in ``'memory'`` mode.

        Parameters
        ----------
        adatas
            AnnData matrices to concatenate with. Each matrix is referred to as
            a “batch”.
        join
            Use intersection (``'inner'``) or union (``'outer'``) of variables.
        batch_key
            Add the batch annotation to :attr:`obs` using this key.
        batch_categories
            Use these as categories for the batch annotation. By default, use increasing numbers.
        index_unique
            Make the index unique by joining the existing index names with the
            batch category, using ``index_unique='-'``, for instance. Provide
            ``None`` to keep existing indices.

        Returns
        -------
        :class:`~anndata.AnnData`
            The concatenated :class:`~anndata.AnnData`, where ``adata.obs[batch_key]``
            stores a categorical variable labeling the batch.

        Notes
        -----

        .. warning::

           If you use ``join='outer'`` this fills 0s for sparse data when
           variables are absent in a batch. Use this with care. Dense data is
           filled with ``NaN``. See the examples.

        Examples
        --------
        Joining on intersection of variables.

        >>> adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['s1', 's2'],
        >>>                   'anno1': ['c1', 'c2']},
        >>>                  {'var_names': ['a', 'b', 'c'],
        >>>                   'annoA': [0, 1, 2]})
        >>> adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['s3', 's4'],
        >>>                   'anno1': ['c3', 'c4']},
        >>>                  {'var_names': ['d', 'c', 'b'],
        >>>                   'annoA': [0, 1, 2]})
        >>> adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['s1', 's2'],
        >>>                   'anno2': ['d3', 'd4']},
        >>>                  {'var_names': ['d', 'c', 'b'],
        >>>                   'annoA': [0, 2, 3],
        >>>                   'annoB': [0, 1, 2]})
        >>>
        >>> adata = adata1.concatenate(adata2, adata3)
        >>> adata
        AnnData object with n_obs × n_vars = 6 × 2
            obs_keys = ['anno1', 'anno2', 'batch']
            var_keys = ['annoA-0', 'annoA-1', 'annoB-2', 'annoA-2']
        >>> adata.X
        array([[2., 3.],
               [5., 6.],
               [3., 2.],
               [6., 5.],
               [3., 2.],
               [6., 5.]], dtype=float32)
        >>> adata.obs
             anno1 anno2 batch
        s1-0    c1   NaN     0
        s2-0    c2   NaN     0
        s3-1    c3   NaN     1
        s4-1    c4   NaN     1
        s1-2   NaN    d3     2
        s2-2   NaN    d4     2
        >>> adata.var.T
                 b  c
        annoA-0  1  2
        annoA-1  2  1
        annoB-2  2  1
        annoA-2  3  2

        Joining on the union of variables.

        >>> adata = adata1.concatenate(adata2, adata3, join='outer')
        >>> adata
        AnnData object with n_obs × n_vars = 6 × 4
            obs_keys = ['anno1', 'anno2', 'batch']
            var_keys = ['annoA-0', 'annoA-1', 'annoB-2', 'annoA-2']
        >>> adata.var.T
        index      a    b    c    d
        annoA-0  0.0  1.0  2.0  NaN
        annoA-1  NaN  2.0  1.0  0.0
        annoB-2  NaN  2.0  1.0  0.0
        annoA-2  NaN  3.0  2.0  0.0
        >>> adata.var_names
        Index(['a', 'b', 'c', 'd'], dtype='object')
        >>> adata.X
        array([[ 1.,  2.,  3., nan],
               [ 4.,  5.,  6., nan],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.]], dtype=float32)
        >>> adata.X.sum(axis=0)
        array([nan, 25., 23., nan], dtype=float32)
        >>> import pandas as pd
        >>> Xdf = pd.DataFrame(adata.X, columns=adata.var_names)
        index    a    b    c    d
        0      1.0  2.0  3.0  NaN
        1      4.0  5.0  6.0  NaN
        2      NaN  3.0  2.0  1.0
        3      NaN  6.0  5.0  4.0
        4      NaN  3.0  2.0  1.0
        5      NaN  6.0  5.0  4.0
        >>> Xdf.sum()
        index
        a     5.0
        b    25.0
        c    23.0
        d    10.0
        dtype: float32
        >>> from numpy import ma
        >>> adata.X = ma.masked_invalid(adata.X)
        >>> adata.X
        masked_array(
          data=[[1.0, 2.0, 3.0, --],
                [4.0, 5.0, 6.0, --],
                [--, 3.0, 2.0, 1.0],
                [--, 6.0, 5.0, 4.0],
                [--, 3.0, 2.0, 1.0],
                [--, 6.0, 5.0, 4.0]],
          mask=[[False, False, False,  True],
                [False, False, False,  True],
                [ True, False, False, False],
                [ True, False, False, False],
                [ True, False, False, False],
                [ True, False, False, False]],
          fill_value=1e+20,
          dtype=float32)
        >>> adata.X.sum(axis=0).data
        array([ 5., 25., 23., 10.], dtype=float32)

        The masked array is not saved but has to be reinstantiated after saving.

        >>> adata.write('./test.h5ad')
        >>> from anndata import read_h5ad
        >>> adata = read_h5ad('./test.h5ad')
        >>> adata.X
        array([[ 1.,  2.,  3., nan],
               [ 4.,  5.,  6., nan],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.]], dtype=float32)

        For sparse data, everything behaves similarly, except that for ``join='outer'``, zeros are added.

        >>> from scipy.sparse import csr_matrix
        >>> adata1 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
        >>>                  {'obs_names': ['s1', 's2'],
        >>>                   'anno1': ['c1', 'c2']},
        >>>                  {'var_names': ['a', 'b', 'c']})
        >>> adata2 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
        >>>                  {'obs_names': ['s3', 's4'],
        >>>                   'anno1': ['c3', 'c4']},
        >>>                  {'var_names': ['d', 'c', 'b']})
        >>> adata3 = AnnData(csr_matrix([[1, 2, 0], [0, 5, 6]]),
        >>>                  {'obs_names': ['s5', 's6'],
        >>>                   'anno2': ['d3', 'd4']},
        >>>                  {'var_names': ['d', 'c', 'b']})
        >>>
        >>> adata = adata1.concatenate(adata2, adata3, join='outer')
        >>> adata.var_names
        Index(['a', 'b', 'c', 'd'], dtype='object')
        >>> adata.X.toarray()
        array([[0., 2., 3., 0.],
               [0., 5., 6., 0.],
               [0., 3., 2., 0.],
               [0., 6., 5., 0.],
               [0., 0., 2., 1.],
               [0., 6., 5., 0.]], dtype=float32)
        """
        if self.isbacked:
            raise ValueError(
                'Currently, concatenate does only work in \'memory\' mode.')

        if len(adatas) == 0:
            return self
        elif len(adatas) == 1 and not isinstance(adatas[0], AnnData):
            adatas = adatas[0]  # backwards compatibility
        all_adatas = (self,) + tuple(adatas)

        # for controlled behavior, make all variable names unique
        printed_info = False
        for i, ad in enumerate(all_adatas):
            if not ad.var_names.is_unique:
                ad.var_names = utils.make_index_unique(ad.var_names)
                if not printed_info:
                    logger.info(
                        'Making variable names unique for controlled concatenation.')
                    printed_info = True

        # define variable names of joint AnnData
        mergers = dict(inner=set.intersection, outer=set.union)
        var_names_reduce = reduce(mergers[join], (set(ad.var_names) for ad in all_adatas))
        # restore order of initial var_names, append non-sortable names at the end
        # see how this was done in the repo at commit state
        # 40a24f
        var_names = []
        for v in all_adatas[0].var_names:
            if v in var_names_reduce:
                var_names.append(v)
                var_names_reduce.remove(v)  # update the set
        var_names = pd.Index(var_names + list(var_names_reduce))

        if batch_categories is None:
            categories = [str(i) for i, _ in enumerate(all_adatas)]
        elif len(batch_categories) == len(all_adatas):
            categories = batch_categories
        else:
            raise ValueError('Provide as many `batch_categories` as `adatas`.')

        out_shape = (sum(a.n_obs for a in all_adatas), len(var_names))

        any_sparse = any(issparse(a.X) for a in all_adatas)
        if join == 'outer':
            if any_sparse:  # not sure whether the lil_matrix is really the best option
                X = sparse.lil_matrix(out_shape, dtype=self.X.dtype)
            else:
                X = np.empty(out_shape, dtype=self.X.dtype)
                X[:] = np.nan
        else:
            Xs = []

        # create layers dict that contains layers shared among all AnnDatas
        layers = OrderedDict()
        shared_layers = [key for key in all_adatas[0].layers.keys()
                         if all([key in ad.layers.keys() for ad in all_adatas])]
        for key in shared_layers:
            layers[key] = []

        # check whether tries to do 'outer' join and layers is non_empty.
        if join == 'outer' and len(shared_layers) > 0:
            logger.info(
                'layers concatenation is not yet available for \'outer\' intersection and will be ignored.')

        # check whether layers are not consistently set in all AnnData objects.
        n_layers = np.array([len(ad.layers.keys()) for ad in all_adatas])
        if join == 'inner' and not all(len(shared_layers) == n_layers):
            logger.info(
                'layers are inconsistent - only layers that are shared among all AnnData objects are included.')

        var = pd.DataFrame(index=var_names)

        if join == 'inner':
            ad_ref = all_adatas[0]
            cols_intersect = set(ad_ref.var.columns)
            for ad in all_adatas[1:]:
                cols_intersect &= set(ad.var.columns)
                cols_intersect = {
                    col for col in cols_intersect
                    if ad_ref.var.loc[var_names, col].equals(ad.var.loc[var_names, col])
                }
                if not cols_intersect: break

        obs_i = 0  # start of next adata’s observations in X
        out_obss = []
        for i, ad in enumerate(all_adatas):
            if join == 'outer':
                # only those names that are actually present in the current AnnData
                vars_intersect = [v for v in var_names if v in ad.var_names]
            else:
                vars_intersect = var_names

            # X
            if join == 'outer':
                # this is pretty slow, I guess sparse matrices shouldn't be
                # constructed like that
                X[obs_i:obs_i+ad.n_obs,
                  var_names.isin(vars_intersect)] = ad[:, vars_intersect].X
            else:
                Xs.append(ad[:, vars_intersect].X)
            obs_i += ad.n_obs

            # layers
            if join == 'inner':
                for key in shared_layers:
                    layers[key].append(ad[:, vars_intersect].layers[key])

            # obs
            obs = ad.obs.copy()
            obs[batch_key] = pd.Categorical(ad.n_obs * [categories[i]], categories)
            if (is_string_dtype(all_adatas[0].obs.index) and not
                is_string_dtype(ad.obs.index)):
                obs.index = obs.index.astype(str)
            if index_unique is not None:
                if not is_string_dtype(ad.obs.index):
                    obs.index = obs.index.astype(str)
                obs.index = obs.index.values + index_unique + categories[i]
            out_obss.append(obs)

            # var
            for c in ad.var.columns:
                if join == 'inner' and c in cols_intersect:
                    if c not in var.columns:
                        var.loc[vars_intersect, c] = ad.var.loc[vars_intersect, c]
                    continue
                new_c = c + (index_unique if index_unique is not None else '-') + categories[i]
                var.loc[vars_intersect, new_c] = ad.var.loc[vars_intersect, c]

        if join == 'inner':
            if any_sparse:
                from scipy.sparse import vstack
                X = vstack(Xs)
            else:
                X = np.concatenate(Xs)

            for key in shared_layers:
                if any(issparse(a.layers[key]) for a in all_adatas):
                    layers[key] = vstack(layers[key])
                else:
                    layers[key] = np.concatenate(layers[key])

        obs = pd.concat(out_obss, sort=True)

        if any_sparse:
            sparse_format = all_adatas[0].X.getformat()
            X = X.asformat(sparse_format)

        new_adata = AnnData(X, obs, var, layers=layers) if join == 'inner' else AnnData(X, obs, var)
        if not obs.index.is_unique:
            logger.info(
                'Or pass `index_unique!=None` to `.concatenate`.')
        return new_adata

    def var_names_make_unique(self, join: str = '-'):
        self.var.index = utils.make_index_unique(self.var.index, join)

    var_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def obs_names_make_unique(self, join: str = '-'):
        self.obs.index = utils.make_index_unique(self.obs.index, join)

    obs_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def _check_uniqueness(self):
        if not self.obs.index.is_unique:
            utils.warn_names_duplicates('obs')
        if not self.var.index.is_unique:
            utils.warn_names_duplicates('var')

    def __contains__(self, key: Any):
        raise AttributeError('AnnData has no attribute __contains__, '
                             'don\'t check `in adata`.')

    def _check_dimensions(self, key=None):
        if key is None:
            key = {'obs', 'var', 'obsm', 'varm'}
        else:
            key = {key}
        if 'obs' in key and len(self._obs) != self._n_obs:
            raise ValueError('Observations annot. `obs` must have number of '
                             'rows of `X` ({}), but has {} rows.'
                             .format(self._n_obs, self._obs.shape[0]))
        if 'var' in key and len(self._var) != self._n_vars:
            raise ValueError('Variables annot. `var` must have number of '
                             'columns of `X` ({}), but has {} rows.'
                             .format(self._n_vars, self._var.shape[0]))
        if 'obsm' in key:
            obsm = self._obsm
            if not all([o.shape[0] == self._n_obs for o in obsm.values()]) and len(obsm.dim_names) != self._n_obs:
                raise ValueError('Observations annot. `obsm` must have number of '
                                'rows of `X` ({}), but has {} rows.'
                                .format(self._n_obs, len(self._obsm)))
        if 'varm' in key:
            varm = self._varm
            if not all([v.shape[0] == self._n_vars for v in varm.values()]) and len(varm.dim_names) != self._n_var:
                raise ValueError('Variables annot. `varm` must have number of '
                                'columns of `X` ({}), but has {} rows.'
                                .format(self._n_vars, len(self._varm)))

    def write_h5ad(
        self,
        filename: Optional[PathLike] = None,
        compression: Optional[str] = None,
        compression_opts: Union[int, Any] = None,
        force_dense: Optional[bool] = None
    ):
        """Write ``.h5ad``-formatted hdf5 file.

        .. note::

            Setting compression to ``'gzip'`` can save disk space but
            will slow down writing and subsequent reading. Prior to
            v0.6.16, this was the default for parameter
            ``compression``.

        Generally, if you have sparse data that are stored as a dense
        matrix, you can dramatically improve performance and reduce
        disk space by converting to a :class:`~scipy.sparse.csr_matrix`::

            from scipy.sparse import csr_matrix
            adata.X = csr_matrix(adata.X)

        Parameters
        ----------
        filename
            Filename of data file. Defaults to backing file.
        compression : ``None``,  {``'gzip'``, ``'lzf'``} (default: ``None``)
            See the h5py :ref:`dataset_compression`.
        compression_opts
            See the h5py :ref:`dataset_compression`.
        force_dense
            Write sparse data as a dense matrix. Defaults to ``True`` if object is
            backed, otherwise to ``False``.
        """
        from ..readwrite.write import _write_h5ad

        if filename is None and not self.isbacked:
            raise ValueError('Provide a filename!')
        if filename is None:
            filename = self.filename
        if force_dense is None:
            force_dense = self.isbacked

        _write_h5ad(filename, self, compression=compression,
                    compression_opts=compression_opts, force_dense=force_dense)

        if self.isbacked:
            self.file.close()

    write = write_h5ad  # a shortcut and backwards compat

    def write_csvs(self, dirname: PathLike, skip_data: bool = True, sep: str = ','):
        """Write annotation to ``.csv`` files.

        It is not possible to recover the full :class:`~anndata.AnnData` from the
        output of this function. Use :meth:`~anndata.AnnData.write` for this.

        Parameters
        ----------
        dirname
            Name of directory to which to export.
        skip_data
             Skip the data matrix :attr:`X`.
        sep
             Separator for the data.
        """
        from ..readwrite.write import write_csvs
        write_csvs(dirname, self, skip_data=skip_data, sep=sep)

    def write_loom(self, filename: PathLike, write_obsm_varm: bool = False):
        """Write ``.loom``-formatted hdf5 file.

        Parameters
        ----------
        filename
            The filename.
        """
        from ..readwrite.write import write_loom
        write_loom(filename, self, write_obsm_varm = write_obsm_varm)

    def write_zarr(
        self,
        store: Union[MutableMapping, PathLike],
        chunks: Union[bool, int, Tuple[int, ...]],
    ):
        """Write a hierarchical Zarr array store.

        Parameters
        ----------
        store
            The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
        chunks
            Chunk shape.
        """
        from ..readwrite.write import write_zarr
        write_zarr(store, self, chunks=chunks)

    def chunked_X(self, chunk_size: Optional[int] = None):
        """Return an iterator over the rows of the data matrix :attr:`X`.

        Parameters
        ----------
        chunk_size
            Row size of a single chunk.
        """
        if chunk_size is None:
            # Should be some adaptive code
            chunk_size = 6000
        start = 0
        n = self.n_obs
        for _ in range(int(n // chunk_size)):
            end = start + chunk_size
            yield (self.X[start:end], start, end)
            start = end
        if start < n:
            yield (self.X[start:n], start, n)

    def chunk_X(
        self,
        select: Union[int, List[int], Tuple[int, ...], np.ndarray] = 1000,
        replace: bool = True,
    ):
        """Return a chunk of the data matrix :attr:`X` with random or specified indices.

        Parameters
        ----------
        select
            If select is an integer, a random chunk of row size = select will be returned.
            If select is a list, tuple or numpy array of integers, then a chunk
            with these indices will be returned.

        replace
            If select is an integer then ``replace=True`` specifies random sampling of indices
            with replacement, ``replace=False`` - without replacement.
        """
        if isinstance(select, int):
            select = select if select < self.n_obs else self.n_obs
            choice = np.random.choice(self.n_obs, select, replace)
        elif isinstance(select, (np.ndarray, list, tuple)):
            choice = np.asarray(select)
        else:
            raise ValueError('select should be int or array')

        reverse = None
        if self.isbacked:
            # h5py can only slice with a sorted list of unique index values
            # so random batch with indices [2, 2, 5, 3, 8, 10, 8] will fail
            # this fixes the problem
            indices, reverse = np.unique(choice, return_inverse=True)
            selection = self.X[indices.tolist()]
        else:
            selection = self.X[choice]

        selection = selection.toarray() if issparse(selection) else selection
        return selection if reverse is None else selection[reverse]

    @staticmethod
    def _args_from_dict(ddata: Mapping[str, Any]):
        """Allows to construct an instance of AnnData from a dictionary.

        Acts as interface for the communication with the hdf5 file.

        In particular, from a dict that has been written using
        ``AnnData._to_dict_fixed_width_arrays``.
        """
        d_true_keys = {}
        # backwards compat
        uns_is_not_key = False
        valid_keys = []
        for keys in AnnData._H5_ALIASES.values():
            valid_keys += keys
        valid_keys += ['raw.X', 'raw.var', 'raw.varm', 'raw.cat']
        for key in ddata.keys():
            # if there is another key then the prdedefined
            # then we are reading the old format
            if key not in valid_keys:
                uns_is_not_key = True

        for true_key, keys in AnnData._H5_ALIASES.items():
            for key in keys:
                if key in ddata:
                    d_true_keys[true_key] = ddata[key]
                    if uns_is_not_key: del ddata[key]
                    break
            else:
                d_true_keys[true_key] = None

        # transform recarray to dataframe
        for true_key, keys in AnnData._H5_ALIASES_NAMES.items():
            if d_true_keys[true_key] is not None:
                for key in keys:
                    if key in d_true_keys[true_key].dtype.names:
                        d_true_keys[true_key] = pd.DataFrame.from_records(
                            d_true_keys[true_key], index=key)
                        break
                d_true_keys[true_key].index = d_true_keys[true_key].index.astype('U')
                # transform to unicode string
                # TODO: this is quite a hack
                for c in d_true_keys[true_key].columns:
                    if is_string_dtype(d_true_keys[true_key][c]):
                        d_true_keys[true_key][c] = pd.Index(
                            d_true_keys[true_key][c]).astype('U').values

        # these are the category fields
        k_to_delete = []
        items = (
            ddata.items() if uns_is_not_key
            else ddata['uns'].items() if 'uns' in ddata else []
        )
        for k, v in items:
            if k.endswith('_categories'):
                k_stripped = k.replace('_categories', '')
                if isinstance(v, (str, int)):  # fix categories with a single category
                    v = [v]
                for ann in ['obs', 'var']:
                    if k_stripped in d_true_keys[ann]:
                        d_true_keys[ann][k_stripped] = pd.Categorical.from_codes(
                            codes=d_true_keys[ann][k_stripped].values,
                            categories=v,
                        )
                k_to_delete.append(k)

        for k in k_to_delete:
            if uns_is_not_key:
                del ddata[k]
            else:
                del ddata['uns'][k]

        # assign the variables
        X = d_true_keys['X']
        obs = d_true_keys['obs']
        obsm = d_true_keys['obsm']
        var = d_true_keys['var']
        varm = d_true_keys['varm']
        layers = d_true_keys['layers']

        raw = None
        if 'raw.X' in ddata:
            raw = {}
            raw['X'] = ddata['raw.X']
            del ddata['raw.X']
            # get the dataframe
            raw['var'] = pd.DataFrame.from_records(
                ddata['raw.var'], index='index')
            del ddata['raw.var']
            raw['var'].index = raw['var'].index.astype('U')
            # transform to unicode string
            for c in raw['var'].columns:
                if is_string_dtype(raw['var'][c]):
                    raw['var'][c] = pd.Index(raw['var'][c]).astype('U').values
            # these are the category fields
            if 'raw.cat' in ddata:  # old h5ad didn't have that field
                for k, v in ddata['raw.cat'].items():
                    if k.endswith('_categories'):
                        k_stripped = k.replace('_categories', '')
                        if isinstance(v, (str, int)):  # fix categories with a single category
                            v = [v]
                        raw['var'][k_stripped] = pd.Categorical.from_codes(
                            codes=raw['var'][k_stripped].values,
                            categories=v)
                del ddata['raw.cat']
        if 'raw.varm' in ddata:
            raw['varm'] = ddata['raw.varm']
            del ddata['raw.varm']
        elif raw is not None:
            raw['varm'] = None

        # the remaining fields are the unstructured annotation
        uns = (
            ddata if uns_is_not_key
            else ddata['uns'] if 'uns' in ddata
            else {}
        )

        return X, obs, var, uns, obsm, varm, layers, raw

    def _to_dict_fixed_width_arrays(self, var_len_str=True):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        self.strings_to_categoricals()
        obs_rec, uns_obs = df_to_records_fixed_width(self._obs, var_len_str)
        var_rec, uns_var = df_to_records_fixed_width(self._var, var_len_str)
        layers = dict(self.layers)
        d = {
            'X': self._X,
            'obs': obs_rec,
            'var': var_rec,
            'obsm': self._obsm,
            'varm': self._varm,
            'layers': layers,
            # add the categories to the unstructured annotation
            'uns': {**self._uns, **uns_obs, **uns_var}}

        if self.raw is not None:
            self.strings_to_categoricals(self.raw._var)
            var_rec, uns_var = df_to_records_fixed_width(self.raw._var, var_len_str)
            d['raw.X'] = self.raw.X
            d['raw.var'] = var_rec
            d['raw.varm'] = self.raw.varm
            d['raw.cat'] = uns_var

        return d

    # --------------------------------------------------------------------------
    # all of the following is for backwards compat
    # --------------------------------------------------------------------------

    @property
    @utils.deprecated('X')
    def data(self):
        return self.X

    @data.setter
    @utils.deprecated('X')
    def data(self, value):
        self.X = value

    @property
    @utils.deprecated('n_obs')
    def n_smps(self):
        return self.n_obs

    @property
    @utils.deprecated('obs')
    def smp(self):
        return self.obs

    @smp.setter
    @utils.deprecated('obs')
    def smp(self, value):
        self.obs = value

    @property
    @utils.deprecated('uns')
    def add(self):
        return self.uns

    @add.setter
    @utils.deprecated('uns')
    def add(self, value):
        self.uns = value

    @property
    @utils.deprecated('obsm')
    def smpm(self):
        return self.obsm

    @smpm.setter
    @utils.deprecated('obsm')
    def smpm(self, value):
        self.obsm = value

    @property
    @utils.deprecated('obs_names')
    def smp_names(self):
        return self.obs_names

    @smp_names.setter
    @utils.deprecated('obs_names')
    def smp_names(self, names):
        self.obs_names = names

    @utils.deprecated('obs_keys')
    def smp_keys(self):
        return self.obs_keys()

    @utils.deprecated('obsm_keys')
    def smpm_keys(self):
        return self.obsm_keys()

    def _clean_up_old_format(self, uns):
        # multicolumn keys
        # all of the rest is only for backwards compat
        if uns and '_obs_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['_obs_keys_multicol'])
            del uns['_obs_keys_multicol']
        elif uns and 'obs_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['obs_keys_multicol'])
            del uns['obs_keys_multicol']
        elif uns and '_smp_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['_smp_keys_multicol'])
            del uns['_smp_keys_multicol']
        elif uns and 'smp_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['smp_keys_multicol'])
            del uns['smp_keys_multicol']
        else:
            _keys_multicol_obs = []
        if uns and '_var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['_var_keys_multicol'])
            del uns['_var_keys_multicol']
        elif uns and 'var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['var_keys_multicol'])
            del uns['var_keys_multicol']
        else:
            _keys_multicol_var = []

        # now, for compat, fill the old multicolumn entries into obsm and varm
        # and remove them from obs and var
        for key in _keys_multicol_obs:
            self._obsm[key] = self._get_multicol_field_obs(key)
        for key in _keys_multicol_var:
            self._varm[key] = self._get_multicol_field_var(key)

    def _get_multicol_field_obs(self, key_multicol):
        return self._get_and_delete_multicol_field('obs', key_multicol)

    def _get_multicol_field_var(self, key_multicol):
        return self._get_and_delete_multicol_field('var', key_multicol)

    def _get_and_delete_multicol_field(self, a, key_multicol):
        keys = []
        for k in getattr(self, a).columns:
            if k.startswith(key_multicol):
                keys.append(k)
        values = getattr(self, a)[keys].values
        getattr(self, a).drop(keys, axis=1, inplace=True)
        return values
