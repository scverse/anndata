"""Main class and helper functions.
"""
import os, sys
import warnings
import logging as logg
from enum import Enum
from collections import Mapping, Sequence, Sized, ChainMap
from functools import reduce
from typing import Union

import numpy as np
from numpy import ma
import pandas as pd
from pandas.core.index import RangeIndex
from pandas.api.types import is_string_dtype, is_categorical
from scipy import sparse
from scipy.sparse import issparse
from scipy.sparse.sputils import IndexMixin
from textwrap import dedent
from natsort import natsorted

from . import h5py
from . import utils

# FORMAT = '%(levelname)s: %(message)s'  # TODO: add a better formatter
FORMAT = '%(message)s'
logg.basicConfig(format=FORMAT, level=logg.INFO, stream=sys.stdout)

_MAIN_NARRATIVE = """\
:class:`~anndata.AnnData` stores a data matrix ``.X`` together with
annotations of observations ``.obs``, variables ``.var`` and
unstructured annotations ``.uns``.

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg"
         style="width: 350px; margin: 10px 0px 15px 20px">

An :class:`~anndata.AnnData` object ``adata`` can be sliced like a pandas
dataframe, for example, ``adata_subset = adata[:,
list_of_variable_names]``. :class:`~anndata.AnnData`'s basic structure is
similar to R's ExpressionSet [Huber15]_. If setting an `.h5ad`-formatted HDF5
backing file ``.filename``, data remains on the disk but is automatically loaded
into memory if needed. See this `blog post
<http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_
for more details.
"""


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sparse.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


class BoundRecArr(np.recarray):
    """A `np.recarray` to which fields can be added using `.['key']`.

    To enable this, it is bound to a instance of AnnData.
    """
    _attr_choices = ['obsm', 'varm']

    def __new__(cls, input_array, parent, attr):
        """
        Parameters
        ----------
        source : `np.ndarray`
            A (structured) numpy array.
        parent : object
            Any object to which the BoundRecArr shall be bound to.
        attr : string
            The name of the attribute as which it appears in parent.
        """
        arr = np.asarray(input_array).view(cls)
        arr._parent = parent
        arr._attr = attr
        return arr

    def __array_finalize__(self, obj):
        if obj is None: return
        self._parent = getattr(obj, '_parent', None)
        self._attr = getattr(obj, '_attr', None)

    def copy(self, order='C'):
        new = super(BoundRecArr, self).copy()
        new._parent = self._parent
        return new

    def flipped(self):
        new_attr = (self._attr_choices[1] if self._attr == self._attr_choices[0]
                    else self._attr_choices[0])
        return BoundRecArr(self, self._parent, new_attr)

    def keys(self):
        return self.dtype.names

    def __setitem__(self, key, arr):
        if arr.ndim == 1:
            raise ValueError('Use adata.obs or adata.var for 1-dimensional arrays.')
        if self.shape[0] != arr.shape[0]:
            raise ValueError('Can only assign an array of same length ({}), '
                             'not of length {}.'
                             .format(self.shape[0], arr.shape[0]))
        # the following always allocates a new array
        # even if the key already exists and dimensions match
        # TODO: one could check for this case
        # dtype
        merged_dtype = []
        found_key = False
        for descr in self.dtype.descr:
            if descr[0] == key:
                merged_dtype.append((key, arr.dtype, arr.shape[1]))
                found_key = True
            else:
                merged_dtype.append(descr)
        if not found_key:
            merged_dtype.append((key, arr.dtype, arr.shape[1]))
        # create new array
        new = np.empty(len(self), dtype=merged_dtype)
        # fill the array
        for name in new.dtype.names:
            if name == key:
                new[name] = arr
            else:
                new[name] = self[name]
        # make it a BoundRecArr
        # TODO: why can we not do this step before filling the array?
        new = BoundRecArr(new, self._parent, self._attr)
        setattr(self._parent, self._attr, new)

    def to_df(self):
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=RangeIndex(0, self.shape[0], name=None))
        for key in self.keys():
            value = self[key]
            for icolumn, column in enumerate(value.T):
                df['{}{}'.format(key, icolumn+1)] = column
        return df


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


def df_to_records_fixed_width(df):
    uns = {}  # unstructured dictionary for storing categories
    names = ['index']
    if is_string_dtype(df.index):
        max_len_index = df.index.map(len).max()
        index = df.index.values.astype('S{}'.format(max_len_index))
    else:
        index = df.index.values
    arrays = [index]
    for k in df.columns:
        names.append(k)
        if is_string_dtype(df[k]):
            c = df[k].astype('category')
            # transform to category only if there are less categories than
            # observations
            if len(c.cat.categories) < len(c):
                uns[k + '_categories'] = c.cat.categories.values
                arrays.append(c.cat.codes)
            else:
                max_len_index = df[k].map(len).max()
                arrays.append(df[k].values.astype('S{}'.format(max_len_index)))
        elif is_categorical(df[k]):
            uns[k + '_categories'] = df[k].cat.categories
            arrays.append(df[k].cat.codes)
        else:
            arrays.append(df[k].values)
    formats = [v.dtype for v in arrays]
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}), uns


def _fix_shapes(X, single_col=False):
    """Fix shapes of array or sparse matrix.

    Flatten 2d array with a single row or column to emulate numpys
    behavior of slicing.
    """
    if X.dtype.names is None and len(X.shape) not in {0, 1, 2}:
        raise ValueError('X needs to be 2-dimensional, not '
                         '{}-dimensional.'.format(len(X.shape)))
    if len(X.shape) == 2:
        n_obs, n_vars = X.shape
        if n_obs == 1 and n_vars == 1:
            X = X[0, 0]
        elif n_obs == 1 or n_vars == 1:
            if issparse(X): X = X.toarray()
            X = X.flatten()
    elif len(X.shape) == 1:
        n_vars = X.shape[0]
        n_obs = 1
    elif len(X.shape) == 1 and single_col:
        n_obs = X.shape[0]
        n_vars = 1
    else:
        n_vars = 1
        n_obs = 1
    return X, n_obs, n_vars


def _normalize_index(index, names):
    def name_idx(i):
        if isinstance(i, str):
            # `where` returns an 1-tuple (1D array) of found indices
            i = np.where(names == i)[0]
            if len(i) == 0:  # returns array of length 0 if nothing is found
                raise IndexError(
                    'Name "{}" is not valid observation/variable name.'
                    .format(index))
            i = i[0]
        return i

    if isinstance(index, slice):
        start = name_idx(index.start)
        stop = name_idx(index.stop)
        # string slices can only be inclusive, so +1 in that case
        if isinstance(index.stop, str):
            stop = None if stop is None else stop + 1
        step = index.step
        return slice(start, stop, step)
    elif isinstance(index, (int, str)):
        return name_idx(index)
    elif isinstance(index, (Sequence, np.ndarray)):
        return np.fromiter(map(name_idx, index), 'int64')
    else:
        raise IndexError('Unknown index {!r} of type {}'
                         .format(index, type(index)))


def _gen_dataframe(anno, length, index_names):
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
            _anno = pd.DataFrame(anno)
    return _anno


class AnnDataFileManager:
    """Backing file manager for AnnData.
    """

    def __init__(self, adata, filename=None, filemode=None):
        self._adata = adata
        self._filename = filename
        self._filemode = filemode
        if filename:
            self.open()

    def __repr__(self):
        if self._filename is None:
            return 'Backing file manager: no file is set.'
        else:
            return 'Backing file manager of file {}.'.format(self._filename)

    def __getitem__(self, key):
        return self._file[key]

    def __setitem__(self, key, value):
        self._file[key] = value

    def __delitem__(self, key):
        del self._file[key]

    @property
    def filename(self):
        return self._filename

    def open(self, filename=None, filemode=None):
        if filename is not None:
            self._filename = filename
        if filemode is not None:
            self._filemode = filemode
        if self._filename is None:
            raise ValueError(
                'Cannot open backing file if backing not initialized.')
        self._file = h5py.File(self._filename, self._filemode)

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        self._file.close()

    def _to_memory_mode(self):
        """Close the backing file, forget filename, *do* change to memory mode."""
        self._adata.__X = self._adata._read('X')[()]
        self._file.close()
        self._file = None
        self._filename = None

    @property
    def isopen(self):
        """State of backing file."""
        if self._file is None:
            return False
        # try accessing the id attribute to see if the file is open
        return bool(self._file.id)


def _init_actual_AnnData(adata_view):
    if adata_view.isbacked:
        raise ValueError(
            'You cannot modify elements of an AnnData view, '
            'but need a copy of the subset.\n\n'
            'Call `adata_subset = adata[index].copy(filename=...)`.')
    adata_view._init_as_actual(adata_view.copy())


class ArrayView(np.ndarray):

    def __new__(cls, input_array, view_args=None):
        arr = np.asarray(input_array).view(cls)
        arr._view_args = view_args
        return arr

    def __array_finalize__(self, obj):
        if obj is None: return
        self._view_args = getattr(obj, '_view_args', None)

    def __setitem__(self, idx, value):
        if self._view_args is None:
            super(type(self), self).__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            _init_actual_AnnData(adata_view)
            getattr(adata_view, attr_name)[idx] = value

    def keys(self):
        # it's a structured array
        return self.dtype.names

    def copy(self, order='C'):
        # we want a conventional array
        return np.array(self)


# the following two definitions are exactly equivalent

class SparseCSRView(sparse.csr_matrix):

    def __init__(self, *args, view_args=None, **kwargs):
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __setitem__(self, idx, value):
        if self._view_args is None:
            super(type(self), self).__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            _init_actual_AnnData(adata_view)
            getattr(adata_view, attr_name)[idx] = value


class SparseCSCView(sparse.csc_matrix):

    def __init__(self, *args, view_args=None, **kwargs):
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __setitem__(self, idx, value):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            _init_actual_AnnData(adata_view)
            getattr(adata_view, attr_name)[idx] = value


class DictView(dict):

    def __init__(self, *args, view_args=None, **kwargs):
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __setitem__(self, idx, value):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            _init_actual_AnnData(adata_view)
            getattr(adata_view, attr_name)[idx] = value


class DataFrameView(pd.DataFrame):

    _metadata = ['_view_args']

    def __init__(self, *args, view_args=None, **kwargs):
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    def __setitem__(self, idx, value):
        if self._view_args is None:
            super().__setitem__(idx, value)
        else:
            adata_view, attr_name = self._view_args
            _init_actual_AnnData(adata_view)
            getattr(adata_view, attr_name)[idx] = value


class Raw(IndexMixin):

    def __init__(self, adata=None, X=None, var=None, varm=None):
        self._adata = adata
        if X is not None:
            if len(X.shape) == 2:
                n_obs, n_vars = X.shape
                if n_obs == 1 and n_vars == 1:
                    X = X[0, 0]
                elif n_obs == 1 or n_vars == 1:
                    if issparse(X): X = X.toarray()
                    X = X.flatten()
            self._X = X
            self._var = var
            self._varm = varm
        else:
            self._X = adata.X.copy()
            self._var = adata.var.copy()
            self._varm = adata.varm.copy()

    @property
    def X(self):
        if self._adata.isbacked:
            if not self._adata.file.isopen: self._adata.file.open()
            X = self._adata.file['raw.X']
            if self._adata.isview: return X[self._adata._oidx, self._adata._vidx]
            else: return X
        else: return self._X

    @property
    def var(self):
        return self._var

    @property
    def n_vars(self):
        return self.X.shape[1]

    @property
    def varm(self):
        return self._varm

    @property
    def var_names(self):
        return self.var.index

    def __getitem__(self, index):
        oidx, vidx = self._normalize_indices(index)
        if self._adata is not None or not self._adata.isbacked: X = self._X[oidx, vidx]
        else: X = self._adata.file['raw.X'][oidx, vidx]
        if isinstance(vidx, (int, np.int64)): vidx = slice(vidx, vidx+1, 1)
        var = self._var.iloc[vidx]
        if self._varm is not None:
            varm = self._varm[vidx]
        else:
            varm = None
        return Raw(self._adata, X=X, var=var, varm=varm)

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
        obs, var = super(Raw, self)._unpack_index(packed_index)
        obs = _normalize_index(obs, self._adata.obs_names)
        var = _normalize_index(var, self.var_names)
        return obs, var


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


class AnnData(IndexMixin):
    __doc__ = dedent("""\
    An annotated data matrix.

    {main_narrative}

    Parameters
    ----------
    X : `np.ndarray`, `sparse.spmatrix`
        A #observations × #variables data matrix. A view of the data is used if the
        data type matches, otherwise, a copy is made.
    obs : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional observation annotation of length #observations.
    var : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional variable annotation of length #variables.
    uns : `dict` or `None`, optional (default: `None`)
        Unstructured annotation for the whole dataset.
    obsm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional observation annotation of length #observations.
    varm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional observation annotation of length #observations.
    dtype : simple `np.dtype`, optional (default: 'float32')
        Data type used for storage.

    See Also
    --------
    read_h5ad
    read_csv
    read_excel
    read_hdf
    read_loom
    read_mtx
    read_text
    read_umi_tools

    Notes
    -----
    Multi-dimensional annotations are stored in ``.obsm`` and ``.varm``.

    :class:`~anndata.AnnData` stores observations (samples) of variables
    (features) in the rows of a matrix. This is the convention of the modern
    classics of stats [Hastie09]_ and Machine Learning [Murphy12]_, the convention of
    dataframes both in R and Python and the established stats and machine
    learning packages in Python (`statsmodels
    <http://www.statsmodels.org/stable/index.html>`_, `scikit-learn
    <http://scikit-learn.org/>`_). It is the opposite of the convention for
    storing genomic data.

    A data matrix is flattened if either #observations (`n_obs`) or #variables
    (`n_vars`) is 1, so that Numpy's slicing behavior is reproduced::

        adata = AnnData(np.ones((2, 2)))
        adata[:, 0].X == adata.X[:, 0]

    Attributes
    ----------
    X
    filename
    isbacked
    isview
    n_obs
    n_vars
    shape
    obs
    obsm
    obs_names
    raw
    var
    varm
    var_names

    Methods
    -------
    concatenate
    copy
    transpose
    obs_names_make_unique
    var_names_make_unique
    write
    write_csvs
    write_loom
    """).format(main_narrative=_MAIN_NARRATIVE)

    _BACKED_ATTRS = ['X', 'raw.X']

    # backwards compat
    _H5_ALIASES = {
        'X': {'X', '_X', 'data', '_data'},
        'obs': {'obs', '_obs', 'smp', '_smp'},
        'var': {'var', '_var'},
        'uns': {'uns'},
        'obsm': {'obsm', '_obsm', 'smpm', '_smpm'},
        'varm': {'varm', '_varm'},
    }

    _H5_ALIASES_NAMES = {
        'obs': {'obs_names', 'smp_names', 'row_names', 'index'},
        'var': {'var_names', 'col_names', 'index'},
    }

    def __init__(
            self, X=None, obs=None, var=None, uns=None,
            obsm=None, varm=None, raw=None,
            dtype='float32', single_col=False,
            filename=None, filemode=None,
            asview=False, oidx=None, vidx=None):
        if asview:
            if not isinstance(X, AnnData):
                raise ValueError('`X` has to be an AnnData object.')
            self._init_as_view(X, oidx, vidx)
        else:
            self._init_as_actual(
                X=X, obs=obs, var=var, uns=uns,
                obsm=obsm, varm=varm, raw=raw,
                dtype=dtype, single_col=single_col,
                filename=filename, filemode=filemode)

    def _init_as_view(self, adata_ref, oidx, vidx):
        self._isview = True
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx
        # the file is the same as of the reference object
        self.file = adata_ref.file
        # views on attributes of adata_ref
        oidx_normalized, vidx_normalized = oidx, vidx
        if isinstance(oidx, (int, np.int64)): oidx_normalized = slice(oidx, oidx+1, 1)
        if isinstance(vidx, (int, np.int64)): vidx_normalized = slice(vidx, vidx+1, 1)
        obs_sub = adata_ref.obs.iloc[oidx_normalized]
        var_sub = adata_ref.var.iloc[vidx_normalized]
        self._obsm = ArrayView(adata_ref.obsm[oidx_normalized], view_args=(self, 'obsm'))
        self._varm = ArrayView(adata_ref.varm[vidx_normalized], view_args=(self, 'varm'))
        # hackish solution here, no copy should be necessary
        uns_new = self._adata_ref._uns.copy()
        # fix _n_obs, _n_vars
        if isinstance(oidx, slice):
            self._n_obs = len(obs_sub.index)
        elif isinstance(oidx, (int, np.int64)):
            self._n_obs = 1
        elif isinstance(oidx, Sized):
            self._n_obs = len(oidx)
        else:
            raise KeyError('Unknown Index type')
        if isinstance(vidx, slice):
            self._n_vars = len(var_sub.index)
        elif isinstance(vidx, (int, np.int64)):
            self._n_vars = 1
        elif isinstance(vidx, Sized):
            self._n_vars = len(vidx)
        else:
            raise KeyError('Unknown Index type')
        # need to do the slicing after setting self._n_obs, self._n_vars
        self._slice_uns_sparse_matrices_inplace(uns_new, self._oidx)
        # fix categories
        self._remove_unused_categories(adata_ref.obs, obs_sub, uns_new)
        self._remove_unused_categories(adata_ref.var, var_sub, uns_new)
        # set attributes
        self._obs = DataFrameView(obs_sub, view_args=(self, 'obs'))
        self._var = DataFrameView(var_sub, view_args=(self, 'var'))
        self._uns = DictView(uns_new, view_args=(self, 'uns'))
        # set data
        if self.isbacked: self._X = None
        else: self._init_X_as_view()
        # set raw, easy, as it's immutable anyways...
        if adata_ref._raw is not None:
            # slicing along variables axis is ignored
            self._raw = adata_ref.raw[oidx]
        else:
            self._raw = None

    def _init_X_as_view(self):
        X = self._adata_ref.X[self._oidx, self._vidx]
        if len(X.shape) == 2:
            n_obs, n_vars = X.shape
            if n_obs == 1 and n_vars == 1:
                X = X[0, 0]
            elif n_obs == 1 or n_vars == 1:
                if issparse(X): X = X.toarray()
                X = X.flatten()
        if isinstance(X, sparse.csr_matrix):
            self._X = SparseCSRView(X, view_args=(self, 'X'))
        elif isinstance(X, sparse.csc_matrix):
            self._X = SparseCSCView(X, view_args=(self, 'X'))
        elif issparse(X):
            raise ValueError('View on non-csr/csc sparse matrices not implemented.')
        else:
            self._X = ArrayView(X, view_args=(self, 'X'))

    def _init_as_actual(
            self, X=None, obs=None, var=None, uns=None,
            obsm=None, varm=None, raw=None,
            dtype='float32', single_col=False,
            filename=None, filemode=None):
        from .readwrite.read import _read_h5ad

        # view attributes
        self._isview = False
        self._adata_ref = None
        self._oidx = None
        self._vidx = None

        # ----------------------------------------------------------------------
        # various ways of initializing the data
        # ----------------------------------------------------------------------

        # init from file
        if filename is None:
            self.file = AnnDataFileManager(self, None)
        else:
            if any((X, obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If initializing from `filename`, '
                    'no further arguments may be passed.')
            self.file = AnnDataFileManager(self, filename, filemode)
            # will read from backing file
            # what is returned is, at this stage, a dict
            # that needs to be processed
            X = _read_h5ad(self, mode=filemode)

        # init from dictionary
        if isinstance(X, Mapping):
            if any((obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If `X` is a dict no further arguments must be provided.')
            X, obs, var, uns, obsm, varm, raw = self._from_dict(X)

        # init from AnnData
        if isinstance(X, AnnData):
            if any((obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If `X` is a dict no further arguments must be provided.')
            X, obs, var, uns, obsm, varm, raw = X.X, X.obs, X.var, X.uns, X.obsm, X.varm, X.raw

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
            # if type doesn't match, a copy is made, otherwise, use a view
            if issparse(X) or isinstance(X, ma.MaskedArray):
                # TODO: maybe use view on data attribute of sparse matrix
                #       as in readwrite.read_10x_h5
                if X.dtype != np.dtype(dtype): X = X.astype(dtype)
            else:  # is np.ndarray
                X = X.astype(dtype, copy=False)
            # data matrix and shape
            self._X, self._n_obs, self._n_vars = _fix_shapes(X)
        else:
            self._X = None
            self._n_obs = len(obs)
            self._n_vars = len(var)

        # annotations
        self._obs = _gen_dataframe(obs, self._n_obs,
                                   ['obs_names', 'row_names', 'smp_names'])
        self._var = _gen_dataframe(var, self._n_vars, ['var_names', 'col_names'])

        # unstructured annotations
        self._uns = uns or {}

        # multi-dimensional array annotations
        if obsm is None:
            try:
                obsm = np.empty(self._n_obs, dtype=[])
            except TypeError:
                raise TypeError(
                    'TypeError: Empty data-type'
                    '--> try installing a more recent numpy version: '
                    '    pip install numpy --upgrade')
        if varm is None: varm = np.empty(self._n_vars, dtype=[])
        self._obsm = BoundRecArr(obsm, self, 'obsm')
        self._varm = BoundRecArr(varm, self, 'varm')

        self._check_dimensions()
        self._check_uniqueness()

        # raw
        if raw is None:
            self._raw = None
        else:
            if isinstance(raw, Raw):
                self._raw = raw
            else:
                # is dictionary from reading the file
                self._raw = Raw(
                    self, X=raw['X'], var=raw['var'], varm=raw['varm'])

        # clean up old formats
        self._clean_up_old_format(uns)

    def __sizeof__(self):
        size = 0
        for attr in ['_X', '_obs', '_var', '_uns', '_obsm', '_varm']:
            s = getattr(self, attr).__sizeof__()
            size += s
        return size

    def _gen_repr(self, n_obs, n_vars):
        if self.isbacked:
            backed_at = 'backed at \'{}\''.format(self.filename)
        else:
            backed_at = ''
        descr = (
            'AnnData object with n_obs × n_vars = {} × {} {}'
            .format(n_obs, n_vars, backed_at))
        for attr in ['obs_keys', 'var_keys', 'uns_keys', 'obsm_keys', 'varm_keys']:
            keys = getattr(self, attr)()
            if len(keys) > 0:
                descr += '\n    {} = {}'.format(attr, keys)
        return descr

    def __repr__(self):
        if self.isview:
            return 'View of ' + self._gen_repr(self.n_obs, self.n_vars)
        else:
            return self._gen_repr(self.n_obs, self.n_vars)

    @property
    def shape(self):
        """Shape of data matrix: (n_obs, n_vars)."""
        return self.n_obs, self.n_vars

    @property
    def X(self) -> Union[np.ndarray, sparse.spmatrix]:
        """Data matrix of shape `n_obs` × `n_vars` (`np.ndarray`, `sp.sparse.spmatrix`)."""
        if self.isbacked:
            if not self.file.isopen: self.file.open()
            X = self.file['X']
            if self.isview: return X[self._oidx, self._vidx]
            else: return X
        else: return self._X

    @X.setter
    def X(self, value):
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
                    if (issparse(value) and not issparse(self._adata_ref._X)
                            or not issparse(value) and issparse(self._adata_ref._X)):
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
    def raw(self):
        """Store raw version of `.X` and `.var` as `.raw.X` and `.raw.X`.

        `.raw` can be initialized by setting ``adata.raw = adata1``
        """
        return self._raw

    @raw.setter
    def raw(self, value):
        if not isinstance(value, AnnData):
            raise ValueError('Can only init raw attribute with an AnnData object.')
        if self.isview:
            self._init_as_actual(self.copy())
        self._raw = Raw(value)

    @property
    def n_obs(self):
        """Number of observations."""
        return self._n_obs

    @property
    def n_vars(self):
        """Number of variables/features."""
        return self._n_vars

    @property
    def obs(self):
        """One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    @obs.setter
    def obs(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        if self.isview: self._adata_ref._obs.iloc[self._oidx] = value
        else: self._obs = value

    @property
    def var(self):
        """One-dimensional annotation of variables/ features (`pd.DataFrame`)."""
        return self._var

    @var.setter
    def var(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_vars:
            raise ValueError('Length does not match.')
        if self.isview: self._adata_ref._var.iloc[self._vidx] = value
        else: self._var = value

    @property
    def uns(self):
        """Unstructured annotation (ordered dictionary)."""
        return self._uns

    @uns.setter
    def uns(self, value):
        if self.isview:
            # here, we directly generate the copy
            adata = self._adata_ref._getitem_copy((self._oidx, self._vidx))
            self._init_as_actual(adata)
        self._uns = value

    @property
    def obsm(self):
        """Multi-dimensional annotation of observations (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_obs`. Is sliced with `data` and `obs` but behaves otherwise like a
        `dict`.
        """
        return self._obsm

    @obsm.setter
    def obsm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        if self.isview:
            self._adata_ref._obsm[self._oidx] = value
        else:
            value = BoundRecArr(value, self, 'obsm')
            self._obsm = value

    @property
    def varm(self):
        """Multi-dimensional annotation of variables/ features (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_vars`. Is sliced with `data` and `var` but behaves otherwise like a
        `dict`.
        """
        return self._varm

    @varm.setter
    def varm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_vars:
            raise ValueError('Length does not match.')
        if self.isview:
            self._adata_ref._varm[self._vidx] = value
        else:
            value = BoundRecArr(value, self, 'varm')
            self._varm = value

    @property
    def obs_names(self):
        """Names of observations (alias for `.obs.index`)."""
        return self.obs.index

    @obs_names.setter
    def obs_names(self, names):
        self._obs.index = names
        if not self._obs.index.is_unique:
            utils.warn_names_duplicates('obs', self._obs)

    @property
    def var_names(self):
        """Names of variables (alias for `.var.index`)."""
        return self._var.index

    @var_names.setter
    def var_names(self, names):
        self._var.index = names
        if not self._var.index.is_unique:
            utils.warn_names_duplicates('var', self._var)

    def obs_keys(self):
        """List keys of observation annotation `.obs`."""
        return self._obs.keys().tolist()

    def var_keys(self):
        """List keys of variable annotation `var`."""
        return self._var.keys().tolist()

    def obsm_keys(self):
        """List keys of observation annotation `obsm`."""
        return list(self._obsm.keys())

    def varm_keys(self):
        """List keys of variable annotation `varm`."""
        return list(self._varm.keys())

    def uns_keys(self):
        """List keys of unstructured annotation."""
        return sorted(list(self._uns.keys()))

    @property
    def isbacked(self):
        """`True` if object is backed on disk, `False` otherwise."""
        return self.filename is not None

    @property
    def isview(self):
        """`True` if object is view of another AnnData object, `False` otherwise."""
        return self._isview

    @property
    def filename(self):
        """Change to backing mode by setting the filename of a `.h5ad` file.

        - Setting the filename writes the stored data to disk.
        - Setting the filename when the filename was previously another name
          moves the backing file from the previous file to the new file. If you
          want to copy the previous file, use copy(filename='new_filename').
        """
        return self.file.filename

    @filename.setter
    def filename(self, filename):
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
                    os.rename(self.filename, filename)
                else:
                    # do nothing
                    return
            else:
                # change from memory to backing-mode
                # write the content of self to disk
                self.write(filename)
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

    def _normalize_indices(self, packed_index):
        # deal with slicing with pd.Series
        if isinstance(packed_index, pd.Series):
            packed_index = packed_index.values
        if isinstance(packed_index, tuple) and len(packed_index) == 2:
            if isinstance(packed_index[1], pd.Series):
                packed_index = packed_index[0], packed_index[1].values
            if isinstance(packed_index[0], pd.Series):
                packed_index = packed_index[0].values, packed_index[1]
        obs, var = super(AnnData, self)._unpack_index(packed_index)
        obs = _normalize_index(obs, self.obs_names)
        var = _normalize_index(var, self.var_names)
        return obs, var

    # TODO: this is not quite complete...
    def __delitem__(self, index):
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

    def __getitem__(self, index):
        """Returns a sliced view of the object."""
        return self._getitem_view(index)
        # return self._getitem_copy(index)

    def _getitem_view(self, index):
        oidx, vidx = self._normalize_indices(index)
        return AnnData(self, oidx=oidx, vidx=vidx, asview=True)

    def _getitem_copy(self, index):
        oidx, vidx = self._normalize_indices(index)
        if isinstance(oidx, (int, np.int64)): oidx = slice(oidx, oidx+1, 1)
        if isinstance(vidx, (int, np.int64)): vidx = slice(vidx, vidx+1, 1)
        if not self.isbacked: X = self._X[oidx, vidx]
        else: X = self.file['X'][oidx, vidx]
        obs_new = self._obs.iloc[oidx]
        self._remove_unused_categories(self._obs, obs_new, self._uns)
        obsm_new = self._obsm[oidx]
        var_new = self._var.iloc[vidx]
        self._remove_unused_categories(self._var, var_new, self._uns)
        varm_new = self._varm[vidx]
        assert obs_new.shape[0] == X.shape[0], (oidx, obs_new)
        assert var_new.shape[0] == X.shape[1], (vidx, var_new)
        uns_new = self._uns.copy()
        self._slice_uns_sparse_matrices_inplace(uns_new, oidx)
        raw_new = None if self.raw is None else self.raw[oidx]
        return AnnData(X, obs_new, var_new, uns_new, obsm_new, varm_new, raw=raw_new)

    def _remove_unused_categories(self, df_full, df_sub, uns):
        from pandas.api.types import is_categorical
        for k in df_full:
            if is_categorical(df_full[k]):
                all_categories = df_full[k].cat.categories
                df_sub[k].cat.remove_unused_categories(inplace=True)
                # also correct the colors...
                if k + '_colors' in uns:
                    uns[k + '_colors'] = np.array(uns[k + '_colors'])[
                        np.where(np.in1d(
                            all_categories, df_sub[k].cat.categories))[0]]

    def _sanitize(self):
        """Transform string arrays to categorical data types, if they store less
        categories than the total number of samples.
        """
        for ann in ['obs', 'var']:
            for key in getattr(self, ann).columns:
                df = getattr(self, ann)
                if is_string_dtype(df[key]):
                    c = pd.Categorical(
                        df[key], categories=natsorted(np.unique(df[key])))
                    if len(c.categories) < len(c):
                        df[key] = c
                        df[key].cat.categories = df[key].cat.categories.astype('U')
                        logg.info(
                            '... storing {} as categorical type'.format(key))
                        logg.info(
                            '    access categories as adata.{}[\'{}\'].cat.categories'
                            .format(ann, key))

    def _slice_uns_sparse_matrices_inplace(self, uns, oidx):
        # slice sparse spatrices of n_obs × n_obs in self.uns
        if not (isinstance(oidx, slice) and
                oidx.start is None and oidx.step is None and oidx.stop is None):
            for k, v in uns.items():
                if isinstance(v, sparse.spmatrix) and v.shape == (
                        self.n_obs, self.n_obs):
                    uns[k] = v.tocsc()[:, self.n_obs].tocsr()[oidx, :]

    def _inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        if not self.isbacked:
            self._X = self._X[:, index]
            self._n_vars = self._X.shape[1]
        else:
            X = self.file['X']
            X = X[:, index]
            self._n_vars = X.shape[1]
            self._set_backed('X', X)
        self._var = self._var.iloc[index]
        # TODO: the following should not be necessary!
        self._varm = BoundRecArr(self._varm[index], self, 'varm')
        return None

    def _inplace_subset_obs(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        if not self.isbacked:
            self._X = self._X[index, :]
            self._n_obs = self._X.shape[0]
        else:
            X = self.file['X']
            X = X[index, :]
            self._n_obs = X.shape[0]
            self._set_backed('X', X)
        self._slice_uns_sparse_matrices_inplace(self._uns, index)
        self._obs = self._obs.iloc[index]
        # TODO: the following should not be necessary!
        self._obsm = BoundRecArr(self._obsm[index], self, 'obsm')
        return None

    def _get_obs_array(self, k):
        """Get an array along the observation dimension by first looking up
        obs_keys and then var_names."""
        x = (self._obs[k] if k in self.obs_keys()
             else self[:, k].data if k in set(self.var_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in obs_keys or var_names.'
                             .format(k))
        return x

    def _get_var_array(self, k):
        """Get an array along the variables dimension by first looking up
        var_keys and then obs_names."""
        x = (self._var[k] if k in self.var_keys()
             else self[k] if k in set(self.obs_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in var_keys or obs_names.'
                             .format(k))
        return x

    def __setitem__(self, index, val):
        if self.isview:
            raise ValueError('Object is view and cannot be accessed with `[]`.')
        obs, var = self._normalize_indices(index)
        if not self.isbacked:
            self._X[obs, var] = val
        else:
            X = self.file['X']
            X[obs, var] = val
            self._set_backed('X', X)

    def __len__(self):
        return self.shape[0]

    def transpose(self):
        """Transpose whole object.

        Data matrix is transposed, observations and variables are interchanged.
        """
        if not self.isbacked: X = self._X
        else: X = self.file['X']
        if sparse.isspmatrix_csr(X):
            return AnnData(X.T.tocsr(), self._var, self._obs, self._uns,
                           self._varm.flipped(), self._obsm.flipped(),
                           filename=self.filename)
        return AnnData(X.T, self._var, self._obs, self._uns,
                       self._varm.flipped(), self._obsm.flipped(),
                       filename=self.filename)

    T = property(transpose)

    def copy(self, filename=None):
        """Full copy, optionally on disk."""
        if not self.isbacked:
            return AnnData(self._X.copy(),
                           self._obs.copy(),
                           self._var.copy(), self._uns.copy(),
                           self._obsm.copy(), self._varm.copy(),
                           raw=None if self._raw is None else self._raw.copy())
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

    def concatenate(self, *adatas, join='inner', batch_key='batch', batch_categories=None, index_unique=None):
        """Concatenate along the observations axis.

        The `.uns` and `.varm` attributes of the passed `adatas` are ignored.

        If you use `join='outer'`, then note that this fills 0s for data that is
        non-present. Use this with care.

        Parameters
        ----------
        adatas : :class:`~anndata.AnnData`
            AnnData matrices to concatenate with.
        join: `str` (default: 'inner')
            Use intersection (``'inner'``) or union (``'outer'``) of variables.
        batch_key : `str` (default: 'batch')
            Add the batch annotation to `.obs` using this key.
        batch_categories : list, optional (default: `range(len(adatas)+1)`)
            Use these as categories for the batch annotation.
        index_unique : `str` or `None`, optional (default: `None`)
            Make the index unique by joining the existing index names with the
            batch category. Provide `None` to keep existing indices.

        Returns
        -------
        adata : :class:`~anndata.AnnData`
            The concatenated AnnData, where `adata.obs[batch_key]` stores a
            categorical variable labeling the batch.

        Examples
        --------
        >>> adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'anno1': ['c1', 'c2']},
        >>>                  {'var_names': ['a', 'b', 'c']})
        >>> adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'anno1': ['c3', 'c4']},
        >>>                  {'var_names': ['b', 'c', 'd']})
        >>> adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'anno2': ['d3', 'd4']},
        >>>                  {'var_names': ['b', 'c', 'd']})
        >>>
        >>> adata = adata1.concatenate(adata2, adata3, index_unique='-')
        >>> adata.X
        [[ 2.  3.]
         [ 5.  6.]
         [ 1.  2.]
         [ 4.  5.]
         [ 1.  2.]
         [ 4.  5.]]
        >>> adata.obs_names
        Index(['0-0', '1-0', '0-1', '1-1', '0-2', '1-2'], dtype='object')
        >>> adata.obs
            anno1 anno2 batch
        0-0    c1   NaN     0
        1-0    c2   NaN     0
        0-1    c3   NaN     1
        1-1    c4   NaN     1
        0-2   NaN    d3     2
        1-2   NaN    d4     2
        """
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
                    logg.info(
                        'Making variable names unique for controlled concatenation.')
                    printed_info = True

        # define variable names of joint AnnData
        mergers = dict(inner=set.intersection, outer=set.union)
        var_names_reduce = reduce(mergers[join], (set(ad.var_names) for ad in all_adatas))
        # restore order of initial var_names, append non-sortable names at the end
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
        mat_cls = sparse.csc_matrix if any_sparse else np.ndarray
        X = mat_cls(out_shape, dtype=self.X.dtype)
        var = pd.DataFrame(index=var_names)

        obs_i = 0  # start of next adata’s observations in X
        out_obss = []
        for i, ad in enumerate(all_adatas):
            vars_intersect = [v for v in var_names if v in ad.var_names]

            # X
            X[obs_i:obs_i+ad.n_obs,
              var_names.isin(vars_intersect)] = ad[:, vars_intersect].X
            obs_i += ad.n_obs

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
            # potential add additional columns
            var.loc[vars_intersect, ad.var.columns] = ad.var.loc[vars_intersect, :]

        obs = pd.concat(out_obss)
        uns = all_adatas[0].uns
        obsm = np.concatenate([ad.obsm for ad in all_adatas])
        varm = self.varm  # TODO

        new_adata = AnnData(X, obs, var, uns, obsm, None, filename=self.filename)
        if not obs.index.is_unique:
            logg.info(
                'Or pass `index_unique!=None` to `.concatenate`.')
        return new_adata

    def var_names_make_unique(self, join='-'):
        self.var.index = utils.make_index_unique(self.var.index, join)

    var_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def obs_names_make_unique(self, join='-'):
        self.obs.index = utils.make_index_unique(self.obs.index, join)

    obs_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def _check_uniqueness(self):
        if not self.obs.index.is_unique:
            utils.warn_names_duplicates('obs', self.obs)
        if not self.var.index.is_unique:
            utils.warn_names_duplicates('var', self.var)

    def __contains__(self, key):
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
        if 'obsm' in key and len(self._obsm) != self._n_obs:
            raise ValueError('Observations annot. `obsm` must have number of '
                             'rows of `X` ({}), but has {} rows.'
                             .format(self._n_obs, len(self._obsm)))
        if 'varm' in key and len(self._varm) != self._n_vars:
            raise ValueError('Variables annot. `varm` must have number of '
                             'columns of `X` ({}), but has {} rows.'
                             .format(self._n_vars, len(self._varm)))

    def write(self, filename=None, compression='gzip', compression_opts=None):
        """Write `.h5ad`-formatted hdf5 file and close a potential backing file.

        Parameters
        ----------
        filename : `str`, optional (default: ``.filename``)
            Filename of data file. Defaults to backing file.
        compression : `None` or {'gzip', 'lzf'}, optional (default: `'gzip'`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        compression_opts : `int`, optional (default: `None`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        """
        from .readwrite.write import _write_h5ad

        if filename is None and not self.isbacked:
            raise ValueError('Provide a filename!')
        if filename is None:
            filename = self.filename
        _write_h5ad(filename, self, compression=compression, compression_opts=compression_opts)
        if self.isbacked:
            self.file.close()

    def write_csvs(self, dirname, skip_data=True, sep=','):
        """Write annotation to `.csv` files.

        It is not possible to recover the full :class:`~anndata.AnnData` from the
        output of this function. Use :func:`~anndata.write` for this.

        Parameters
        ----------
        dirname : `str`
            Name of directory to which to export.
        skip_data : `bool`, optional (default: True)
             Skip the data matrix `.X`.
        sep : `str`, optional (default: ',')
             Separator for the data.
        """
        from .readwrite.write import write_csvs
        write_csvs(dirname, self, skip_data=skip_data, sep=sep)

    def write_loom(self, filename):
        """Write `.loom`-formatted hdf5 file.

        Parameters
        ----------
        filename : `str`
            The filename.
        """
        from .readwrite.write import write_loom
        write_loom(filename, self)

    @staticmethod
    def _from_dict(ddata):
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
        valid_keys += ['raw.X', 'raw.var', 'raw.varm']
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
        from pandas.api.types import is_string_dtype
        from pandas import Index

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
                        d_true_keys[true_key][c] = Index(
                            d_true_keys[true_key][c]).astype('U').values

        # these are the category fields
        k_to_delete = []
        items = (ddata.items() if uns_is_not_key
                else ddata['uns'].items() if 'uns' in ddata else [])
        for k, v in items:
            if k.endswith('_categories'):
                k_stripped = k.replace('_categories', '')
                if k_stripped in d_true_keys['obs']:
                    d_true_keys['obs'][k_stripped] = pd.Categorical.from_codes(
                        codes=d_true_keys['obs'][k_stripped].values,
                        categories=v)
                if k_stripped in d_true_keys['var']:
                    d_true_keys['var'][k_stripped] = pd.Categorical.from_codes(
                        codes=d_true_keys['var'][k_stripped].values,
                        categories=v)
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
                    raw['var'][c] = Index(raw['var'][c]).astype('U').values
        if 'raw.varm' in ddata:
            raw['varm'] = ddata['raw.varm']
            del ddata['raw.varm']
        elif raw is not None:
            raw['varm'] = None

        # the remaining fields are the unstructured annotation
        uns = (ddata if uns_is_not_key
               else ddata['uns'] if 'uns' in ddata else {})

        return X, obs, var, uns, obsm, varm, raw

    def _to_dict_fixed_width_arrays(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        obs_rec, uns_obs = df_to_records_fixed_width(self._obs)
        var_rec, uns_var = df_to_records_fixed_width(self._var)
        d = {
            'X': self._X,
            'obs': obs_rec,
            'var': var_rec,
            'obsm': self._obsm,
            'varm': self._varm,
            # add the categories to the unstructured annotation
            'uns': {**self._uns, **uns_obs, **uns_var}}

        if self.raw is not None:
            # we ignore categorical data types here
            # they should never occur
            var_rec, uns_var = df_to_records_fixed_width(self.raw._var)
            if len(uns_var) > 0:
                warnings.warn('Categorical dtypes in `.raw.var` are cast to integer.')
            d['raw.X'] = self.raw.X
            d['raw.var'] = var_rec
            d['raw.varm'] = self.raw.varm

        return d

    # --------------------------------------------------------------------------
    # all of the following is for backwards compat
    # --------------------------------------------------------------------------

    @property
    def data(self):
        """Deprecated access to X."""
        print('DEPRECATION WARNING: use attribute `.X` instead of `.data`, '
              '`.data` will be removed in the future.')
        return self.X

    @data.setter
    def data(self, value):
        print('DEPRECATION WARNING: use attribute `.X` instead of `.data`, '
              '`.data` will be removed in the future.')
        self.X = value

    @property
    def n_smps(self):
        """Deprecated access to `.n_obs`."""
        return self.n_obs

    @property
    def smp(self):
        """Deprecated access to `.obs`."""
        return self.obs

    # for backwards compat
    @smp.setter
    def smp(self, value):
        self.obs = value

    @property
    def add(self):
        """Deprecated access to `.uns`, remains for backwards compatibility."""
        # FutureWarning and DeprecationWarning are not visible by default, use print
        print('DEPRECATION WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        return self.uns

    @add.setter
    def add(self, value):
        print('DEPRECATION WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        self.uns = value

    @property
    def smpm(self):
        """Multi-dimensional annotation of observations (mutable structured `np.ndarray`)        """
        return self.obsm

    @smpm.setter
    def smpm(self, value):
        self.obsm = value

    @property
    def smp_names(self):
        return self.obs_names

    @smp_names.setter
    def smp_names(self, names):
        self.obs_names = names

    def smp_keys(self):
        """List keys of observation annotation `obs`."""
        return self.obs_keys()

    def smpm_keys(self):
        """List keys of observation annotation `obsm`."""
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
