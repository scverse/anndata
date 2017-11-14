# Author: P. Angerer
#         Alex Wolf (http://falexwolf.de)
"""Store an annotated data matrix.
"""
from collections import Mapping, Sequence
from collections import OrderedDict
from enum import Enum

import numpy as np
from numpy import ma
import pandas as pd
from pandas.core.index import RangeIndex
from scipy import sparse as sp
from scipy.sparse.sputils import IndexMixin
from numpy.lib.recfunctions import merge_arrays

import logging as logg
import warnings


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sp.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


class BoundRecArr(np.recarray):
    """A np.recarray to which fields can be added using [].

    To enable this, it is bound to a instance of AnnData.

    Parameters
    ----------
    source : `np.ndarray`
        A (structured) numpy array.
    parent : object
        Any object to which the BoundRecArr shall be bound to.
    attr : string
        The name of the attribute as which it appears in parent.
    """
    _attr_choices = ['smpm', 'varm']

    def __new__(cls, input_array, parent, attr):
        arr = np.asarray(input_array).view(cls)
        arr._parent = parent
        arr._attr = attr
        return arr

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self._parent = getattr(obj, '_parent', None)
        self._attr = getattr(obj, '_attr', None)

    def copy(self):
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
            raise ValueError('Use adata.smp or adata.var for 1-dimensional arrays.')
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
    from pandas.api.types import is_string_dtype, is_categorical
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
            uns[k + '_categories'] = c.cat.categories.values
            arrays.append(c.cat.codes)
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
        # TODO: int immutable, copy of references to ints in X.shape
        # only valid until accidental change
        n_smps, n_vars = X.shape
        if n_smps == 1 and n_vars == 1:
            X = X[0, 0]
        elif n_smps == 1 or n_vars == 1:
            if sp.issparse(X): X = X.toarray()
            X = X.flatten()
    elif len(X.shape) == 1:
        n_vars = X.shape[0]
        n_smps = 1
    elif len(X.shape) == 1 and single_col:
        n_smps = X.shape[0]
        n_vars = 1
    else:
        n_vars = 1
        n_smps = 1
    return X, n_smps, n_vars


def _gen_dataframe(anno, length, index_name1, index_name2):
    if anno is None or len(anno) == 0:
        _anno = pd.DataFrame(index=RangeIndex(0, length, name=None))
    elif index_name1 in anno:
        _anno = pd.DataFrame(anno, index=anno[index_name1],
                                columns=[k for k in anno.keys() if k != index_name1])
    elif index_name2 in anno:
        _anno = pd.DataFrame(anno, index=anno[index_name2],
                                columns=[k for k in anno.keys() if k != index_name2])
    else:
        _anno = pd.DataFrame(anno)
    return _anno


class AnnData(IndexMixin):
    """Store an annotated data matrix.

    Stores a data matrix ``.data`` as an array or sparse matrix, annotation of
    samples ``.smp`` and variables ``.var`` as dataframes and unstructured
    annotation ``.uns`` as dictionary.

    Parameters
    ----------
    data : `np.ndarray`, `sp.spmatrix`, `np.ma.MaskedArray` or `dict`
        If not a `dict`, this is a #samples × #variables data matrix. A view of
        the data is used if the data type matches, otherwise, a copy is made. If
        a `dict`, it must contain `data` as `'data'` and optionally arrays for
        `'row_names'` (`'smp_names'`) and `'col_names'` (`'var_names'`).
    smp : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional sample annotation of length #samples.
    var : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional variable annotation of length #variables.
    uns : `dict` or `None`, optional (default: `None`)
        Unstructured annotation for the whole dataset.
    smpm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional sample annotation of length #samples.
    varm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional sample annotation of length #samples.
    dtype : simple `np.dtype`, optional (default: 'float32')
        Data type used for storage.
    single_col : `bool`, optional (default: `False`)
        Interpret one-dimensional input array as column.

    Notes
    -----
    A data matrix is flattened if either #samples (`n_smps`) or #variables
    (`n_vars`) is 1, so that Numpy's slicing behavior is reproduced::

        ad = AnnData(np.ones((2, 2)))
        ad[:, 0].X == ad.X[:, 0]

    Examples
    --------
    >>> adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
    >>>                  {'smp_names': ['s1', 's2'],
    >>>                   'anno1': ['c1', 'c2']},
    >>>                  {'var_names': ['a', 'b', 'c']})
    >>> adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
    >>>                  {'smp_names': ['s3', 's4'],
    >>>                   'anno1': ['c3', 'c4']},
    >>>                  {'var_names': ['b', 'c', 'd']})
    >>> adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
    >>>                  {'smp_names': ['s5', 's6'],
    >>>                   'anno2': ['d3', 'd4']},
    >>>                  {'var_names': ['b', 'c', 'd']})
    >>>
    >>> adata = adata1.concatenate([adata2, adata3])
    >>> adata.X
    [[ 2.  3.]
     [ 5.  6.]
     [ 1.  2.]
     [ 4.  5.]
     [ 1.  2.]
     [ 4.  5.]]
    >>> adata.smp
       anno1 anno2 batch
    s1    c1   NaN     0
    s2    c2   NaN     0
    s3    c3   NaN     1
    s4    c4   NaN     1
    s5   NaN    d3     2
    s6   NaN    d4     2
    """

    def __init__(self, data, smp=None, var=None, uns=None, smpm=None, varm=None,
                 dtype='float32', single_col=False):

        # generate from a dictionary
        if isinstance(data, Mapping):
            if any((smp, var, uns, smpm, varm)):
                raise ValueError(
                    'If `data` is a dict no further arguments must be provided.')
            data, smp, var, uns, smpm, varm = self._from_dict(data)

        # check data type of data
        for s_type in StorageType:
            if isinstance(data, s_type.value):
                break
        else:
            class_names = ', '.join(c.__name__ for c in StorageType.classes())
            raise ValueError('data needs to be of one of {}, not {}'
                             .format(class_names, type(data)))

        # if type doesn't match, a copy is made, otherwise, use a view
        if sp.issparse(data) or isinstance(data, ma.MaskedArray):
            # TODO: maybe use view on data attribute of sparse matrix
            #       as in readwrite.read_10x_h5
            if data.dtype != np.dtype(dtype): data = data.astype(dtype)
        else:  # is np.ndarray
            data = data.astype(dtype, copy=False)

        # data matrix and shape
        self._data, self._n_smps, self._n_vars = _fix_shapes(data)

        # annotations
        self._smp = _gen_dataframe(smp, self._n_smps, 'smp_names', 'row_names')
        self._var = _gen_dataframe(var, self._n_vars, 'var_names', 'col_names')

        # unstructured annotations
        self._uns = uns or {}

        # multi-dimensional array annotations
        if smpm is None: smpm = np.empty(self._n_smps, dtype=[])
        if varm is None: varm = np.empty(self._n_vars, dtype=[])
        self._smpm = BoundRecArr(smpm, self, 'smpm')
        self._varm = BoundRecArr(varm, self, 'varm')

        self._check_dimensions()

        # clean up old formats
        self._clean_up_old_format(uns)

    def __repr__(self):
        return ('AnnData object with n_smps × n_vars= {} × {}\n'
                '    smp_keys = {}\n'
                '    var_keys = {}\n'
                '    uns_keys = {}\n'
                '    smpm_keys = {}\n'
                '    varm_keys = {}'
                .format(self._n_smps, self._n_vars,
                        self.smp_keys(), self.var_keys(),
                        self.uns_keys(),
                        self.smpm_keys(), self.varm_keys()))

    @property
    def data(self):
        """Data matrix of shape `n_smps` × `n_vars` (`np.ndarray`, `sp.spmatrix`, `np.ma.MaskedArray`)."""
        return self._data

    @data.setter
    def data(self, value):
        if not self._data.shape == value.shape:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self._data.shape))
        self._data = value

    @property
    def X(self):
        """Shortcut for data matrix `data`."""
        return self._data

    @X.setter
    def X(self, value):
        if not self._data.shape == value.shape:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self._data.shape))
        self._data = value

    @property
    def n_smps(self):
        """Number of samples (rows)."""
        return self._n_smps

    @property
    def n_vars(self):
        """Number of variables/ features (columns)."""
        return self._n_vars

    @property
    def smp(self):
        """One-dimensional annotation of samples (`pd.DataFrame`)."""
        return self._smp

    @smp.setter
    def smp(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_smps:
            raise ValueError('Length does not match.')
        self._smp = value

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
        self._var = value

    @property
    def uns(self):
        """Unstructured annotation (ordered dictionary)."""
        return self._uns

    @uns.setter
    def uns(self, value):
        self._uns = value

    @property
    def add(self):
        """Deprecated, remains for backwards compatibility."""
        # FutureWarning and DeprecationWarning are not visible by default, use print
        print('WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        return self._uns

    @add.setter
    def add(self, value):
        print('WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        self._uns = value

    @property
    def smpm(self):
        """Multi-dimensional annotation of samples (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_smps`. Is sliced with `data` and `smp` but behaves otherwise like a
        `dict`.
        """
        return self._smpm

    @smpm.setter
    def smpm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_smps:
            raise ValueError('Length does not match.')
        value = BoundRecArr(value, self, 'smpm')
        self._smpm = value

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
        value = BoundRecArr(value, self, 'varm')
        self._varm = value

    @property
    def smp_names(self):
        """Index for samples (`smp.index`)."""
        return self._smp.index.values

    @smp_names.setter
    def smp_names(self, names):
        self._smp.index = names

    @property
    def var_names(self):
        """Index for variables (`var.index`)."""
        return self._var.index.values

    @var_names.setter
    def var_names(self, names):
        self._var.index = names

    def smp_keys(self):
        """List keys of sample annotation `smp`."""
        return self._smp.keys().tolist()

    def var_keys(self):
        """List keys of variable annotation `var`."""
        return self._var.keys().tolist()

    def smpm_keys(self):
        """List keys of sample annotation `smpm`."""
        return list(self._smpm.keys())

    def varm_keys(self):
        """List keys of variable annotation `varm`."""
        return list(self._varm.keys())

    def uns_keys(self):
        """List keys of unstructured annotation."""
        return sorted(list(self._uns.keys()))

    def _normalize_indices(self, packed_index):
        smp, var = super(AnnData, self)._unpack_index(packed_index)
        smp = self._normalize_index(smp, self.smp_names)
        var = self._normalize_index(var, self.var_names)
        return smp, var

    def _normalize_index(self, index, names):
        def name_idx(i):
            if isinstance(i, str):
                # `where` returns an 1-tuple (1D array) of found indices
                i = np.where(names == i)[0]
                if len(i) == 0:  # returns array of length 0 if nothing is found
                    raise IndexError('Name "{}" is not valid variable or sample index.'
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
        elif isinstance(index, (int, str)):
            start = name_idx(index)
            stop = start + 1
            step = 1
        elif isinstance(index, (Sequence, np.ndarray)):
            return np.fromiter(map(name_idx, index), 'int64')
        else:
            raise IndexError('Unknown index {!r} of type {}'
                             .format(index, type(index)))

        return slice(start, stop, step)

    def __delitem__(self, index):
        smp, var = self._normalize_indices(index)
        del self._data[smp, var]
        if var == slice(None):
            del self._smp.iloc[smp, :]
        if smp == slice(None):
            del self._var.iloc[var, :]

    def __getitem__(self, index):
        # Note: this cannot be made inplace
        # http://stackoverflow.com/questions/31916617/using-keyword-arguments-in-getitem-method-in-python
        smp, var = self._normalize_indices(index)
        data = self._data[smp, var]
        smp_new = self._smp.iloc[smp]
        smpm_new = self._smpm[smp]
        var_new = self._var.iloc[var]
        varm_new = self._varm[var]
        assert smp_new.shape[0] == data.shape[0], (smp, smp_new)
        assert var_new.shape[0] == data.shape[1], (var, var_new)
        uns_new = self._uns.copy()
        # slice sparse spatrices of n_smps × n_smps in self._uns
        if not (isinstance(smp, slice) and
                smp.start is None and smp.step is None and smp.stop is None):
            raised_warning = False
            for k, v in self._uns.items():  # TODO: make sure this really works as expected
                if isinstance(v, sp.spmatrix) and v.shape == (self._n_smps, self._n_smps):
                    uns_new[k] = v.tocsc()[:, smp].tocsr()[smp, :]
                    if not raised_warning:
                        logg.warn('Slicing adjacency matrices can be dangerous. '
                                  'Consider recomputing the data graph.')
                        raised_warning = True
        adata = AnnData(data, smp_new, var_new, uns_new, smpm_new, varm_new)
        return adata

    def _inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        self._data = self._data[:, index]
        self._n_vars = self._data.shape[1]
        self._var = self._var.iloc[index]
        # TODO: the following should not be necessary!
        self._varm = BoundRecArr(self._varm[index], self, 'varm')
        return None

    def _inplace_subset_smp(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        self._data = self._data[index, :]
        raised_warning = False
        # TODO: solve this in a better way, also for var
        for k, v in self._uns.items():
            if isinstance(v, sp.spmatrix) and v.shape == (self._n_smps, self._n_smps):
                self._uns[k] = v.tocsc()[:, index].tocsr()[index, :]
                if not raised_warning:
                    logg.warn('Slicing adjacency matrices can be dangerous. '
                              'Consider recomputing the data graph.')
                    raised_warning = True
        self._n_smps = self._data.shape[0]
        self._smp = self._smp.iloc[index]
        # TODO: the following should not be necessary!
        self._smpm = BoundRecArr(self._smpm[index], self, 'smpm')
        return None

    def _get_smp_array(self, k):
        """Get an array along the sample dimension by first looking up
        smp_keys and then var_names."""
        x = (self._smp[k] if k in self.smp_keys()
             else self[:, k].data if k in set(self.var_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in smp_keys or var_names.'
                             .format(k))
        return x

    def _get_var_array(self, k):
        """Get an array along the variables dimension by first looking up
        var_keys and then smp_names."""
        x = (self._var[k] if k in self.var_keys()
             else self[k] if k in set(self.smp_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in var_keys or smp_names.'
                             .format(k))
        return x

    def __setitem__(self, index, val):
        smp, var = self._normalize_indices(index)
        self._data[smp, var] = val

    def __len__(self):
        return self._data.shape[0]

    def transpose(self):
        """Transpose whole object.

        Data matrix is transposed, samples and variables are interchanged.
        """
        if sp.isspmatrix_csr(self._data):
            return AnnData(self._data.T.tocsr(), self._var, self._smp, self._uns,
                           self._varm.flipped(), self._smpm.flipped())
        return AnnData(self._data.T, self._var, self._smp, self._uns,
                       self._varm.flipped(), self._smpm.flipped())

    T = property(transpose)

    def copy(self):
        """Full copy with memory allocated."""
        return AnnData(self._data.copy(), self._smp.copy(), self._var.copy(), self._uns.copy(),
                       self._smpm.copy(), self._varm.copy())

    def concatenate(self, adatas):
        """Concatenate along the samples axis after intersecting the variables names.

        The `.var`, `.varm`, and `.uns` attributes of the passed adatas are ignored.

        Parameters
        ----------
        adatas : AnnData or list of AnnData
            AnnData matrices to concatenate with.

        Returns
        -------
        adata : AnnData
            The concatenated AnnData, where `adata.smp['batch']` stores a
            categorical variable labeling the batch.
        """
        if isinstance(adatas, AnnData): adatas = [adatas]
        joint_variables = self.var_names
        for adata2 in adatas:
            joint_variables = np.intersect1d(
                joint_variables, adata2.var_names, assume_unique=True)
        adatas_to_concat = []
        categories = [str(i) for i in range(len(adatas)+1)]
        for i, ad in enumerate([self] + adatas):
            ad = ad[:, joint_variables]
            ad.smp['batch'] = pd.Categorical(
                ad.n_smps*[categories[i]], categories=categories)
            adatas_to_concat.append(ad)
        Xs = [ad.X for ad in adatas_to_concat]
        if sp.issparse(self.X):
            from scipy.sparse import vstack
            X = vstack(Xs)
        else:
            X = np.concatenate(Xs)
        smp = pd.concat([ad.smp for ad in adatas_to_concat])
        smpm = np.concatenate([ad.smpm for ad in adatas_to_concat])
        var = adatas_to_concat[0].var
        varm = adatas_to_concat[0].varm
        uns = adatas_to_concat[0].uns
        return AnnData(X, smp, var, uns, smpm, varm)

    def __contains__(self, key):
        raise AttributeError('AnnData has no attribute __contains__, don\'t check `in adata`.')

    def _check_dimensions(self, key=None):
        if key is None:
            key = {'smp', 'var', 'smpm', 'varm'}
        else:
            key = {key}
        if 'smp' in key and len(self._smp) != self._n_smps:
            raise ValueError('Sample annotation `smp` needs to have the same amount of '
                             'rows as data ({}), but has {} rows'
                             .format(self._n_smps, self._smp.shape[0]))
        if 'var' in key and len(self._var) != self._n_vars:
            raise ValueError('Variable annotation `var` needs to have the same amount of '
                             'columns as data  ({}), but has {} rows'
                             .format(self._n_vars, self._var.shape[0]))
        if 'smpm' in key and len(self._smpm) != self._n_smps:
            raise ValueError('Sample annotation `smpm` needs to have the same amount of '
                             'rows as data ({}), but has {} rows'
                             .format(self._n_smps, self._smp.shape[0]))
        if 'varm' in key and len(self._varm) != self._n_vars:
            raise ValueError('Variable annotation `varm` needs to have the same amount of '
                             'columns as data ({}), but has {} rows'
                             .format(self._n_vars, self._var.shape[0]))

    def _to_dict_dataframes(self):
        d = {'data': pd.DataFrame(self._data, index=self._smp.index),
             'smp': self._smp,
             'var': self._var,
             'smpm': self._smpm.to_df(),
             'varm': self._varm.to_df()}
        d = merge_dicts(d, self._uns)
        return d

    def _to_dict_fixed_width_arrays(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        smp_rec, uns_smp = df_to_records_fixed_width(self._smp)
        var_rec, uns_var = df_to_records_fixed_width(self._var)
        d = {'_data': self._data,
             '_smp': smp_rec,
             '_var': var_rec,
             '_smpm': self._smpm,
             '_varm': self._varm}
        d = merge_dicts(d, self._uns, uns_smp, uns_var)
        return d

    def _from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.

        In particular, from a dict that has been written using
        ``AnnData._to_dict_fixed_width_arrays``.
        """
        # data matrix
        if '_data' in ddata:
            data = ddata['_data']
            del ddata['_data']
        elif 'data' in ddata:
            data = ddata['data']
            del ddata['data']
        elif '_X' in ddata:
            data = ddata['_X']
            del ddata['_X']
        elif 'X' in ddata:
            data = ddata['X']
            del ddata['X']
        # simple annotation
        if ('_smp' in ddata and isinstance(ddata['_smp'], (np.ndarray, pd.DataFrame))
            and '_var' in ddata and isinstance(ddata['_var'], (np.ndarray, pd.DataFrame))):
            smp = ddata['_smp']
            del ddata['_smp']
            var = ddata['_var']
            del ddata['_var']
        elif ('smp' in ddata and isinstance(ddata['smp'], (np.ndarray, pd.DataFrame))
              and 'var' in ddata and isinstance(ddata['var'], (np.ndarray, pd.DataFrame))):
            smp = ddata['smp']
            del ddata['smp']
            var = ddata['var']
            del ddata['var']
        else:
            smp, var = OrderedDict(), OrderedDict()
            if 'row_names' in ddata:
                smp['smp_names'] = ddata['row_names']
                del ddata['row_names']
            elif 'smp_names' in ddata:
                smp['smp_names'] = ddata['smp_names']
                del ddata['smp_names']
            if 'col_names' in ddata:
                var['var_names'] = ddata['col_names']
                del ddata['col_names']
            elif 'var_names' in ddata:
                var['var_names'] = ddata['var_names']
                del ddata['var_names']
            smp = merge_dicts(smp, ddata.get('row', {}), ddata.get('smp', {}))
            var = merge_dicts(var, ddata.get('col', {}), ddata.get('var', {}))
            for k in ['row', 'smp', 'col', 'var']:
                if k in ddata:
                    del ddata[k]

        # transform recarray to dataframe
        if isinstance(smp, np.ndarray) and isinstance(var, np.ndarray):

            if 'smp_names' in smp.dtype.names:
                smp = pd.DataFrame.from_records(smp, index='smp_names')
            elif 'row_names' in smp.dtype.names:
                smp = pd.DataFrame.from_records(smp, index='row_names')
            elif 'index' in smp.dtype.names:
                smp = pd.DataFrame.from_records(smp, index='index')
            else:
                smp = pd.DataFrame.from_records(smp)
            smp.index = smp.index.astype('U')

            if 'var_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='var_names')
            elif 'col_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='col_names')
            elif 'index' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='index')
            else:
                var = pd.DataFrame.from_records(var)
            var.index = var.index.astype('U')

            # these are the category fields
            k_to_delete = []
            for k in ddata.keys():
                if k.endswith('_categories'):
                    k_stripped = k.replace('_categories', '')
                    if k_stripped in smp:
                        smp[k_stripped] = pd.Categorical.from_codes(
                            codes=smp[k_stripped].values,
                            categories=ddata[k])
                    if k_stripped in var:
                        var[k_stripped] = pd.Categorical.from_codes(
                            codes=smp[k_stripped].values,
                            categories=ddata[k])
                    k_to_delete.append(k)

            for k in k_to_delete:
                del ddata[k]

        # multicolumn annotation
        if 'smpm' in ddata:
            smpm = ddata['smpm']
            del ddata['smpm']
        elif '_smpm' in ddata:
            smpm = ddata['_smpm']
            del ddata['_smpm']
        else:
            smpm = None
        if 'varm' in ddata:
            varm = ddata['varm']
            del ddata['varm']
        elif '_varm' in ddata:
            varm = ddata['_varm']
            del ddata['_varm']
        else:
            varm = None

        # the remaining fields are the unstructured annotation
        uns = ddata

        return data, smp, var, uns, smpm, varm

    def _clean_up_old_format(self, uns):
        # multicolumn keys
        # all of the rest is only for backwards compat
        if uns and '_smp_keys_multicol' in uns:
            _keys_multicol_smp = list(uns['_smp_keys_multicol'])
            del uns['_smp_keys_multicol']
        elif uns and 'smp_keys_multicol' in uns:
            _keys_multicol_smp = list(uns['smp_keys_multicol'])
            del uns['smp_keys_multicol']
        else:
            _keys_multicol_smp = []
        if uns and '_var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['_var_keys_multicol'])
            del uns['_var_keys_multicol']
        elif uns and 'var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['var_keys_multicol'])
            del uns['var_keys_multicol']
        else:
            _keys_multicol_var = []

        # now, for compat, fill the old multicolumn entries into smpm and varm
        # and remove them from smp and var
        for key in _keys_multicol_smp:
            self._smpm[key] = self._get_multicol_field_smp(key)
        for key in _keys_multicol_var:
            self._varm[key] = self._get_multicol_field_var(key)

    def _get_multicol_field_smp(self, key_multicol):
        return self._get_and_delete_multicol_field('smp', key_multicol)

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


def merge_dicts(*ds):
    """Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = ds[0]
    for d in ds[1:]:
        result.update(d)
    return result
