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
    source : np.ndarray
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
        # - TODO: one could check for this case
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
        # - TODO: why can we not do this step before filling the array?
        new = BoundRecArr(new, self._parent, self._attr)
        setattr(self._parent, self._attr, new)

    def to_dataframe(self):
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


class AnnData(IndexMixin):
    """Store an annotated data matrix.

    Stores a data matrix ``.X`` as an array or sparse matrix, annotation of
    samples ``.smp`` and variables ``.var`` as dataframes and unstructured
    annotation ``.uns`` as dictionary.

    Parameters
    ----------
    data : np.ndarray, np.ma.MaskedArray, sp.spmatrix, dictionary
        A #samples × #variables data matrix `X`. If a dictionary, the dictionary
        must contain `X` as ``'X'`` and optionally a field ``'row_names'``
        (``'smp_names'``) that stores sample names; similar for ``'col_names'``
        (``'var_names'``).
    smp : pd.DataFrame, dictionary, np.ndarray or None, optional (default: ``None``)
        A #samples × #sample_keys array containing sample names (`index`) and
        other sample annotation in the columns.
    var : pd.DataFrame, dictionary, np.ndarray or None, optional (default: ``None``)
        The same as `smp`, but of shape #variables × #variable_keys for
        annotation of variables.
    uns : dictionary or None, optional (default: ``None``)
        Unstructured annotation for the whole dataset.
    smpm : np.ndarray, optional (default: ``None``)
        A structured numpy array with length #samples, intended for storing
        multidimensional annotation of numeric data types.
    varm : np.ndarray, optional (default: ``None``)
        A structured numpy array with length #variables, intended for storing
        multidimensional annotation of numeric data types.
    dtype : simple `np.dtype`, optional (default: ``'float32'``)
        Convert data matrix, smpm and varm to this data type upon
        initialization.
    single_col : bool, optional (default: ``False``)
        Interpret one-dimensional input array as column.

    Notes
    -----
    A data matrix is flattened if either #samples (`n_smps`) or #variables
    (`n_vars`) is 1, so that Numpy's behavior is reproduced (``adata[:, 0].X ==
    adata.X[:, 0]``).
    """

    def __init__(self, data, smp=None, var=None, uns=None, smpm=None, varm=None,
                 dtype='float32', single_col=False):
        if isinstance(data, Mapping):
            if any((smp, var, uns, smpm, varm)):
                raise ValueError(
                    'If `data` is a dict no further arguments must be provided.')
            X, smp, var, uns, smpm, varm = self.from_dict(data)
        else:
            X = data

        # check data type of X
        for s_type in StorageType:
            if isinstance(X, s_type.value):
                self.storage_type = s_type
                break
        else:
            class_names = ', '.join(c.__name__ for c in StorageType.classes())
            raise ValueError('X needs to be of one of the following types [{}] not {}'
                             .format(class_names, type(X)))

        # type conversion: if type doesn't match, a copy is made
        if sp.issparse(X) or isinstance(X, ma.MaskedArray):
            # TODO: maybe use view on data attribute of sparse matrix
            #       as in readwrite.read_10x_h5
            if X.dtype != np.dtype(dtype): X = X.astype(dtype)
        else:  # is np.ndarray
            X = X.astype(dtype, copy=False)

        if X.dtype.names is None and len(X.shape) not in {0, 1, 2}:
            raise ValueError('X needs to be 2-dimensional, not '
                             '{}D'.format(len(X.shape)))

        # fix shapes
        self.X = X
        if len(self.X.shape) == 2:
            # TODO: int immutable, copy of references to ints in self.X.shape
            # only valid until accidental change
            self.n_smps, self.n_vars = self.X.shape
            # flatten to emulate numpys behavior upon slicing
            if self.n_smps == 1 and self.n_vars == 1:
                self.X = self.X[0, 0]
            elif self.n_smps == 1 or self.n_vars == 1:
                if sp.issparse(self.X): self.X = self.X.toarray()
                self.X = self.X.flatten()
        elif len(self.X.shape) == 1 and single_col:
            self.n_smps = self.X.shape[0]
            self.n_vars = 1
        elif len(self.X.shape) == 1:
            self.n_vars = self.X.shape[0]
            self.n_smps = 1
        else:
            self.n_vars = 1
            self.n_smps = 1

        # annotations
        if smp is None:
            self.smp = pd.DataFrame(index=RangeIndex(0, self.n_smps, name=None))
        elif 'smp_names' in smp:
            self.smp = pd.DataFrame(smp, index=smp['smp_names'],
                                    columns=[k for k in smp.keys() if k != 'smp_names'])
        elif 'row_names' in smp:
            self.smp = pd.DataFrame(smp, index=smp['row_names'],
                                    columns=[k for k in smp.keys() if k != 'row_names'])
        elif smp is not None and len(smp) > 0:
            self.smp = pd.DataFrame(smp)
        else:
            self.smp = pd.DataFrame(index=RangeIndex(0, self.n_smps, name=None))

        if var is None:
            self.var = pd.DataFrame(index=RangeIndex(0, self.n_vars, name=None))
        elif 'var_names' in var:
            self.var = pd.DataFrame(var, index=var['var_names'],
                                    columns=[k for k in var.keys() if k != 'var_names'])
        elif 'col_names' in var:
            self.var = pd.DataFrame(var, index=var['col_names'],
                                    columns=[k for k in var.keys() if k != 'col_names'])
        elif var is not None and len(var) > 0:
            self.var = pd.DataFrame(var)
        else:
            self.var = pd.DataFrame(index=RangeIndex(0, self.n_vars, name=None))

        # unstructured annotations
        self.uns = uns or {}

        # multi-dimensional array annotations
        if smpm is None: smpm = np.empty(self.n_smps, dtype=[])
        if varm is None: varm = np.empty(self.n_vars, dtype=[])
        self.smpm = BoundRecArr(smpm, self, 'smpm')
        self.varm = BoundRecArr(varm, self, 'varm')

        self._check_dimensions()

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
            self.smpm[key] = self._get_multicol_field_smp(key)
        for key in _keys_multicol_var:
            self.varm[key] = self._get_multicol_field_var(key)

    def __setattr__(self, key, value):
        if key in {'smp', 'var'}:
            if not isinstance(value, pd.DataFrame):
                raise ValueError('Can only assign pandas.DataFrame.')
        if key in {'uns'}:
            if not isinstance(value, dict):
                raise ValueError('Can only assign dictionary.')
        object.__setattr__(self, key, value)
        if key in {'smp', 'var'}:
            self._check_dimensions(key)

    def __contains__(self, key):
        raise AttributeError('AnnData has no attribute __contains__, don\'t check `in adata`.')

    def __repr__(self):
        return ('AnnData object with n_smps × n_vars= {} × {}\n'
                '    smp_keys = {}\n'
                '    var_keys = {}\n'
                '    uns_keys = {}\n'
                '    smpm_keys = {}\n'
                '    varm_keys = {}'
                .format(self.n_smps, self.n_vars,
                        self.smp_keys(), self.var_keys(),
                        sorted(list(self.uns.keys())),
                        self.smpm_keys(), self.varm_keys()))

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

    def smp_keys(self):
        """Return keys of sample annotation ``smp``."""
        return self.smp.keys().tolist()

    def var_keys(self):
        """Return keys of variable annotation ``var``."""
        return self.var.keys().tolist()

    def smpm_keys(self):
        """Return keys of sample annotation ``smpm``."""
        return list(self.smpm.keys())

    def varm_keys(self):
        """Return keys of variable annotation ``varm``."""
        return list(self.varm.keys())

    @property
    def add(self):
        warnings.warn('Use `.self.uns` instead of `.add`, `.add` will bre removed in the future.', UserWarning)
        return self.uns

    @property
    def smp_names(self):
        """Samples index."""
        return self.smp.index.values

    @smp_names.setter
    def smp_names(self, names):
        self.smp.index = names

    @property
    def var_names(self):
        """Variables index."""
        return self.var.index.values

    @var_names.setter
    def var_names(self, names):
        self.var.index = names

    def uns_keys(self):
        """Return keys of unstructured annotation."""
        return self.uns.keys()

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
        del self.X[smp, var]
        if var == slice(None):
            del self.smp.iloc[smp, :]
        if smp == slice(None):
            del self.var.iloc[var, :]

    def __getitem__(self, index):
        # Note: this cannot be made inplace
        # http://stackoverflow.com/questions/31916617/using-keyword-arguments-in-getitem-method-in-python
        smp, var = self._normalize_indices(index)
        X = self.X[smp, var]
        smp_new = self.smp.iloc[smp]
        smpm_new = self.smpm[smp]
        var_new = self.var.iloc[var]
        varm_new = self.varm[var]
        assert smp_new.shape[0] == X.shape[0], (smp, smp_new)
        assert var_new.shape[0] == X.shape[1], (var, var_new)
        uns_new = self.uns.copy()
        # slice sparse spatrices of n_smps × n_smps in self.uns
        if not (isinstance(smp, slice) and
                smp.start is None and smp.step is None and smp.stop is None):
            raised_warning = False
            for k, v in self.uns.items():  # TODO: make sure this really works as expected
                if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                    uns_new[k] = v.tocsc()[:, smp].tocsr()[smp, :]
                    if not raised_warning:
                        logg.warn('Slicing adjacency matrices can be dangerous. '
                                  'Consider recomputing the data graph.')
                        raised_warning = True
        adata = AnnData(X, smp_new, var_new, uns_new, smpm_new, varm_new)
        return adata

    def inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        self.X = self.X[:, index]
        self.n_vars = self.X.shape[1]
        self.var = self.var.iloc[index]
        self.varm = self.varm[index]
        return None

    def inplace_subset_smp(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        self.X = self.X[index, :]
        raised_warning = False
        # TODO: solve this in a better way, also for var
        for k, v in self.uns.items():
            if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                self.uns[k] = v.tocsc()[:, index].tocsr()[index, :]
                if not raised_warning:
                    logg.warn('Slicing adjacency matrices can be dangerous. '
                              'Consider recomputing the data graph.')
                    raised_warning = True
        self.n_smps = self.X.shape[0]
        self.smp = self.smp.iloc[index]
        self.smpm = self.smpm[index]
        return None

    def get_smp_array(self, k):
        """Get an array along the sample dimension by first looking up
        smp_keys and then var_names."""
        x = (self.smp[k] if k in self.smp_keys()
             else self[:, k].X if k in set(self.var_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in smp_keys or var_names.'
                             .format(k))
        return x

    def get_var_array(self, k):
        """Get an array along the variables dimension by first looking up
        var_keys and then smp_names."""
        x = (self.var[k] if k in self.var_keys()
             else self[k] if k in set(self.smp_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in var_keys or smp_names.'
                             .format(k))
        return x

    def __setitem__(self, index, val):
        smp, var = self._normalize_indices(index)
        self.X[smp, var] = val

    def __len__(self):
        return self.X.shape[0]

    def transpose(self):
        """Return a transposed view of the object.

        Sample axis (rows) and variable axis are interchanged. No additional memory.
        """
        if sp.isspmatrix_csr(self.X):
            return AnnData(self.X.T.tocsr(), self.var, self.smp, self.uns,
                           self.varm.flipped(), self.smpm.flipped())
        return AnnData(self.X.T, self.var, self.smp, self.uns,
                       self.varm.flipped(), self.smpm.flipped())

    T = property(transpose)

    def copy(self):
        """Full copy with memory allocated."""
        return AnnData(self.X.copy(), self.smp.copy(), self.var.copy(), self.uns.copy(),
                       self.smpm.copy(), self.varm.copy())

    def _check_dimensions(self, key=None):
        if key is None:
            key = {'smp', 'var', 'smpm', 'varm'}
        else:
            key = {key}
        if 'smp' in key and len(self.smp) != self.n_smps:
            raise ValueError('Sample annotation `smp` needs to have the same amount of '
                             'rows as data ({}), but has {} rows'
                             .format(self.n_smps, self.smp.shape[0]))
        if 'var' in key and len(self.var) != self.n_vars:
            raise ValueError('Variable annotation `var` needs to have the same amount of '
                             'columns as data  ({}), but has {} rows'
                             .format(self.n_vars, self.var.shape[0]))
        if 'smpm' in key and len(self.smpm) != self.n_smps:
            raise ValueError('Sample annotation `smpm` needs to have the same amount of '
                             'rows as data ({}), but has {} rows'
                             .format(self.n_smps, self.smp.shape[0]))
        if 'varm' in key and len(self.varm) != self.n_vars:
            raise ValueError('Variable annotation `varm` needs to have the same amount of '
                             'columns as data ({}), but has {} rows'
                             .format(self.n_vars, self.var.shape[0]))

    def to_dict_dataframes(self):
        d = {'X': pd.DataFrame(self.X, index=self.smp.index),
             'smp': self.smp,
             'var': self.var,
             'smpm': self.smpm.to_dataframe(),
             'varm': self.varm.to_dataframe()}
        d = merge_dicts(d, self.uns)
        return d

    def to_dict_fixed_width_arrays(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        smp_rec, uns_smp = df_to_records_fixed_width(self.smp)
        var_rec, uns_var = df_to_records_fixed_width(self.var)
        d = {'_X': self.X,
             '_smp': smp_rec,
             '_var': var_rec,
             '_smpm': self.smpm,
             '_varm': self.varm}
        d = merge_dicts(d, self.uns, uns_smp, uns_var)
        return d

    def from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.

        In particular, from a dict that has been written using
        ``AnnData.to_dict_fixed_width_arrays``.
        """
        # data matrix
        if '_X' in ddata:
            X = ddata['_X']
            del ddata['_X']
        elif 'X' in ddata:
            X = ddata['X']
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

        return X, smp, var, uns, smpm, varm


def merge_dicts(*ds):
    """Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = ds[0]
    for d in ds[1:]:
        result.update(d)
    return result
