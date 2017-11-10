# Author: Alex Wolf (http://falexwolf.de)
#         P. Angerer
"""Store an annotated data matrix.
"""
from collections import Mapping, Sequence
from collections import OrderedDict
from enum import Enum

import numpy as np
from numpy import ma
import pandas as pd
from scipy import sparse as sp
from scipy.sparse.sputils import IndexMixin

import logging as logg
import warnings


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sp.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


def _find_corresponding_multicol_key(key, keys_multicol):
    """Find the corresponding multicolumn key."""
    for mk in keys_multicol:
        if key.startswith(mk) and 'of' in key:
            return mk
    return None


def _gen_keys_from_multicol_key(key_multicol, n_keys):
    """Generates single-column keys from multicolumn key."""
    keys = [('{}{:03}of{:03}')
            .format(key_multicol, i+1, n_keys) for i in range(n_keys)]
    return keys


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
    smp :  pd.DataFrame, np.ndarray, dictionary or None, optional (default: None)
        A #samples × #sample_keys array containing sample names (`index`) and
        other sample annotation in the columns. A passed dict is converted to a
        record array.
    var : np.ndarray, dict or None, optional (default: None)
        The same as `smp`, but of shape `n_vars` × `n_var_keys` for annotation of
        variables.
    uns : dict or None, optional (default: None)
        Unstructured annotation for the whole dataset.
    dtype : simple `np.dtype`, optional (default: ``'float32'``)
        Convert data matrix to this type upon initialization.
    single_col : bool, optional (default: ``False``)
        Interpret one-dimensional input array as column.

    Notes
    -----
    A data matrix is flattened if either #samples (`n_smps`) or #variables
    (`n_vars`) is 1, so that Numpy's behavior is reproduced (``adata[:, 0].X ==
    adata.X[:, 0]``).
    """

    def __init__(self, data, smp=None, var=None, uns=None, dtype='float32', single_col=False):
        if isinstance(data, Mapping):
            if any((smp, var, uns)):
                raise ValueError('If `data` is a dict no further arguments must be provided.')
            X, smp, var, uns = self.from_dict(data)
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

        # multicolumn keys
        smp_keys_multicol = None
        if uns and '_smp_keys_multicol' in uns:
            self._keys_multicol_smp = uns['_smp_keys_multicol']
            del uns['_smp_keys_multicol']
        elif uns and 'smp_keys_multicol' in uns:
            self._keys_multicol_smp = uns['smp_keys_multicol']
            del uns['smp_keys_multicol']
        else:
            self._keys_multicol_smp = []
        var_keys_multicol = None
        if uns and '_var_keys_multicol' in uns:
            self._keys_multicol_var = uns['_var_keys_multicol']
            del uns['_var_keys_multicol']
        elif uns and 'var_keys_multicol' in uns:
            self._keys_multicol_var = uns['var_keys_multicol']
            del uns['var_keys_multicol']
        else:
            self._keys_multicol_var = []

        # annotations
        if smp is None:
            self.smp = pd.DataFrame(index=np.arange(self.n_smps))
        elif 'smp_names' in smp:
            self.smp = pd.DataFrame(smp, index=smp['smp_names'],
                                    columns=[k for k in smp.keys() if k != 'smp_names'])
        elif 'row_names' in smp:
            self.smp = pd.DataFrame(smp, index=smp['row_names'],
                                    columns=[k for k in smp.keys() if k != 'row_names'])
        else:
            self.smp = pd.DataFrame(smp)
        if var is None:
            self.var = pd.DataFrame(index=np.arange(self.n_vars))
        elif 'var_names' in var:
            self.var = pd.DataFrame(var, index=var['var_names'],
                                    columns=[k for k in var.keys() if k != 'var_names'])
        elif 'col_names' in var:
            self.var = pd.DataFrame(var, index=var['col_names'],
                                    columns=[k for k in var.keys() if k != 'col_names'])
        else:
            self.var = pd.DataFrame(var)
        self._check_dimensions()

        # unstructured annotations
        self.uns = uns or {}

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
        raise AttributeError("AnnData has no attribute __contains__, don't check `in adata`.")

    def __repr__(self):
        return ('AnnData object with n_smps × n_vars= {} × {}\n'
                '    smp_keys = {}\n'
                '    var_keys = {}\n'
                '    uns_keys = {}'
                .format(self.n_smps, self.n_vars,
                        self.smp_keys(), self.var_keys(), sorted(list(self.uns.keys()))))

    def set_multicol_field_smp(self, key_multicol, values):
        self._set_multicol_field('smp', key_multicol, values)

    def set_multicol_field_var(self, key_multicol, values):
        self._set_multicol_field('var', key_multicol, values)

    def _set_multicol_field(self, a, key_multicol, values):
        if key_multicol not in getattr(self, '_keys_multicol_' + a):
            getattr(self, '_keys_multicol_' + 'smp').append(key_multicol)
        keys = _gen_keys_from_multicol_key(key_multicol, len(values[0]))
        # remove all fields from the array that are not among the new keys
        keys_set = set(keys)
        for k in getattr(self, a).columns:
            if k.startswith(key_multicol) and k not in keys_set:
                getattr(self, a).drop(k, axis='column', inplace=True)
        for ik, k in enumerate(keys):
            getattr(self, a)[k] = values[:, ik]

    def get_multicol_field_smp(self, key_multicol):
        return self._get_multicol_field('smp', key_multicol)

    def get_multicol_field_var(self, key_multicol):
        return self._get_multicol_field('var', key_multicol)

    def _get_multicol_field(self, a, key_multicol):
        if key_multicol not in getattr(self, '_keys_multicol_' + a):
            raise KeyError('Only multicolumn keys are valid.')
        keys = []
        for k in getattr(self, a).columns:
            if k.startswith(key_multicol):
                keys.append(k)
        return getattr(self, a)[keys].values

    def smp_keys(self):
        """Return keys of sample annotation."""
        return self._keys('smp')

    def var_keys(self):
        """Return keys of variable annotation."""
        return self._keys('var')

    def _keys(self, a):
        keys = []
        for key in getattr(self, a).columns:
            mk = _find_corresponding_multicol_key(
                key, getattr(self, '_keys_multicol_' + a))
            if mk is None: keys.append(key)
            elif mk not in keys: keys.append(mk)
        return keys

    @property
    def add(self):
        warnings.warn('Use `.self.uns` instead of `.add` in the future.', DeprecationWarning)
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
        """Return keys of unstional unstructured annotation."""
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
        smp_ann = self.smp.iloc[smp]
        var_ann = self.var.iloc[var]
        assert smp_ann.shape[0] == X.shape[0], (smp, smp_ann)
        assert var_ann.shape[0] == X.shape[1], (var, var_ann)
        uns_ann = self.uns.copy()
        # slice sparse spatrices of n_smps × n_smps in self.uns
        if not (isinstance(smp, slice) and
                smp.start is None and smp.step is None and smp.stop is None):
            raised_warning = False
            for k, v in self.uns.items():  # TODO: make sure this really works as expected
                if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                    uns_ann[k] = v.tocsc()[:, smp].tocsr()[smp, :]
                    if not raised_warning:
                        logg.warn('Slicing adjacency matrices can be dangerous. '
                                  'Consider recomputing the data graph.')
                        raised_warning = True
        adata = AnnData(X, smp_ann, var_ann, uns_ann)
        return adata

    def inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        self.X = self.X[:, index]
        self.var = self.var.iloc[index]
        self.n_vars = self.X.shape[1]
        return None

    def inplace_subset_smp(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        self.X = self.X[index, :]
        self.smp = self.smp.iloc[index]
        raised_warning = False
        for k, v in self.uns.items():
            if isinstance(v, sp.spmatrix) and v.shape == (self.n_smps, self.n_smps):
                self.uns[k] = v.tocsc()[:, index].tocsr()[index, :]
                if not raised_warning:
                    logg.warn('Slicing adjacency matrices can be dangerous. '
                              'Consider recomputing the data graph.')
                    raised_warning = True
        self.n_smps = self.X.shape[0]
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
            return AnnData(self.X.T.tocsr(), self.var, self.smp, self.uns)
        return AnnData(self.X.T, self.var, self.smp, self.uns)

    T = property(transpose)

    def copy(self):
        """Full copy with memory allocated."""
        return AnnData(self.X.copy(), self.smp.copy(), self.var.copy(), self.uns.copy())

    def _check_dimensions(self, key=None):
        if key is None:
            key = {'smp', 'var'}
        else:
            key = {key}
        if 'smp' in key and len(self.smp) != self.n_smps:
            raise ValueError('Sample annotation needs to have the same amount of '
                             'rows as data has ({}), but has {} rows'
                             .format(self.n_smps, self.smp.shape[0]))
        if 'var' in key and len(self.var) != self.n_vars:
            raise ValueError('Feature annotation needs to have the same amount of '
                             'rows as data has columns ({}), but has {} rows'
                             .format(self.n_vars, self.var.shape[0]))

    def to_dict(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        from pandas.api.types import is_string_dtype, is_categorical
        for c in self.smp.columns:
            if is_string_dtype(self.smp[c]):
                cat = self.smp[c].astype('category')
                self.uns[c + '_categories'] = cat.cat.categories.values
                self.smp[c] = cat.cat.codes
            if is_categorical(self.smp[c]):
                self.uns[c + '_categories'] = self.smp[c].categories
        d = {'_X': self.X,
             '_smp': self.smp.to_records(),
             '_var': self.var.to_records()}
        for k, v in self.uns.items():
            d[k] = v
        d['_smp_keys_multicol'] = self._keys_multicol_smp
        d['_var_keys_multicol'] = self._keys_multicol_var
        return d

    def from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.
        """
        uns = dict(ddata.items())
        del ddata
        if '_X' in uns:
            X = uns['_X']
            del uns['_X']
        elif 'X' in uns:
            X = uns['X']
            del uns['X']
        if ('_smp' in uns and isinstance(uns['_smp'], (np.ndarray, pd.DataFrame))
            and '_var' in uns and isinstance(uns['_var'], (np.ndarray, pd.DataFrame))):
            smp = uns['_smp']
            del uns['_smp']
            var = uns['_var']
            del uns['_var']
        elif ('smp' in uns and isinstance(uns['smp'], (np.ndarray, pd.DataFrame))
              and 'var' in uns and isinstance(uns['var'], (np.ndarray, pd.DataFrame))):
            smp = uns['smp']
            del uns['smp']
            var = uns['var']
            del uns['var']
        else:
            smp, var = OrderedDict(), OrderedDict()
            if 'row_names' in uns:
                smp['smp_names'] = uns['row_names']
                del uns['row_names']
            elif 'smp_names' in uns:
                smp['smp_names'] = uns['smp_names']
                del uns['smp_names']

            if 'col_names' in uns:
                var['var_names'] = uns['col_names']
                del uns['col_names']
            elif 'var_names' in uns:
                var['var_names'] = uns['var_names']
                del uns['var_names']

            smp = merge_dicts(smp, uns.get('row', {}), uns.get('smp', {}))
            var = merge_dicts(var, uns.get('col', {}), uns.get('var', {}))
            for k in ['row', 'smp', 'col', 'var']:
                if k in uns:
                    del uns[k]

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

            if 'var_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='var_names')
            elif 'col_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='col_names')
            elif 'index' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='index')
            else:
                var = pd.DataFrame.from_records(var)

            k_to_delete = []
            for k in uns.keys():
                if k.endswith('_categories'):
                    k_stripped = k.replace('_categories', '')
                    if k_stripped in smp:
                        smp[k_stripped] = pd.Categorical.from_codes(
                            codes=smp[k_stripped].values,
                            categories=uns[k])
                    if k_stripped in var:
                        var[k_stripped] = pd.Categorical.from_codes(
                            codes=smp[k_stripped].values,
                            categories=uns[k])
                    k_to_delete.append(k)

            for k in k_to_delete:
                del uns[k]

        return X, smp, var, uns


def merge_dicts(*ds):
    """Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    Note
    ----
    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = ds[0]
    for d in ds[1:]:
        result.update(d)
    return result
